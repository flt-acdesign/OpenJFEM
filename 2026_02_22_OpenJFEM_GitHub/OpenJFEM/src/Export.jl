# ============================================================================
# Export.jl — Export functions for JFEM
#
# This file is included at the top level (not a module) and has access to
# WriteVTK and JSON packages from the parent scope.
#
# Contains:
#   sanitize!(d)                  — NaN/Inf cleaner for JSON export
#   build_jfem_element_tables     — builds sorted element tables for JFEM binary
#   build_jfem_constraint_tables  — builds CELAS/RBE2/RBE3 tables for JFEM v3
#   collect_spc_data              — collects SPC constraints for a subcase
#   collect_point_loads           — collects FORCE/MOMENT cards for a subcase
#   collect_jfem_subcase_data     — collects per-subcase JFEM binary data
#   export_vtk_subcase            — exports a single subcase to VTK format
#   export_json                   — exports aggregated results to JSON
#   export_jfem_binary            — exports mesh + results to JFEM binary format (v3)
#   export_card_inventory         — exports card inventory to JSON
# ============================================================================

function sanitize!(d)
    if d isa Dict
        for (k,v) in d; d[k] = sanitize!(v); end
    elseif d isa Vector
        for i in eachindex(d); d[i] = sanitize!(d[i]); end
    elseif d isa Float64
        if isnan(d) || isinf(d); return 0.0; end
    end
    return d
end

function build_jfem_element_tables(model, id_map)
    jfem_node_ids = sort(collect(keys(id_map)))
    pshells = model["PSHELLs"]
    pbarls  = model["PBARLs"]
    prods_m = model["PRODs"]

    # Helper to look up property by PID (handles both string and int keys)
    function find_prop(pdict, pid)
        p = get(pdict, string(pid), nothing)
        if p === nothing; p = get(pdict, pid, nothing); end
        return p
    end

    jfem_quads = Tuple{Int,Int,Vector{Int},Float32}[]   # eid, pid, nodes, thickness
    jfem_trias = Tuple{Int,Int,Vector{Int},Float32}[]
    for (id, el) in model["CSHELLs"]
        eid = parse(Int, id)
        if !haskey(el, "NODES"); continue; end
        nids = el["NODES"]
        if !all(n -> haskey(id_map, n), nids); continue; end
        pid = get(el, "PID", 0)
        prop = find_prop(pshells, pid)
        t = Float32(prop !== nothing ? get(prop, "T", 0.0) : 0.0)
        if length(nids) == 4
            push!(jfem_quads, (eid, pid, nids, t))
        elseif length(nids) == 3
            push!(jfem_trias, (eid, pid, nids, t))
        end
    end
    sort!(jfem_quads, by=x->x[1])
    sort!(jfem_trias, by=x->x[1])

    jfem_bars = Tuple{Int,Int,Int,Int,Float32}[]   # eid, pid, ga, gb, area
    for (id, bar) in model["CBARs"]
        eid = parse(Int, id)
        if !haskey(bar, "GA"); continue; end
        ga, gb = bar["GA"], bar["GB"]
        if !haskey(id_map, ga) || !haskey(id_map, gb); continue; end
        pid = get(bar, "PID", 0)
        prop = find_prop(pbarls, pid)
        a = Float32(prop !== nothing ? get(prop, "A", 0.0) : 0.0)
        push!(jfem_bars, (eid, pid, ga, gb, a))
    end
    sort!(jfem_bars, by=x->x[1])

    jfem_rods = Tuple{Int,Int,Int,Int,Float32}[]   # eid, pid, ga, gb, area
    for (id, rod) in model["CRODs"]
        eid = parse(Int, id)
        if !haskey(rod, "GA"); continue; end
        ga, gb = rod["GA"], rod["GB"]
        if !haskey(id_map, ga) || !haskey(id_map, gb); continue; end
        pid = get(rod, "PID", 0)
        prop = find_prop(prods_m, pid)
        a = Float32(prop !== nothing ? get(prop, "A", 0.0) : 0.0)
        push!(jfem_rods, (eid, pid, ga, gb, a))
    end
    sort!(jfem_rods, by=x->x[1])

    return jfem_node_ids, jfem_quads, jfem_trias, jfem_bars, jfem_rods
end

# --- v3 extension: constraint tables and per-subcase load data ---

# Coordinate transform for force/moment direction vectors (duplicated from solver/helpers.jl to avoid cross-module dependency)
function _export_coord_transform(model, cid, vec)
    if cid == 0; return vec; end
    if !haskey(model["CORDs"], string(cid)); return vec; end
    cord = model["CORDs"][string(cid)]
    R = hcat(cord["U"], cord["V"], cord["W"])
    return R * vec
end

function build_jfem_constraint_tables(model, id_map)
    # --- CELAS1 springs ---
    jfem_celas = Tuple{Int,Int,Int,Int,Int,Float32}[]   # eid, g1, c1, g2, c2, stiffness
    pelases = get(model, "PELASs", Dict())
    for (id, el) in get(model, "CELASs", Dict())
        eid = parse(Int, id)
        g1 = get(el, "G1", 0); c1 = get(el, "C1", 0)
        g2 = get(el, "G2", 0); c2 = get(el, "C2", 0)
        pid = get(el, "PID", 0)
        pelas = get(pelases, string(pid), nothing)
        if pelas === nothing; pelas = get(pelases, pid, nothing); end
        K_stiff = Float32(pelas !== nothing ? get(pelas, "K", 0.0) : 0.0)
        push!(jfem_celas, (eid, g1, c1, g2, c2, K_stiff))
    end
    sort!(jfem_celas, by=x->x[1])

    # --- RBE2 rigid body elements ---
    jfem_rbe2s = []   # (eid, gn, cm, slave_nids::Vector{Int})
    for (id, rbe) in get(model, "RBE2s", Dict())
        eid = parse(Int, id)
        gn = rbe["GN"]; cm = Int(rbe["CM"])
        slaves = Int.(rbe["GM"])
        push!(jfem_rbe2s, (eid=eid, gn=gn, cm=cm, slaves=slaves))
    end
    sort!(jfem_rbe2s, by=x->x.eid)

    # --- RBE3 interpolation elements ---
    jfem_rbe3s = []   # (eid, refgrid, refc, dep_grids::Vector{Int})
    for (id, rbe) in get(model, "RBE3s", Dict())
        eid = parse(Int, id)
        refgrid = rbe["REFGRID"]; refc = Int(rbe["REFC"])
        deps = Int.(rbe["DEP_GRIDS"])
        push!(jfem_rbe3s, (eid=eid, refgrid=refgrid, refc=refc, deps=deps))
    end
    sort!(jfem_rbe3s, by=x->x.eid)

    return jfem_celas, jfem_rbe2s, jfem_rbe3s
end

function collect_spc_data(model, spc_id)
    # Returns Dict{Int, Int} mapping nid → dof_mask (e.g., 123456)
    spc_nodes = Dict{Int, Set{Int}}()
    if isnothing(spc_id); return Dict{Int,Int}(); end

    # Resolve SPCADD
    sets = Set{Int}()
    sid = Int(spc_id)
    if haskey(model["SPCADDs"], sid)
        union!(sets, model["SPCADDs"][sid])
    else
        push!(sets, sid)
    end

    # Collect SPC1 entries matching the set
    for spc in model["SPC1s"]
        if Int(spc["SID"]) in sets
            for n in spc["NODES"]
                if !haskey(spc_nodes, n); spc_nodes[n] = Set{Int}(); end
                for ch in spc["C"]
                    if isdigit(ch); push!(spc_nodes[n], parse(Int, string(ch))); end
                end
            end
        end
    end

    # Convert to integer mask: Set{1,2,3} → 123
    result = Dict{Int,Int}()
    for (nid, dofs) in spc_nodes
        mask = 0
        for d in sort(collect(dofs))
            mask = mask * 10 + d
        end
        result[nid] = mask
    end
    return result
end

function _collect_point_loads_recursive(model, sid, scale, forces_acc, moments_acc)
    # Collect FORCE cards
    for frc in get(model, "FORCEs", [])
        if Int(frc["SID"]) == sid
            gid = frc["GID"]
            global_dir = _export_coord_transform(model, Int(frc["CID"]), frc["Dir"])
            fvec = global_dir .* (frc["Mag"] * scale)
            if !haskey(forces_acc, gid); forces_acc[gid] = zeros(3); end
            forces_acc[gid] .+= fvec
        end
    end

    # Collect MOMENT cards
    for mom in get(model, "MOMENTs", [])
        if Int(mom["SID"]) == sid
            gid = mom["GID"]
            global_dir = _export_coord_transform(model, Int(mom["CID"]), mom["Dir"])
            mvec = global_dir .* (mom["Mag"] * scale)
            if !haskey(moments_acc, gid); moments_acc[gid] = zeros(3); end
            moments_acc[gid] .+= mvec
        end
    end

    # Recurse through LOAD combos
    for c in get(model, "LOAD_COMBOS", [])
        if Int(c["SID"]) == sid
            for sub in c["COMPS"]
                _collect_point_loads_recursive(model, Int(sub["LID"]), scale * c["S"] * sub["S"], forces_acc, moments_acc)
            end
        end
    end
end

function collect_point_loads(model, load_id)
    forces_acc = Dict{Int, Vector{Float64}}()
    moments_acc = Dict{Int, Vector{Float64}}()
    if isnothing(load_id); return forces_acc, moments_acc; end
    _collect_point_loads_recursive(model, Int(load_id), 1.0, forces_acc, moments_acc)
    return forces_acc, moments_acc
end

function collect_jfem_subcase_data(u, sub_res, id_map, jfem_node_ids, jfem_quads, jfem_trias, jfem_bars, jfem_rods; model=nothing, spc_id=nothing, load_id=nothing)
    safe_f32(x) = (v = Float32(x); isnan(v) || isinf(v) ? Float32(0) : v)
    nNodes_jfem = length(jfem_node_ids)
    nQuads_jfem = length(jfem_quads)
    nTrias_jfem = length(jfem_trias)
    nBars_jfem  = length(jfem_bars)
    nRods_jfem  = length(jfem_rods)

    # Displacements: 6 per node in jfem_node_ids order
    disp_jfem = Vector{Float32}(undef, nNodes_jfem * 6)
    for (i, nid) in enumerate(jfem_node_ids)
        idx = id_map[nid]
        for k in 1:6
            disp_jfem[(i-1)*6+k] = safe_f32(u[(idx-1)*6+k])
        end
    end

    # Shell results: 7 per shell (fx, fy, fxy, mx, my, mxy, vonmises)
    nShells = nQuads_jfem + nTrias_jfem
    shell_jfem = zeros(Float32, nShells * 7)
    shell_force_map = Dict{Int,Any}()
    for f in sub_res["forces"]["quad4"]; shell_force_map[f["eid"]] = f; end
    for f in sub_res["forces"]["tria3"]; shell_force_map[f["eid"]] = f; end
    shell_stress_map = Dict{Int,Any}()
    for s in sub_res["stresses"]["quad4"]; shell_stress_map[s["eid"]] = s; end
    for s in sub_res["stresses"]["tria3"]; shell_stress_map[s["eid"]] = s; end

    for (i, (eid, _, _, _)) in enumerate(jfem_quads)
        base = (i-1) * 7
        f = get(shell_force_map, eid, nothing)
        s = get(shell_stress_map, eid, nothing)
        if f !== nothing
            shell_jfem[base+1] = safe_f32(f["fx"]);  shell_jfem[base+2] = safe_f32(f["fy"])
            shell_jfem[base+3] = safe_f32(f["fxy"]); shell_jfem[base+4] = safe_f32(f["mx"])
            shell_jfem[base+5] = safe_f32(f["my"]);  shell_jfem[base+6] = safe_f32(f["mxy"])
        end
        if s !== nothing
            vm = max(s["z1"]["von_mises"], s["z2"]["von_mises"])
            shell_jfem[base+7] = safe_f32(vm)
        end
    end
    for (i, (eid, _, _, _)) in enumerate(jfem_trias)
        base = (nQuads_jfem + i - 1) * 7
        f = get(shell_force_map, eid, nothing)
        s = get(shell_stress_map, eid, nothing)
        if f !== nothing
            shell_jfem[base+1] = safe_f32(f["fx"]);  shell_jfem[base+2] = safe_f32(f["fy"])
            shell_jfem[base+3] = safe_f32(f["fxy"]); shell_jfem[base+4] = safe_f32(f["mx"])
            shell_jfem[base+5] = safe_f32(f["my"]);  shell_jfem[base+6] = safe_f32(f["mxy"])
        end
        if s !== nothing
            vm = max(s["z1"]["von_mises"], s["z2"]["von_mises"])
            shell_jfem[base+7] = safe_f32(vm)
        end
    end

    # Bar results: 7 per bar (axial, shear_1, shear_2, torque, moment_a1, moment_a2, bar_vonmises)
    bar_jfem = zeros(Float32, nBars_jfem * 7)
    bar_force_map = Dict{Int,Any}()
    for f in sub_res["forces"]["cbar"]; bar_force_map[f["eid"]] = f; end
    bar_stress_map = Dict{Int,Any}()
    for s in sub_res["stresses"]["cbar"]; bar_stress_map[s["eid"]] = s; end

    for (i, (eid, _, _, _, _)) in enumerate(jfem_bars)
        base = (i-1) * 7
        f = get(bar_force_map, eid, nothing)
        s = get(bar_stress_map, eid, nothing)
        if f !== nothing
            bar_jfem[base+1] = safe_f32(f["axial"]);    bar_jfem[base+2] = safe_f32(f["shear_1"])
            bar_jfem[base+3] = safe_f32(f["shear_2"]);  bar_jfem[base+4] = safe_f32(f["torque"])
            bar_jfem[base+5] = safe_f32(f["moment_a1"]); bar_jfem[base+6] = safe_f32(f["moment_a2"])
        end
        if s !== nothing
            vm = abs(s["axial"])
            for pk in ["p1","p2","p3","p4"]
                vm = max(vm, abs(get(s["end_a"], pk, 0.0)), abs(get(s["end_b"], pk, 0.0)))
            end
            bar_jfem[base+7] = safe_f32(vm)
        end
    end

    # Rod results: 2 per rod (axial, torque)
    rod_jfem = zeros(Float32, nRods_jfem * 2)
    rod_force_map = Dict{Int,Any}()
    for f in sub_res["forces"]["crod"]; rod_force_map[f["eid"]] = f; end

    for (i, (eid, _, _, _, _)) in enumerate(jfem_rods)
        base = (i-1) * 2
        f = get(rod_force_map, eid, nothing)
        if f !== nothing
            rod_jfem[base+1] = safe_f32(f["axial"]); rod_jfem[base+2] = safe_f32(f["torque"])
        end
    end

    # --- v3: SPC, forces, moments per subcase ---
    spc_data = Tuple{Int32, UInt32}[]
    force_data = Tuple{Int32, Float32, Float32, Float32}[]
    moment_data = Tuple{Int32, Float32, Float32, Float32}[]

    if model !== nothing
        # Collect SPC constraints
        spc_map = collect_spc_data(model, spc_id)
        for (nid, mask) in sort(collect(spc_map), by=x->x[1])
            push!(spc_data, (Int32(nid), UInt32(mask)))
        end

        # Collect point forces and moments
        forces_dict, moments_dict = collect_point_loads(model, load_id)
        for (nid, fvec) in sort(collect(forces_dict), by=x->x[1])
            if norm(fvec) > 1e-30
                push!(force_data, (Int32(nid), safe_f32(fvec[1]), safe_f32(fvec[2]), safe_f32(fvec[3])))
            end
        end
        for (nid, mvec) in sort(collect(moments_dict), by=x->x[1])
            if norm(mvec) > 1e-30
                push!(moment_data, (Int32(nid), safe_f32(mvec[1]), safe_f32(mvec[2]), safe_f32(mvec[3])))
            end
        end
    end

    return (disp=disp_jfem, shell=shell_jfem, bar=bar_jfem, rod=rod_jfem,
            spc=spc_data, forces=force_data, moments=moment_data)
end

function export_vtk_subcase(filename, output_dir, sid, model, id_map, X, u, stresses)
    base_name = basename(filename)
    vtk_base = replace(base_name, ".bdf" => "") * "_Subcase_$sid"
    vtk_path = joinpath(output_dir, vtk_base)
    points = zeros(3, length(id_map))
    disp = zeros(3, length(id_map))
    for (nid, idx) in id_map
         points[:, idx] = X[idx, :]
         disp[:, idx] = u[(idx-1)*6+1:(idx-1)*6+3]
    end
    cells = MeshCell[]
    data_vonmises = Float64[]
    for (id, el) in model["CSHELLs"]
        if !haskey(el, "NODES"); continue; end
        eid = parse(Int, id); nids = [get(id_map, n, 0) for n in el["NODES"]]; if 0 in nids; continue; end
        if length(nids) == 3
            push!(cells, MeshCell(VTKCellTypes.VTK_TRIANGLE, nids)); push!(data_vonmises, get(stresses, eid, 0.0))
        elseif length(nids) == 4
            push!(cells, MeshCell(VTKCellTypes.VTK_QUAD, nids)); push!(data_vonmises, get(stresses, eid, 0.0))
        end
    end
    for (id, bar) in model["CBARs"]
         if !haskey(bar, "GA"); continue; end
         eid = parse(Int, id); nids = [get(id_map, bar["GA"], 0), get(id_map, bar["GB"], 0)]; if 0 in nids; continue; end
         push!(cells, MeshCell(VTKCellTypes.VTK_LINE, nids)); push!(data_vonmises, get(stresses, eid, 0.0))
    end
    for (id, rod) in model["CRODs"]
         if !haskey(rod, "GA"); continue; end
         eid = parse(Int, id); nids = [get(id_map, rod["GA"], 0), get(id_map, rod["GB"], 0)]; if 0 in nids; continue; end
         push!(cells, MeshCell(VTKCellTypes.VTK_LINE, nids)); push!(data_vonmises, get(stresses, eid, 0.0))
    end
    if !isempty(cells)
        vtk = vtk_grid(vtk_path, points, cells)
        vtk["Displacement", VTKPointData()] = disp
        vtk["VonMises_Stress", VTKCellData()] = data_vonmises
        vtk_save(vtk)
        println("  VTK saved: $vtk_path.vtu")
    end
end

function export_json(filename, output_dir, global_results)
    base_name = basename(filename)
    json_name = replace(base_name, ".bdf" => "") * ".JU.JSON"
    json_path = joinpath(output_dir, json_name)
    println("\n>>> Exporting AGGREGATED JSON: $json_path")
    sanitize!(global_results)
    open(json_path, "w") do f; JSON.print(f, global_results, 4); end
end

function export_jfem_binary(filename, output_dir, id_map, X, jfem_node_ids, jfem_quads, jfem_trias, jfem_bars, jfem_rods, jfem_subcases_data; jfem_celas=[], jfem_rbe2s=[], jfem_rbe3s=[])
    base_name = basename(filename)
    jfem_name = replace(base_name, ".bdf" => "") * ".jfem"
    jfem_path = joinpath(output_dir, jfem_name)
    nNodes_jfem = length(jfem_node_ids)
    nQuads_jfem = length(jfem_quads)
    nTrias_jfem = length(jfem_trias)
    nBars_jfem  = length(jfem_bars)
    nRods_jfem  = length(jfem_rods)
    nCelas_jfem = length(jfem_celas)
    nRBE2_jfem  = length(jfem_rbe2s)
    nRBE3_jfem  = length(jfem_rbe3s)
    println("\n>>> Exporting JFEM binary (v3): $jfem_path")
    open(jfem_path, "w") do io
        # Magic: 'JFEM'
        write(io, UInt8('J')); write(io, UInt8('F')); write(io, UInt8('E')); write(io, UInt8('M'))
        # Header (v3: extended with constraint counts)
        write(io, UInt32(3))                          # version 3
        write(io, UInt32(nNodes_jfem))
        write(io, UInt32(nQuads_jfem))
        write(io, UInt32(nTrias_jfem))
        write(io, UInt32(nBars_jfem))
        write(io, UInt32(nRods_jfem))
        write(io, UInt32(length(jfem_subcases_data)))  # nSubcases
        write(io, UInt32(nCelas_jfem))                 # v3: nCelas
        write(io, UInt32(nRBE2_jfem))                  # v3: nRBE2
        write(io, UInt32(nRBE3_jfem))                  # v3: nRBE3

        # Node table: nid(i32), x(f32), y(f32), z(f32)
        for nid in jfem_node_ids
            idx = id_map[nid]
            write(io, Int32(nid))
            write(io, Float32(X[idx, 1])); write(io, Float32(X[idx, 2])); write(io, Float32(X[idx, 3]))
        end

        # CQUAD4 table: eid(i32), pid(i32), g1-g4(i32), thickness(f32)
        for (eid, pid, nodes, t) in jfem_quads
            write(io, Int32(eid)); write(io, Int32(pid))
            for n in nodes; write(io, Int32(n)); end
            write(io, t)
        end

        # CTRIA3 table: eid(i32), pid(i32), g1-g3(i32), thickness(f32)
        for (eid, pid, nodes, t) in jfem_trias
            write(io, Int32(eid)); write(io, Int32(pid))
            for n in nodes; write(io, Int32(n)); end
            write(io, t)
        end

        # CBAR table: eid(i32), pid(i32), ga(i32), gb(i32), area(f32)
        for (eid, pid, ga, gb, a) in jfem_bars
            write(io, Int32(eid)); write(io, Int32(pid)); write(io, Int32(ga)); write(io, Int32(gb))
            write(io, a)
        end

        # CROD table: eid(i32), pid(i32), ga(i32), gb(i32), area(f32)
        for (eid, pid, ga, gb, a) in jfem_rods
            write(io, Int32(eid)); write(io, Int32(pid)); write(io, Int32(ga)); write(io, Int32(gb))
            write(io, a)
        end

        # v3: CELAS table: eid(i32), g1(i32), c1(i32), g2(i32), c2(i32), stiffness(f32), pad(f32)
        for (eid, g1, c1, g2, c2, K_stiff) in jfem_celas
            write(io, Int32(eid)); write(io, Int32(g1)); write(io, Int32(c1))
            write(io, Int32(g2)); write(io, Int32(c2)); write(io, K_stiff); write(io, Float32(0))
        end

        # v3: RBE2 table (variable-length): eid(i32), gn(i32), cm(i32), nSlaves(u32), [slave_nid(i32) × nSlaves]
        for rbe in jfem_rbe2s
            write(io, Int32(rbe.eid)); write(io, Int32(rbe.gn)); write(io, Int32(rbe.cm))
            write(io, UInt32(length(rbe.slaves)))
            for s in rbe.slaves; write(io, Int32(s)); end
        end

        # v3: RBE3 table (variable-length): eid(i32), refgrid(i32), refc(i32), nDep(u32), [dep_nid(i32) × nDep]
        for rbe in jfem_rbe3s
            write(io, Int32(rbe.eid)); write(io, Int32(rbe.refgrid)); write(io, Int32(rbe.refc))
            write(io, UInt32(length(rbe.deps)))
            for d in rbe.deps; write(io, Int32(d)); end
        end

        # Per-subcase data
        for sc in jfem_subcases_data
            write(io, UInt32(sc.sid))
            write(io, sc.disp)    # nNodes * 6 Float32
            write(io, sc.shell)   # (nQuads + nTrias) * 7 Float32
            write(io, sc.bar)     # nBars * 7 Float32
            write(io, sc.rod)     # nRods * 2 Float32

            # v3: SPC data
            write(io, UInt32(length(sc.spc)))
            for (nid, mask) in sc.spc
                write(io, nid); write(io, mask)
            end

            # v3: Applied forces
            write(io, UInt32(length(sc.forces)))
            for (nid, fx, fy, fz) in sc.forces
                write(io, nid); write(io, fx); write(io, fy); write(io, fz)
            end

            # v3: Applied moments
            write(io, UInt32(length(sc.moments)))
            for (nid, mx, my, mz) in sc.moments
                write(io, nid); write(io, mx); write(io, my); write(io, mz)
            end
        end
    end
    println("  JFEM v3: $(nNodes_jfem) nodes, $(nQuads_jfem)Q+$(nTrias_jfem)T shells, $(nBars_jfem) bars, $(nRods_jfem) rods, $(nCelas_jfem) springs, $(nRBE2_jfem) RBE2, $(nRBE3_jfem) RBE3, $(length(jfem_subcases_data)) subcases")
end

function export_card_inventory(cards, output_dir, filename)
    processed_card_types = Set([
        "GRID", "CORD2R", "CORD1R", "CORD2C", "CORD2S",
        "CTRIA3", "CQUAD4", "CBAR", "CROD", "CONROD", "CELAS1",
        "RBE2", "RBE3",
        "PSHELL", "PBARL", "PBAR", "PBAR*", "PROD", "PCOMP", "PELAS",
        "MAT1", "MAT2", "MAT8",
        "FORCE", "MOMENT", "PLOAD4", "PLOAD2", "PLOAD1", "GRAV",
        "SPC1", "SPC", "SPCADD", "MPC", "MPCADD", "LOAD",
        "CONM2",
        "PARAM"
    ])
    card_counts = Dict{String,Int}()
    unprocessed_cards = Dict{String,Int}()
    for (cname, clist) in cards
        card_counts[cname] = length(clist)
        if !(cname in processed_card_types)
            unprocessed_cards[cname] = length(clist)
        end
    end
    inv_json = Dict(
        "card_counts" => Dict(card_counts),
        "processed_card_types" => sort(collect(processed_card_types)),
        "unprocessed_cards" => Dict(unprocessed_cards)
    )
    inv_path = joinpath(output_dir, replace(basename(filename), ".bdf" => "") * ".CARDS.JSON")
    open(inv_path, "w") do f; JSON.print(f, inv_json, 4); end
    println(">>> Card inventory exported: $inv_path")
    if !isempty(unprocessed_cards)
        println("    WARNING: $(length(unprocessed_cards)) unprocessed card type(s):")
        for (cname, cnt) in sort(collect(unprocessed_cards), by=x->x[1])
            println("      $cname: $cnt")
        end
    end
end
