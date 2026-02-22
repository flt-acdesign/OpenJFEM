# loads.jl — Load resolution (FORCE, MOMENT, PLOAD4, GRAV, PLOAD1, LOAD combos)

function resolve_loads(model, sid, scale, id_map, elem_map, node_coords, F_acc)
    raw_forces = Dict{Int, Vector{Float64}}()
    add_force = (gid, vec) -> begin
        if !haskey(raw_forces, gid); raw_forces[gid] = zeros(6); end
        raw_forces[gid] .+= vec
    end

    for frc in model["FORCEs"]; if Int(frc["SID"]) == sid
        global_dir = get_coord_transform(model, Int(frc["CID"]), frc["Dir"])
        add_force(frc["GID"], zeros(6)); raw_forces[frc["GID"]][1:3] .+= global_dir * frc["Mag"] * scale
    end; end

    for mom in model["MOMENTs"]; if Int(mom["SID"]) == sid
        global_dir = get_coord_transform(model, Int(mom["CID"]), mom["Dir"])
        add_force(mom["GID"], zeros(6)); raw_forces[mom["GID"]][4:6] .+= global_dir * mom["Mag"] * scale
    end; end

    for pload in model["PLOAD4s"]; if Int(pload["SID"]) == sid
        eid = pload["EID"]
        if haskey(model["CSHELLs"], string(eid))
            el_def = model["CSHELLs"][string(eid)]
            nids = [get(id_map, n, 0) for n in el_def["NODES"]]
            if !any(x->x==0, nids)
                Xc = node_coords[nids, :]
                v1 = Xc[2,:] - Xc[1,:]; v2 = Xc[3,:] - Xc[1,:]
                normal_vec = (length(nids) == 4) ? cross(Xc[3,:] - Xc[1,:], Xc[4,:] - Xc[2,:]) : cross(v1, v2)
                tf = 0.5 * norm(normal_vec) * pload["P"] * scale
                f_node = normalize(normal_vec) * (tf / length(nids))
                for idx in nids; dof = (idx-1)*6; F_acc[dof+1:dof+3] .+= f_node; end
            end
        end
    end; end

    # --- GRAV (gravity/acceleration body forces) ---
    for grav in get(model, "GRAVs", [])
        if Int(grav["SID"]) == sid
            dir_raw = grav["N"]
            n_dir = norm(dir_raw)
            grav_dir = n_dir > 1e-30 ? dir_raw ./ n_dir : [0.0, 0.0, 0.0]
            accel = grav["A"] * scale
            grav_vec = grav_dir .* accel

            # CONM2 concentrated masses
            for (_, cm) in get(model, "CONM2s", Dict())
                gid = cm["GID"]
                if haskey(id_map, gid)
                    f_mass = cm["M"] .* grav_vec
                    idx = id_map[gid]; dof = (idx-1)*6
                    F_acc[dof+1:dof+3] .+= f_mass
                end
            end

            # Shell element mass: rho * t * area, distributed to element nodes
            mats_m = model["MATs"]
            for (_, el) in get(model, "CSHELLs", Dict())
                if !haskey(el, "NODES"); continue; end
                nids = el["NODES"]; nn = length(nids)
                pid = string(get(el, "PID", 0))
                prop = get(model["PSHELLs"], pid, nothing)
                if prop === nothing; continue; end
                mid = string(get(prop, "MID", 0))
                mat = get(mats_m, mid, nothing)
                if mat === nothing; continue; end
                rho = get(mat, "RHO", 0.0)
                if rho <= 0; continue; end
                idxs = [get(id_map, n, 0) for n in nids]
                if any(x->x==0, idxs); continue; end
                Xc = node_coords[idxs, :]
                if nn == 3
                    v1 = Xc[2,:] - Xc[1,:]; v2 = Xc[3,:] - Xc[1,:]
                    area = 0.5 * norm(cross(v1, v2))
                else
                    d13 = Xc[3,:] - Xc[1,:]; d24 = Xc[4,:] - Xc[2,:]
                    area = 0.5 * norm(cross(d13, d24))
                end
                total_mass = rho * get(prop, "T", 0.0) * area
                f_per_node = (total_mass / nn) .* grav_vec
                for nid in nids
                    if haskey(id_map, nid)
                        dof = (id_map[nid]-1)*6
                        F_acc[dof+1:dof+3] .+= f_per_node
                    end
                end
            end

            # Bar element mass: rho * A * L, distributed to 2 nodes
            for (_, bar) in get(model, "CBARs", Dict())
                pid = string(get(bar, "PID", 0))
                prop = get(model["PBARLs"], pid, nothing)
                if prop === nothing; continue; end
                mid = string(get(prop, "MID", 0))
                mat = get(mats_m, mid, nothing)
                if mat === nothing; continue; end
                rho = get(mat, "RHO", 0.0)
                if rho <= 0; continue; end
                ga, gb = bar["GA"], bar["GB"]
                if !haskey(id_map, ga) || !haskey(id_map, gb); continue; end
                i1, i2 = id_map[ga], id_map[gb]
                L = norm(node_coords[i2,:] - node_coords[i1,:])
                total_mass = rho * get(prop, "A", 0.0) * L
                f_per_node = (total_mass / 2) .* grav_vec
                for (nid, idx) in [(ga, i1), (gb, i2)]
                    dof = (idx-1)*6
                    F_acc[dof+1:dof+3] .+= f_per_node
                end
            end

            # Rod element mass: rho * A * L (CROD and CONROD)
            for rodset in [get(model, "CRODs", Dict()), get(model, "CONRODs", Dict())]
                for (_, rod) in rodset
                    local A_rod, mid_rod
                    if haskey(rod, "MID")  # CONROD
                        A_rod = get(rod, "A", 0.0)
                        mid_rod = string(rod["MID"])
                    else  # CROD
                        pid = string(get(rod, "PID", 0))
                        prop = get(get(model, "PRODs", Dict()), pid, nothing)
                        if prop === nothing; continue; end
                        A_rod = get(prop, "A", 0.0)
                        mid_rod = string(get(prop, "MID", 0))
                    end
                    mat = get(mats_m, mid_rod, nothing)
                    if mat === nothing; continue; end
                    rho = get(mat, "RHO", 0.0)
                    if rho <= 0; continue; end
                    ga, gb = rod["GA"], rod["GB"]
                    if !haskey(id_map, ga) || !haskey(id_map, gb); continue; end
                    i1, i2 = id_map[ga], id_map[gb]
                    L = norm(node_coords[i2,:] - node_coords[i1,:])
                    total_mass = rho * A_rod * L
                    f_per_node = (total_mass / 2) .* grav_vec
                    for idx in [i1, i2]
                        dof = (idx-1)*6
                        F_acc[dof+1:dof+3] .+= f_per_node
                    end
                end
            end
        end
    end

    # --- PLOAD1 (distributed load on bar elements) ---
    for pload in get(model, "PLOAD1s", [])
        if Int(pload["SID"]) == sid
            eid = pload["EID"]
            bar = nothing
            for (bid, b) in get(model, "CBARs", Dict())
                if parse(Int, bid) == eid; bar = b; break; end
            end
            if bar === nothing; continue; end
            ga, gb = bar["GA"], bar["GB"]
            if !haskey(id_map, ga) || !haskey(id_map, gb); continue; end
            i1, i2 = id_map[ga], id_map[gb]
            p1g = SVector{3}(node_coords[i1,1], node_coords[i1,2], node_coords[i1,3])
            p2g = SVector{3}(node_coords[i2,1], node_coords[i2,2], node_coords[i2,3])
            L = norm(p2g - p1g)
            if L < 1e-9; continue; end
            vx = normalize(p2g - p1g)
            vref = resolve_bar_vref(bar, p1g, id_map, node_coords)
            vz = normalize(cross(vx, vref))
            vy = cross(vz, vx)

            # PLOAD1 TYPE: 1=FX, 2=FY, 3=FZ, 4=MX, 5=MY, 6=MZ (element local)
            ltype = pload["LOAD_TYPE"]
            x1_frac = pload["X1"]; p1_val = pload["P1"] * scale
            x2_frac = pload["X2"]; p2_val = pload["P2"] * scale
            seg_L = (x2_frac - x1_frac) * L
            if seg_L < 1e-12; continue; end

            # Equivalent nodal forces for linear distribution p(x) = p1 + (p2-p1)*x/L_seg
            f_a_mag = seg_L * (2*p1_val + p2_val) / 6
            f_b_mag = seg_L * (p1_val + 2*p2_val) / 6

            if ltype >= 1 && ltype <= 3
                # Force load in element local direction
                local_dir = ltype == 1 ? vx : (ltype == 2 ? vy : vz)
                # Transverse moments (only for forces perpendicular to bar axis)
                if ltype >= 2  # transverse
                    m_a_mag = seg_L^2 * (3*p1_val + 2*p2_val) / 60
                    m_b_mag = -seg_L^2 * (2*p1_val + 3*p2_val) / 60
                    if ltype == 2  # FY → moment about Z (vz)
                        F_acc[(i1-1)*6+1:(i1-1)*6+3] .+= f_a_mag .* local_dir
                        F_acc[(i2-1)*6+1:(i2-1)*6+3] .+= f_b_mag .* local_dir
                        F_acc[(i1-1)*6+4:(i1-1)*6+6] .+= m_a_mag .* vz
                        F_acc[(i2-1)*6+4:(i2-1)*6+6] .+= m_b_mag .* vz
                    else  # FZ → moment about Y (vy, negated sign)
                        F_acc[(i1-1)*6+1:(i1-1)*6+3] .+= f_a_mag .* local_dir
                        F_acc[(i2-1)*6+1:(i2-1)*6+3] .+= f_b_mag .* local_dir
                        F_acc[(i1-1)*6+4:(i1-1)*6+6] .+= (-m_a_mag) .* vy
                        F_acc[(i2-1)*6+4:(i2-1)*6+6] .+= (-m_b_mag) .* vy
                    end
                else  # axial (ltype==1): no moments
                    F_acc[(i1-1)*6+1:(i1-1)*6+3] .+= f_a_mag .* local_dir
                    F_acc[(i2-1)*6+1:(i2-1)*6+3] .+= f_b_mag .* local_dir
                end
            elseif ltype >= 4 && ltype <= 6
                # Moment load in element local direction
                local_dir = ltype == 4 ? vx : (ltype == 5 ? vy : vz)
                F_acc[(i1-1)*6+4:(i1-1)*6+6] .+= f_a_mag .* local_dir
                F_acc[(i2-1)*6+4:(i2-1)*6+6] .+= f_b_mag .* local_dir
            end
        end
    end

    # Apply all point forces directly to their grids.
    for (gid, f_vec) in raw_forces
        if haskey(id_map, gid); idx = id_map[gid]; dof = (idx-1)*6; F_acc[dof+1:dof+6] .+= f_vec; end
    end

    for c in model["LOAD_COMBOS"]
        if Int(c["SID"]) == sid
            for sub in c["COMPS"]; resolve_loads(model, Int(sub["LID"]), scale * c["S"] * sub["S"], id_map, elem_map, node_coords, F_acc); end
        end
    end
end

# RBE3 rigid body coefficient
function _rbe3_rb_coeff(comp_dof::Int, ref_dof::Int, dx::Float64, dy::Float64, dz::Float64)
    if comp_dof == ref_dof; return 1.0; end
    if comp_dof <= 3 && ref_dof >= 4
        if comp_dof == 1
            if ref_dof == 5; return dz; end
            if ref_dof == 6; return -dy; end
        elseif comp_dof == 2
            if ref_dof == 4; return -dz; end
            if ref_dof == 6; return dx; end
        elseif comp_dof == 3
            if ref_dof == 4; return dy; end
            if ref_dof == 5; return -dx; end
        end
    end
    return 0.0
end
