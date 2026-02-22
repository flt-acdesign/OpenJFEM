# constraints.jl — RBE2, RBE3, MPC constraint assembly and DOF elimination

# Extracted from assemble_stiffness: processes all constraint elements and
# redistributes stiffness triplets for dependent DOFs.
# Returns: (rbe3_map, I_idx, J_idx, V_val)  — rbe3_map is the merged constraint map,
# and triplet arrays may be replaced if constraints exist.
function assemble_constraints(model, id_map, node_coords, I_idx, J_idx, V_val)
    rbe2s = get(model, "RBE2s", Dict())

    # --- RBE2 (rigid body element - MPC constraint via DOF elimination) ---
    rbe2_map = Dict{Int, Vector{Tuple{Int, Float64}}}()
    n_rbe2 = 0
    n_rbe2_dep = 0
    for (id, rbe) in rbe2s
        gn = rbe["GN"]   # master node
        if !haskey(id_map, gn); continue; end
        i_master = id_map[gn]
        cm_digits = [parse(Int, string(ch)) for ch in string(rbe["CM"]) if isdigit(ch)]
        p_m = SVector{3}(node_coords[i_master,1], node_coords[i_master,2], node_coords[i_master,3])

        for gs in rbe["GM"]  # slave nodes
            if !haskey(id_map, gs); continue; end
            i_slave = id_map[gs]
            p_s = SVector{3}(node_coords[i_slave,1], node_coords[i_slave,2], node_coords[i_slave,3])
            dx, dy, dz = p_s[1]-p_m[1], p_s[2]-p_m[2], p_s[3]-p_m[3]

            for c in cm_digits
                slave_dof = (i_slave-1)*6 + c
                pairs = Tuple{Int,Float64}[]

                if c <= 3  # translation DOF
                    push!(pairs, ((i_master-1)*6 + c, 1.0))
                    if c == 1
                        if abs(dz) > 1e-12; push!(pairs, ((i_master-1)*6 + 5, dz)); end
                        if abs(dy) > 1e-12; push!(pairs, ((i_master-1)*6 + 6, -dy)); end
                    elseif c == 2
                        if abs(dx) > 1e-12; push!(pairs, ((i_master-1)*6 + 6, dx)); end
                        if abs(dz) > 1e-12; push!(pairs, ((i_master-1)*6 + 4, -dz)); end
                    else
                        if abs(dy) > 1e-12; push!(pairs, ((i_master-1)*6 + 4, dy)); end
                        if abs(dx) > 1e-12; push!(pairs, ((i_master-1)*6 + 5, -dx)); end
                    end
                else  # rotation DOF (4,5,6)
                    push!(pairs, ((i_master-1)*6 + c, 1.0))
                end

                rbe2_map[slave_dof] = pairs
                n_rbe2_dep += 1
            end
            n_rbe2 += 1
        end
    end
    if n_rbe2 > 0
        log_msg("[SOLVER] RBE2: $(length(rbe2s)) elements, $n_rbe2 master-slave pairs, $n_rbe2_dep dependent DOFs (MPC elimination)")
    end

    # --- RBE3 (weighted interpolation constraint via DOF elimination) ---
    rbe3s = get(model, "RBE3s", Dict())
    rbe3_map = Dict{Int, Vector{Tuple{Int, Float64}}}()
    n_rbe3 = 0
    for (id, rbe) in rbe3s
        ref_gid = rbe["REFGRID"]
        ref_idx = get(id_map, ref_gid, 0)
        if ref_idx == 0; continue; end

        refc_digits = sort!([parse(Int, string(ch)) for ch in string(Int(rbe["REFC"])) if isdigit(ch)])
        comps_digits = sort!([parse(Int, string(ch)) for ch in string(Int(rbe["COMPS"])) if isdigit(ch)])
        if isempty(refc_digits) || isempty(comps_digits); continue; end

        dep_grids = rbe["DEP_GRIDS"]
        wt = Float64(get(rbe, "WT", 1.0))
        p_ref = SVector{3}(node_coords[ref_idx,1], node_coords[ref_idx,2], node_coords[ref_idx,3])

        ind_list = Tuple{Int,SVector{3,Float64}}[]
        for dg in dep_grids
            di = get(id_map, dg, 0)
            if di > 0
                p_i = SVector{3}(node_coords[di,1], node_coords[di,2], node_coords[di,3])
                push!(ind_list, (di, p_i))
            end
        end
        if isempty(ind_list); continue; end

        n_ref = length(refc_digits)
        n_comp = length(comps_digits)

        A = zeros(n_ref, n_ref)
        grid_Ri = Tuple{Int, Matrix{Float64}}[]

        for (di, p_i) in ind_list
            dx = p_i[1] - p_ref[1]
            dy = p_i[2] - p_ref[2]
            dz = p_i[3] - p_ref[3]

            R_i = zeros(n_comp, n_ref)
            for (jj, cdof) in enumerate(comps_digits)
                for (cc, rdof) in enumerate(refc_digits)
                    R_i[jj, cc] = _rbe3_rb_coeff(cdof, rdof, dx, dy, dz)
                end
            end

            A .+= wt .* (R_i' * R_i)
            push!(grid_Ri, (di, R_i))
        end

        A_inv = pinv(A, rtol=1e-10)

        for (di, R_i) in grid_Ri
            C_i = A_inv * (wt .* R_i')

            for (cc, rdof) in enumerate(refc_digits)
                ref_dof = (ref_idx - 1) * 6 + rdof
                if !haskey(rbe3_map, ref_dof)
                    rbe3_map[ref_dof] = Tuple{Int,Float64}[]
                end
                for (jj, cdof) in enumerate(comps_digits)
                    coeff = C_i[cc, jj]
                    if abs(coeff) > 1e-15
                        ind_dof = (di - 1) * 6 + cdof
                        push!(rbe3_map[ref_dof], (ind_dof, coeff))
                    end
                end
            end
        end
        n_rbe3 += 1
    end

    # --- Explicit MPC constraints ---
    mpc_cards = get(model, "MPCs", [])
    mpc_map = Dict{Int, Vector{Tuple{Int, Float64}}}()
    n_mpc_explicit = 0
    for mpc in mpc_cards
        terms = mpc["TERMS"]
        if length(terms) < 2; continue; end
        dep_g = terms[1]["G"]; dep_c = terms[1]["C"]; dep_a = terms[1]["A"]
        if !haskey(id_map, dep_g) || dep_c < 1 || dep_c > 6; continue; end
        if abs(dep_a) < 1e-30; continue; end
        dep_dof = (id_map[dep_g]-1)*6 + dep_c
        pairs = Tuple{Int,Float64}[]
        for i in 2:length(terms)
            t = terms[i]
            if !haskey(id_map, t["G"]) || t["C"] < 1 || t["C"] > 6; continue; end
            ind_dof = (id_map[t["G"]]-1)*6 + t["C"]
            coeff = -t["A"] / dep_a
            push!(pairs, (ind_dof, coeff))
        end
        if !isempty(pairs)
            mpc_map[dep_dof] = pairs
            n_mpc_explicit += 1
        end
    end
    if n_mpc_explicit > 0
        log_msg("[SOLVER] MPC: $n_mpc_explicit explicit MPC constraints")
    end

    # Merge RBE2, RBE3, and explicit MPC constraint maps
    merge!(rbe3_map, rbe2_map, mpc_map)
    n_mpc_total = length(rbe3_map)

    if n_rbe3 > 0
        log_msg("[SOLVER] RBE3: $n_rbe3 elements, $(length(rbe3_map) - length(rbe2_map) - length(mpc_map)) dependent DOFs")
    end

    if n_mpc_total > 0
        log_msg("[SOLVER] MPC elimination: $n_mpc_total total dependent DOFs (RBE2: $(length(rbe2_map)), RBE3: $(length(rbe3_map) - length(rbe2_map) - length(mpc_map)), MPC: $(length(mpc_map)))")
        # Process triplets: redistribute entries involving dependent DOFs
        n_orig = length(I_idx)
        new_I = Int[]; new_J = Int[]; new_V = Float64[]
        sizehint!(new_I, n_orig + n_orig ÷ 10)
        sizehint!(new_J, n_orig + n_orig ÷ 10)
        sizehint!(new_V, n_orig + n_orig ÷ 10)

        for k in 1:n_orig
            i = I_idx[k]; j = J_idx[k]; v = V_val[k]
            i_dep = haskey(rbe3_map, i)
            j_dep = haskey(rbe3_map, j)

            if !i_dep && !j_dep
                push!(new_I, i); push!(new_J, j); push!(new_V, v)
            elseif i_dep && !j_dep
                for (ind_dof, coeff) in rbe3_map[i]
                    push!(new_I, ind_dof); push!(new_J, j); push!(new_V, v * coeff)
                end
            elseif !i_dep && j_dep
                for (ind_dof, coeff) in rbe3_map[j]
                    push!(new_I, i); push!(new_J, ind_dof); push!(new_V, v * coeff)
                end
            else
                for (ind_i, ci) in rbe3_map[i]
                    for (ind_j, cj) in rbe3_map[j]
                        push!(new_I, ind_i); push!(new_J, ind_j); push!(new_V, v * ci * cj)
                    end
                end
            end
        end
        I_idx = new_I; J_idx = new_J; V_val = new_V
        log_msg("[SOLVER] MPC: Triplets redistributed: $n_orig → $(length(I_idx))")
    end

    return rbe3_map, I_idx, J_idx, V_val
end
