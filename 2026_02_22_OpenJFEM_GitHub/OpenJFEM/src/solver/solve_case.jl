# solve_case.jl — Subcase solver orchestrator

function solve_case(K, ndof, model, id_map, X, load_id, spc_id, node_R; max_elem_stiff=0.0, rbe3_map=Dict{Int,Vector{Tuple{Int,Float64}}}(), snorm_normals=Dict{Int,SVector{3,Float64}}())
    F_applied = zeros(ndof)
    if !isnothing(load_id)
        elem_map = Dict{Int, Any}()
        F_global_accum = zeros(ndof)
        resolve_loads(model, Int(load_id), 1.0, id_map, elem_map, X, F_global_accum)
        for i in 1:length(id_map)
            idx = (i-1)*6
            F_applied[idx+1:idx+3] = node_R[i]' * view(F_global_accum, idx+1:idx+3)
            F_applied[idx+4:idx+6] = node_R[i]' * view(F_global_accum, idx+4:idx+6)
        end
    end

    # RBE3 force redistribution
    if !isempty(rbe3_map)
        n_force_redist = 0
        for (dep_dof, pairs) in rbe3_map
            f_dep = F_applied[dep_dof]
            if abs(f_dep) > 1e-30
                for (ind_dof, coeff) in pairs
                    F_applied[ind_dof] += f_dep * coeff
                end
                F_applied[dep_dof] = 0.0
                n_force_redist += 1
            end
        end
        if n_force_redist > 0
            log_msg("[SOLVER] RBE3: Redistributed forces from $n_force_redist dependent DOFs")
        end
    end

    # Store spc_id in model for apply_bc_and_solve to access
    model["_spc_id"] = spc_id

    # Apply boundary conditions and solve
    u_global, fixed_dofs = apply_bc_and_solve(K, ndof, model, id_map, F_applied, node_R, rbe3_map, max_elem_stiff)

    # Clean up temporary key
    delete!(model, "_spc_id")

    R = K * u_global - F_applied

    results_json = Dict("displacements" => [], "spc_forces" => [], "forces" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => [], "conrod" => [], "celas1" => []), "stresses" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => [], "conrod" => [], "celas1" => []), "strains" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => [], "conrod" => [], "celas1" => []))

    u_out = zeros(ndof)

    sorted_nodes = sort(collect(keys(id_map)))
    for nid in sorted_nodes
        idx = id_map[nid]; base = (idx-1)*6
        u_loc = view(u_global, base+1:base+6)

        t_glob = node_R[idx] * u_loc[1:3]
        r_glob = node_R[idx] * u_loc[4:6]

        u_out[base+1:base+3] = t_glob
        u_out[base+4:base+6] = r_glob

        push!(results_json["displacements"], Dict("grid_id" => nid, "t1" => t_glob[1], "t2" => t_glob[2], "t3" => t_glob[3], "r1" => r_glob[1], "r2" => r_glob[2], "r3" => r_glob[3]))

        if (base+1) in fixed_dofs || (base+2) in fixed_dofs || (base+3) in fixed_dofs
             r_loc = view(R, base+1:base+6)
             r_reac_glob = vcat(node_R[idx] * r_loc[1:3], node_R[idx] * r_loc[4:6])
             push!(results_json["spc_forces"], Dict("grid_id" => nid, "t1" => r_reac_glob[1], "t2" => r_reac_glob[2], "t3" => r_reac_glob[3], "r1" => r_reac_glob[4], "r2" => r_reac_glob[5], "r3" => r_reac_glob[6]))
        end
    end

    stresses = Dict{Int, Float64}()

    # Stress recovery for all element types
    recover_shell_stresses!(model, id_map, X, node_R, u_global, snorm_normals, stresses, results_json)
    recover_bar_stresses!(model, id_map, X, node_R, u_global, stresses, results_json)
    recover_rod_stresses!(model, id_map, X, node_R, u_global, stresses, results_json)
    recover_spring_forces!(model, id_map, u_global, stresses, results_json)

    return u_out, stresses, results_json
end
