# boundary_conditions.jl — SPC application, AUTOSPC, and linear solve

function apply_bc_and_solve(K, ndof, model, id_map, F_applied, node_R, rbe3_map, max_elem_stiff)
    log_msg("[SOLVER] Processing Boundary Conditions...")

    fixed_dofs = Set{Int}()

    # Fix RBE3 dependent DOFs
    for dep_dof in keys(rbe3_map)
        push!(fixed_dofs, dep_dof)
    end
    if !isempty(rbe3_map)
        log_msg("[SOLVER] RBE3: Fixed $(length(rbe3_map)) dependent DOFs")
    end

    sets = Set{Int}()
    spc_id = get(model, "_spc_id", nothing)
    if !isnothing(spc_id)
         sid = Int(spc_id)
         if haskey(model["SPCADDs"], sid); union!(sets, model["SPCADDs"][sid]); else; push!(sets, sid); end
    end
    for spc in model["SPC1s"]
        if Int(spc["SID"]) in sets
            for n in spc["NODES"]
                idx = get(id_map, n, 0)
                if idx > 0
                    for c in spc["C"]
                        push!(fixed_dofs, (idx-1)*6 + parse(Int, c))
                    end
                end
            end
        end
    end

    # AUTOSPC Phase 1: absolute threshold
    autospc_threshold = 1e-6
    n_autospc = 0
    n_autospc_trans = 0
    n_autospc_rot = 0
    for i in 1:ndof
        if !(i in fixed_dofs) && abs(K[i,i]) < autospc_threshold
            push!(fixed_dofs, i)
            n_autospc += 1
            dof_local = mod(i-1, 6) + 1
            if dof_local <= 3; n_autospc_trans += 1; else; n_autospc_rot += 1; end
        end
    end

    # AUTOSPC Phase 2: per-grid relative threshold
    autospc_rel_threshold = 1e-4
    n_grids = ndof ÷ 6
    for igrid in 1:n_grids
        base = (igrid-1)*6
        trans_free = [base+k for k in 1:3 if !(base+k in fixed_dofs)]
        if length(trans_free) < 2; continue; end
        trans_diag = [K[d,d] for d in trans_free]
        max_td = maximum(trans_diag)
        if max_td < 1e-12; continue; end
        for (j, d) in enumerate(trans_free)
            if trans_diag[j] / max_td < autospc_rel_threshold
                push!(fixed_dofs, d)
                n_autospc += 1
                n_autospc_trans += 1
            end
        end
    end
    log_msg("[SOLVER] AUTOSPC: Constrained $n_autospc DOFs ($n_autospc_trans trans + $n_autospc_rot rot, threshold=$autospc_threshold/rel=$autospc_rel_threshold)")

    F_norm = norm(F_applied)
    F_max = maximum(abs.(F_applied))
    n_nonzero = count(x -> abs(x) > 1e-10, F_applied)
    log_msg("[SOLVER] Force vector: |F|=$(F_norm), max=$(F_max), nonzero DOFs=$n_nonzero")
    log_msg("[SOLVER] Slicing Matrix (Reducing System)...")
    all_dofs = 1:ndof
    free_dofs = sort(collect(setdiff(all_dofs, fixed_dofs)))
    log_msg("[SOLVER] Fixed DOFs: $(length(fixed_dofs)), Free DOFs: $(length(free_dofs))")

    K_ff = K[free_dofs, free_dofs]
    F_ff = F_applied[free_dofs]

    n_free = length(free_dofs)
    if n_free <= 500000
        log_msg("[SOLVER] Using Direct Solver (Cholesky) for $n_free DOFs...")
        u_ff = try
            F_chol = cholesky(Symmetric(K_ff))

            L_sparse = sparse(F_chol.L)
            L_diag = abs.(diag(L_sparse))
            K_diag = [abs(K_ff[i,i]) for i in 1:n_free]
            pivot_ratios = zeros(n_free)
            for i in 1:n_free
                if K_diag[i] > 1e-30
                    pivot_ratios[i] = L_diag[i]^2 / K_diag[i]
                else
                    pivot_ratios[i] = 1.0
                end
            end
            perm = F_chol.p
            sing_threshold = 1e-7
            singular_local = findall(pivot_ratios .< sing_threshold)
            n_sing = length(singular_local)
            if n_sing > 0
                singular_original = perm[singular_local]
                singular_global = free_dofs[singular_original]
                n_sing_trans = count(d -> mod(d-1,6)+1 <= 3, singular_global)
                n_sing_rot = n_sing - n_sing_trans
                log_msg("[SOLVER] Post-factorization singularity: $n_sing DOFs ($n_sing_trans trans + $n_sing_rot rot, threshold=$sing_threshold)")

                for d in singular_global; push!(fixed_dofs, d); end
                free_dofs = sort(collect(setdiff(1:ndof, fixed_dofs)))
                K_ff = K[free_dofs, free_dofs]
                F_ff = F_applied[free_dofs]
                n_free = length(free_dofs)
                log_msg("[SOLVER] Re-solving with $(length(fixed_dofs)) fixed, $n_free free DOFs")
                F_chol = cholesky(Symmetric(K_ff))
            end

            F_chol \ F_ff
        catch e
            shift_val = max_elem_stiff * 1e-12
            log_msg("[SOLVER] Cholesky failed: $e. Running factorization AUTOSPC (shift=$shift_val)...")

            F_chol_probe = cholesky(Symmetric(K_ff); shift=shift_val)
            L_sparse = sparse(F_chol_probe.L)
            L_diag = abs.(diag(L_sparse))
            L_median = median(L_diag)
            pivot_threshold = sqrt(shift_val) * 10.0
            small_pivot_mask = L_diag .< pivot_threshold
            n_mechanism = count(small_pivot_mask)

            if n_mechanism > 0
                perm = F_chol_probe.p
                mechanism_local = findall(small_pivot_mask)
                mechanism_original = perm[mechanism_local]
                mechanism_global = free_dofs[mechanism_original]
                n_mech_trans = count(d -> mod(d-1,6)+1 <= 3, mechanism_global)
                n_mech_rot = n_mechanism - n_mech_trans
                log_msg("[SOLVER] Factorization AUTOSPC: Found $n_mechanism DOFs ($n_mech_trans trans + $n_mech_rot rot), threshold=$pivot_threshold, L_median=$L_median")

                for d in mechanism_global; push!(fixed_dofs, d); end
                free_dofs = sort(collect(setdiff(1:ndof, fixed_dofs)))
                K_ff = K[free_dofs, free_dofs]
                F_ff = F_applied[free_dofs]
                n_free = length(free_dofs)
                log_msg("[SOLVER] Rebuilt system: Fixed=$(length(fixed_dofs)), Free=$n_free")
            end

            local u_result
            try
                F_chol_clean = cholesky(Symmetric(K_ff))
                u_result = F_chol_clean \ F_ff
                log_msg("[SOLVER] Clean Cholesky succeeded after AUTOSPC ($n_mechanism mechanism DOFs)")
            catch e2
                log_msg("[SOLVER] Clean Cholesky failed: $e2. Using LU factorization...")
                F_lu = lu(K_ff)
                u_result = F_lu \ F_ff
            end
            u_result
        end
    else
        log_msg("[SOLVER] Computing Preconditioner (Smoothed Aggregation)...")
        ml = smoothed_aggregation(K_ff)
        P = aspreconditioner(ml)
        log_msg("[SOLVER] Solving Linear System (CG + AMG)...")
        u_ff = try
            cg(K_ff, F_ff; reltol=1e-8, maxiter=5000, Pl=P)
        catch e
            log_msg("[SOLVER] CG Failed. Trying MINRES...")
            minres(K_ff, F_ff; reltol=1e-8, maxiter=5000)
        end
    end

    # Report residual
    r_solve = K_ff * u_ff - F_ff
    log_msg("[SOLVER] Residual: |r|=$(norm(r_solve)), |r|/|F|=$(norm(r_solve)/max(norm(F_ff),1e-30))")

    log_msg("[SOLVER] Post-Processing...")
    u_global = zeros(ndof)
    u_global[free_dofs] = u_ff

    # RBE3 displacement recovery
    n_rbe3_recovered = 0
    for (dep_dof, pairs) in rbe3_map
        u_avg = 0.0
        for (ind_dof, coeff) in pairs
            u_avg += coeff * u_global[ind_dof]
        end
        u_global[dep_dof] = u_avg
        n_rbe3_recovered += 1
    end
    if n_rbe3_recovered > 0
        log_msg("[SOLVER] RBE3: Recovered $n_rbe3_recovered dependent DOFs")
    end

    return u_global, fixed_dofs
end
