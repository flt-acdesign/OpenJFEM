# assembly.jl — Global stiffness matrix assembly

function assemble_stiffness(model)
    log_msg("[SOLVER] Indexing...")
    ids = sort(collect(keys(model["GRIDs"])), by=x->parse(Int,x))
    n_nodes = length(ids)
    id_map = Dict(parse(Int, k)=>i for (i,k) in enumerate(ids))
    ndof = n_nodes * 6

    node_R = Vector{Matrix{Float64}}(undef, n_nodes)
    node_coords = zeros(n_nodes, 3)

    for (sid, g) in model["GRIDs"]
        idx = id_map[g["ID"]]
        node_coords[idx, :] = g["X"]
        cid = g["CD"]
        if cid == 0
            node_R[idx] = Matrix(1.0I, 3, 3)
        elseif haskey(model["CORDs"], string(cid))
            c = model["CORDs"][string(cid)]
            node_R[idx] = hcat(c["U"], c["V"], c["W"])
        else
             node_R[idx] = Matrix(1.0I, 3, 3)
        end
    end

    snorm_normals = compute_snorm_normals(model, id_map, node_coords)

    cshells = model["CSHELLs"]
    cbars = model["CBARs"]
    crods = get(model, "CRODs", Dict())
    rbe2s = get(model, "RBE2s", Dict())
    pshells = model["PSHELLs"]; pbarls = model["PBARLs"]; mats = model["MATs"]
    k6rot = get(model, "PARAM_K6ROT", 100.0)

    nt = Threads.maxthreadid()
    log_msg("[SOLVER] Computing Element Stiffness ($(Threads.nthreads()) threads)...")

    # Convert id_map Dict to dense Vector
    max_nid = maximum(keys(id_map))
    id_vec = zeros(Int, max_nid)
    for (nid, idx) in id_map
        id_vec[nid] = idx
    end

    # Convert snorm_normals Dict to arrays
    snorm_vec = fill(SVector(0.0, 0.0, 0.0), n_nodes)
    snorm_has = falses(n_nodes)
    for (idx, nrm) in snorm_normals
        snorm_vec[idx] = nrm
        snorm_has[idx] = true
    end

    shell_list = collect(values(cshells))
    n_shells = length(shell_list)

    # Pass 1: count QUAD4 and TRIA3 elements
    n_q4 = 0; n_t3 = 0
    for ei in 1:n_shells
        el = shell_list[ei]
        pid = string(el["PID"])
        if !haskey(pshells, pid); continue; end
        prop = pshells[pid]
        mid = string(prop["MID"])
        if !haskey(mats, mid); continue; end
        nids = el["NODES"]; n = length(nids)
        valid = true
        for k in 1:n
            nid = nids[k]
            if nid < 1 || nid > max_nid || id_vec[nid] == 0; valid = false; break; end
        end
        if !valid; continue; end
        if n == 4; n_q4 += 1; elseif n == 3; n_t3 += 1; end
    end

    # Pre-allocate QUAD4 flat arrays
    q4_idx     = Matrix{Int}(undef, n_q4, 4)
    q4_h       = Vector{Float64}(undef, n_q4)
    q4_br      = Vector{Float64}(undef, n_q4)
    q4_Eref    = Vector{Float64}(undef, n_q4)
    q4_Cm_flat = zeros(3, 3, n_q4)
    q4_Cb_flat = zeros(3, 3, n_q4)
    q4_Cs_flat = zeros(2, 2, n_q4)
    q4_Bmb_flat = zeros(3, 3, n_q4)
    q4_has_Bmb = falses(n_q4)

    # Pre-allocate TRIA3 arrays
    t3_idx     = Matrix{Int}(undef, n_t3, 3)
    t3_h       = Vector{Float64}(undef, n_t3)
    t3_br      = Vector{Float64}(undef, n_t3)
    t3_tst     = Vector{Float64}(undef, n_t3)
    t3_Eref    = Vector{Float64}(undef, n_t3)
    t3_Cm      = Vector{Matrix{Float64}}(undef, n_t3)
    t3_Cb      = Vector{Matrix{Float64}}(undef, n_t3)
    t3_Cs      = Vector{Matrix{Float64}}(undef, n_t3)
    t3_Bmb     = Vector{Union{Nothing, Matrix{Float64}}}(undef, n_t3)

    # Pass 2: fill arrays
    iq4 = 0; it3 = 0
    for ei in 1:n_shells
        el = shell_list[ei]
        pid = string(el["PID"])
        if !haskey(pshells, pid); continue; end
        prop = pshells[pid]
        mid = string(prop["MID"])
        if !haskey(mats, mid); continue; end
        mat = mats[mid]

        nids = el["NODES"]
        n = length(nids)
        h = prop["T"]
        br = get(prop, "BEND_RATIO", 1.0)
        tst = get(prop, "TS_T", 5.0/6.0)
        is_pcomp_clt = get(prop, "TYPE", "") == "PCOMP_CLT" && haskey(prop, "Cm")
        is_ortho = !is_pcomp_clt && get(mat, "TYPE", "") == "MAT8" && haskey(mat, "E1") && haskey(mat, "E2")
        is_mat2  = !is_pcomp_clt && !is_ortho && get(mat, "TYPE", "") == "MAT2" && haskey(mat, "G11")

        valid = true
        for k in 1:n
            nid = nids[k]
            if nid < 1 || nid > max_nid || id_vec[nid] == 0; valid = false; break; end
        end
        if !valid; continue; end

        local Cm_e::Matrix{Float64}, Cb_e::Matrix{Float64}, Cs_e::Matrix{Float64}
        local Bmb_e::Union{Nothing, Matrix{Float64}}
        local E_ref_e::Float64

        if is_pcomp_clt
            Cm_e = prop["Cm"]; Cb_e = prop["Cb"]; Cs_e = prop["Cs"]
            Bmb_e = get(prop, "Bmb", nothing)
            E_ref_e = n == 4 ? get(prop, "E_ref", mat["E"]) : mat["G"]
        elseif is_ortho
            E1 = mat["E1"]; E2 = mat["E2"]; nu12 = mat["NU12"]; G12 = mat["G12"]
            nu21 = nu12 * E2 / max(E1, 1e-30)
            denom = 1.0 - nu12 * nu21
            Q11 = E1/denom; Q22 = E2/denom; Q12 = nu12*E2/denom; Q66 = G12
            el_theta = deg2rad(Float64(get(el, "THETA", 0.0)))
            Q16 = 0.0; Q26 = 0.0
            if abs(el_theta) > 1e-10
                ct = cos(el_theta); st = sin(el_theta); c2 = ct^2; s2 = st^2
                Q11r = Q11*c2^2 + 2*(Q12+2*Q66)*c2*s2 + Q22*s2^2
                Q22r = Q11*s2^2 + 2*(Q12+2*Q66)*c2*s2 + Q22*c2^2
                Q12r = (Q11+Q22-4*Q66)*c2*s2 + Q12*(c2^2+s2^2)
                Q16 = (Q11-Q12-2*Q66)*ct*st*c2 + (Q12-Q22+2*Q66)*ct*st*s2
                Q26 = (Q11-Q12-2*Q66)*ct*st*s2 + (Q12-Q22+2*Q66)*ct*st*c2
                Q66r = (Q11+Q22-2*Q12-2*Q66)*c2*s2 + Q66*(c2^2+s2^2)
                Q11 = Q11r; Q22 = Q22r; Q12 = Q12r; Q66 = Q66r
            end
            Cm_e = h .* [Q11 Q12 Q16; Q12 Q22 Q26; Q16 Q26 Q66]
            Cb_e = br * (h^3/12.0) .* [Q11 Q12 Q16; Q12 Q22 Q26; Q16 Q26 Q66]
            G1Z = get(mat, "G1Z", 0.0); G2Z = get(mat, "G2Z", 0.0)
            if G1Z <= 0.0; G1Z = G12; end; if G2Z <= 0.0; G2Z = G12; end
            if abs(el_theta) > 1e-10
                ct = cos(el_theta); st = sin(el_theta)
                Cs_e = tst*h .* [ct^2*G1Z+st^2*G2Z ct*st*(G1Z-G2Z); ct*st*(G1Z-G2Z) st^2*G1Z+ct^2*G2Z]
            else
                Cs_e = tst*h .* [G1Z 0.0; 0.0 G2Z]
            end
            Bmb_e = nothing
            E_ref_e = n == 4 ? max(E1, E2) : G12
        elseif is_mat2
            G11 = mat["G11"]; G12m = mat["G12"]; G13 = mat["G13"]
            G22 = mat["G22"]; G23 = mat["G23"]; G33 = mat["G33"]
            Cm_e = h .* [G11 G12m G13; G12m G22 G23; G13 G23 G33]
            Cb_e = br * (h^3/12.0) .* [G11 G12m G13; G12m G22 G23; G13 G23 G33]
            Cs_e = tst * h .* [G33 0.0; 0.0 G33]
            Bmb_e = nothing
            E_ref_e = n == 4 ? max(G11, G22) : G33
        else
            E_val = mat["E"]; nu_val = mat["NU"]
            const_mem = E_val * h / (1 - nu_val^2)
            Cm_e = const_mem .* [1.0 nu_val 0.0; nu_val 1.0 0.0; 0.0 0.0 (1-nu_val)/2]
            const_bend = br * (E_val * h^3) / (12 * (1 - nu_val^2))
            Cb_e = const_bend .* [1.0 nu_val 0.0; nu_val 1.0 0.0; 0.0 0.0 (1-nu_val)/2]
            G_val = E_val / (2*(1+nu_val))
            k_shear = tst * G_val * h
            Cs_e = k_shear .* [1.0 0.0; 0.0 1.0]
            Bmb_e = nothing
            E_ref_e = n == 4 ? E_val : G_val
        end

        if n == 4
            iq4 += 1
            i1 = id_vec[nids[1]]; i2 = id_vec[nids[2]]; i3 = id_vec[nids[3]]; i4 = id_vec[nids[4]]
            q4_idx[iq4,1] = i1; q4_idx[iq4,2] = i2; q4_idx[iq4,3] = i3; q4_idx[iq4,4] = i4
            q4_h[iq4] = h; q4_br[iq4] = br; q4_Eref[iq4] = E_ref_e
            for j in 1:3, i in 1:3
                q4_Cm_flat[i,j,iq4] = Cm_e[i,j]
                q4_Cb_flat[i,j,iq4] = Cb_e[i,j]
            end
            for j in 1:2, i in 1:2
                q4_Cs_flat[i,j,iq4] = Cs_e[i,j]
            end
            if Bmb_e !== nothing
                q4_has_Bmb[iq4] = true
                for j in 1:3, i in 1:3; q4_Bmb_flat[i,j,iq4] = Bmb_e[i,j]; end
            end
        elseif n == 3
            it3 += 1
            i1 = id_vec[nids[1]]; i2 = id_vec[nids[2]]; i3 = id_vec[nids[3]]
            t3_idx[it3,1] = i1; t3_idx[it3,2] = i2; t3_idx[it3,3] = i3
            t3_h[it3] = h; t3_br[it3] = br; t3_tst[it3] = tst; t3_Eref[it3] = E_ref_e
            t3_Cm[it3] = Cm_e; t3_Cb[it3] = Cb_e; t3_Cs[it3] = Cs_e; t3_Bmb[it3] = Bmb_e
        end
    end

    log_msg("[SOLVER] Pre-extracted $n_q4 QUAD4 + $n_t3 TRIA3 elements (from $n_shells total)")

    # Convert node_R to flat 3D array
    node_R_flat = zeros(3, 3, n_nodes)
    for i in 1:n_nodes
        for r in 1:3, c in 1:3
            node_R_flat[r, c, i] = node_R[i][r, c]
        end
    end

    # --- PARALLEL QUAD4 ASSEMBLY ---
    per_thread_ws = [FEM.create_quad4_workspace() for _ in 1:nt]

    all_I = Vector{Int}(undef, n_q4 * 576)
    all_J = Vector{Int}(undef, n_q4 * 576)
    all_V = Vector{Float64}(undef, n_q4 * 576)

    prev_blas_threads = LinearAlgebra.BLAS.get_num_threads()
    LinearAlgebra.BLAS.set_num_threads(1)

    sep_T       = [zeros(24,24) for _ in 1:nt]
    sep_tmp     = [zeros(24,24) for _ in 1:nt]
    sep_global  = [zeros(24,24) for _ in 1:nt]
    sep_dofs    = [Vector{Int}(undef, 24) for _ in 1:nt]
    sep_lc      = [zeros(4,2) for _ in 1:nt]
    sep_Cm      = [zeros(3,3) for _ in 1:nt]
    sep_Cb      = [zeros(3,3) for _ in 1:nt]
    sep_Cs      = [zeros(2,2) for _ in 1:nt]
    sep_Bmb     = [zeros(3,3) for _ in 1:nt]

    Threads.@threads :static for ei in 1:n_q4
        tid = Threads.threadid()

        i1 = q4_idx[ei,1]; i2 = q4_idx[ei,2]; i3 = q4_idx[ei,3]; i4 = q4_idx[ei,4]

        p1 = SVector{3}(node_coords[i1,1], node_coords[i1,2], node_coords[i1,3])
        p2 = SVector{3}(node_coords[i2,1], node_coords[i2,2], node_coords[i2,3])
        p3 = SVector{3}(node_coords[i3,1], node_coords[i3,2], node_coords[i3,3])
        p4 = SVector{3}(node_coords[i4,1], node_coords[i4,2], node_coords[i4,3])

        v1, v2, v3 = shell_element_frame_fast(p1, p2, p3, p4, 4)

        # SNORM adjustment
        n_avg = SVector(0.0, 0.0, 0.0); nc = 0
        for idx in (i1, i2, i3, i4)
            if snorm_has[idx]; n_avg = n_avg + snorm_vec[idx]; nc += 1; end
        end
        if nc > 0
            n_avg_s = n_avg / nc; len_s = norm(n_avg_s)
            if len_s > 1e-12
                v3n = SVector{3}(n_avg_s / len_s)
                if dot(v3n, v3) < 0.0; v3n = -v3n; end
                v1p = v1 - dot(v1, v3n) * v3n; v1l = norm(v1p)
                if v1l > 1e-12
                    v1n = SVector{3}(v1p / v1l)
                else
                    v2p = v2 - dot(v2, v3n) * v3n; v1n = SVector{3}(normalize(v2p))
                end
                v1, v2, v3 = v1n, SVector{3}(cross(v3n, v1n)), v3n
            end
        end

        lc = sep_lc[tid]
        c = (p1 + p2 + p3 + p4) / 4.0
        lc[1,1] = dot(p1-c, v1); lc[1,2] = dot(p1-c, v2)
        lc[2,1] = dot(p2-c, v1); lc[2,2] = dot(p2-c, v2)
        lc[3,1] = dot(p3-c, v1); lc[3,2] = dot(p3-c, v2)
        lc[4,1] = dot(p4-c, v1); lc[4,2] = dot(p4-c, v2)

        Cm_local = sep_Cm[tid]; Cb_local = sep_Cb[tid]; Cs_local = sep_Cs[tid]
        for j in 1:3, i in 1:3; Cm_local[i,j] = q4_Cm_flat[i,j,ei]; Cb_local[i,j] = q4_Cb_flat[i,j,ei]; end
        for j in 1:2, i in 1:2; Cs_local[i,j] = q4_Cs_flat[i,j,ei]; end
        Bmb_local = nothing
        if q4_has_Bmb[ei]
            Bmb_local = sep_Bmb[tid]
            for j in 1:3, i in 1:3; Bmb_local[i,j] = q4_Bmb_flat[i,j,ei]; end
        end

        ws_stiff = per_thread_ws[tid]
        Ke_t = FEM.stiffness_quad4_matrices(lc, Cm_local, Cb_local, Cs_local,
                    q4_h[ei], q4_Eref[ei]; bend_ratio=q4_br[ei], k6rot=k6rot, Bmb=Bmb_local, ws=ws_stiff)

        Rel_t = @SMatrix [v1[1] v1[2] v1[3]; v2[1] v2[2] v2[3]; v3[1] v3[2] v3[3]]
        T_t = sep_T[tid]; fill!(T_t, 0.0)
        @inbounds @fastmath for k in 1:4
            idx = k == 1 ? i1 : k == 2 ? i2 : k == 3 ? i3 : i4
            base = (k-1)*6
            for rr in 1:3, cc in 1:3
                val = Rel_t[rr,1]*node_R_flat[1,cc,idx] + Rel_t[rr,2]*node_R_flat[2,cc,idx] + Rel_t[rr,3]*node_R_flat[3,cc,idx]
                T_t[base+rr, base+cc] = val
                T_t[base+3+rr, base+3+cc] = val
            end
        end
        tmp_t = sep_tmp[tid]; fill!(tmp_t, 0.0)
        out_t = sep_global[tid]; fill!(out_t, 0.0)
        @inbounds @fastmath for jj in 1:24, ll in 1:24
            val = T_t[ll, jj]
            if val != 0.0
                for ii in 1:24; tmp_t[ii, jj] += Ke_t[ii, ll] * val; end
            end
        end
        @inbounds @fastmath for jj in 1:24, ll in 1:24
            val = tmp_t[ll, jj]
            if val != 0.0
                for ii in 1:24; out_t[ii, jj] += T_t[ll, ii] * val; end
            end
        end

        dofs = sep_dofs[tid]
        for k in 1:4
            idx = k == 1 ? i1 : k == 2 ? i2 : k == 3 ? i3 : i4
            b = (idx-1)*6
            for d in 1:6; dofs[(k-1)*6+d] = b+d; end
        end
        base = (ei-1)*576; cnt = 0
        for cc in 1:24, rr in 1:24
            cnt += 1
            all_I[base+cnt] = dofs[rr]
            all_J[base+cnt] = dofs[cc]
            all_V[base+cnt] = out_t[rr,cc]
        end
    end

    LinearAlgebra.BLAS.set_num_threads(prev_blas_threads)

    I_idx = Vector{Int}(undef, 0); J_idx = Vector{Int}(undef, 0); V_val = Vector{Float64}(undef, 0)
    est_total = length(all_V) + n_t3*324 + length(cbars)*144 + length(crods)*144
    sizehint!(I_idx, est_total); sizehint!(J_idx, est_total); sizehint!(V_val, est_total)
    append!(I_idx, all_I); append!(J_idx, all_J); append!(V_val, all_V)
    all_I = Int[]; all_J = Int[]; all_V = Float64[]

    # --- SEQUENTIAL TRIA3 ASSEMBLY ---
    lc_buf = zeros(3, 2)
    T_buf = zeros(18, 18)
    dofs_t3 = Vector{Int}(undef, 18)
    for ei in 1:n_t3
        i1 = t3_idx[ei,1]; i2 = t3_idx[ei,2]; i3 = t3_idx[ei,3]
        p1 = SVector{3}(node_coords[i1,1], node_coords[i1,2], node_coords[i1,3])
        p2 = SVector{3}(node_coords[i2,1], node_coords[i2,2], node_coords[i2,3])
        p3 = SVector{3}(node_coords[i3,1], node_coords[i3,2], node_coords[i3,3])
        v1, v2, v3 = shell_element_frame_fast(p1, p2, p3, SVector{3}(0.0,0.0,0.0), 3)

        # SNORM adjustment
        n_avg = SVector(0.0, 0.0, 0.0); nc = 0
        for idx in (i1, i2, i3)
            if snorm_has[idx]; n_avg = n_avg + snorm_vec[idx]; nc += 1; end
        end
        if nc > 0
            n_avg_s = n_avg / nc; len_s = norm(n_avg_s)
            if len_s > 1e-12
                v3n = SVector{3}(n_avg_s / len_s)
                if dot(v3n, v3) < 0.0; v3n = -v3n; end
                v1p = v1 - dot(v1, v3n) * v3n; v1l = norm(v1p)
                if v1l > 1e-12
                    v1n = SVector{3}(v1p / v1l)
                else
                    v2p = v2 - dot(v2, v3n) * v3n; v1n = SVector{3}(normalize(v2p))
                end
                v1, v2, v3 = v1n, SVector{3}(cross(v3n, v1n)), v3n
            end
        end

        c = (p1 + p2 + p3) / 3.0
        lc_buf[1,1] = dot(p1-c, v1); lc_buf[1,2] = dot(p1-c, v2)
        lc_buf[2,1] = dot(p2-c, v1); lc_buf[2,2] = dot(p2-c, v2)
        lc_buf[3,1] = dot(p3-c, v1); lc_buf[3,2] = dot(p3-c, v2)
        Ke_loc = FEM.stiffness_tria3_matrices(lc_buf, t3_Cm[ei], t3_Cb[ei], t3_Cs[ei],
                    t3_h[ei], t3_Eref[ei]; bend_ratio=t3_br[ei], ts_t=t3_tst[ei], k6rot=k6rot, Bmb=t3_Bmb[ei])

        Rel_t = @SMatrix [v1[1] v1[2] v1[3]; v2[1] v2[2] v2[3]; v3[1] v3[2] v3[3]]
        fill!(T_buf, 0.0)
        for k in 1:3
            idx = k == 1 ? i1 : k == 2 ? i2 : i3
            TR = Rel_t * node_R[idx]
            base = (k-1)*6
            T_buf[base+1:base+3, base+1:base+3] = TR
            T_buf[base+4:base+6, base+4:base+6] = TR
        end
        Ke = T_buf' * Ke_loc * T_buf

        for k in 1:3
            idx = k == 1 ? i1 : k == 2 ? i2 : i3
            b = (idx-1)*6
            for d in 1:6; dofs_t3[(k-1)*6+d] = b+d; end
        end
        for cc in 1:18, rr in 1:18
            push!(I_idx, dofs_t3[rr]); push!(J_idx, dofs_t3[cc]); push!(V_val, Ke[rr,cc])
        end
    end

    log_msg("[SOLVER] Shell assembly: $n_q4 QUAD4 (parallel) + $n_t3 TRIA3 (sequential), NZ=$(length(I_idx))")
    T_buf = zeros(24, 24)

    # --- CBARS ---
    for (id, bar) in cbars
        pid = string(bar["PID"])
        if !haskey(pbarls, pid); continue; end
        prop = pbarls[pid]
        mid = string(prop["MID"])
        if !haskey(mats, mid); continue; end
        mat = mats[mid]

        if !haskey(id_map, bar["GA"]) || !haskey(id_map, bar["GB"]); continue; end
        i1, i2 = id_map[bar["GA"]], id_map[bar["GB"]]

        p1 = SVector{3}(node_coords[i1,1], node_coords[i1,2], node_coords[i1,3])
        p2 = SVector{3}(node_coords[i2,1], node_coords[i2,2], node_coords[i2,3])

        wa, wb, has_offset, p1_eff, p2_eff = bar_offsets_and_endpoints(bar, p1, p2)
        L = norm(p2_eff - p1_eff)
        if L < 1e-9; continue; end
        vx = normalize(p2_eff - p1_eff)
        v_ref = resolve_bar_vref(bar, p1, id_map, node_coords)

        if norm(v_ref) < 1e-6
             v_ref = SVector(0.0,0.0,1.0)
             if abs(dot(vx, v_ref)) > 0.9; v_ref = SVector(0.0,1.0,0.0); end
        end
        vz = normalize(cross(vx, v_ref))
        vy = cross(vz, vx)
        Rel_t = vcat(vx', vy', vz')

        fill!(view(T_buf, 1:12, 1:12), 0.0)
        TR1 = Rel_t * node_R[i1]
        TR2 = Rel_t * node_R[i2]
        T_buf[1:3, 1:3] = TR1; T_buf[4:6, 4:6] = TR1
        T_buf[7:9, 7:9] = TR2; T_buf[10:12, 10:12] = TR2
        if has_offset
            S_wa = skew3(wa); S_wb = skew3(wb)
            T_buf[1:3, 4:6] = -Rel_t * S_wa * node_R[i1]
            T_buf[7:9, 10:12] = -Rel_t * S_wb * node_R[i2]
        end

        Iy = get(prop, "I2", get(prop, "I", 0.0))
        Iz = get(prop, "I1", get(prop, "I", 0.0))
        if Iy == 0.0; Iy = Iz; end
        if Iz == 0.0; Iz = Iy; end
        Iyz = Float64(get(prop, "I12", 0.0))
        K1 = get(prop, "K1", 0.0); K2 = get(prop, "K2", 0.0)
        As_y = (K1 > 0.0) ? K1 * prop["A"] : Inf
        As_z = (K2 > 0.0) ? K2 * prop["A"] : Inf
        Ke_loc = FEM.stiffness_frame3d(L, prop["A"], Iy, Iz, prop["J"], mat["E"], mat["G"]; As_y=As_y, As_z=As_z, I12=Iyz)
        T_sub = view(T_buf, 1:12, 1:12)
        Ke = T_sub' * Ke_loc * T_sub

        dofs = [(i1-1)*6+k for k in 1:6]
        append!(dofs, [(i2-1)*6+k for k in 1:6])
        for c in 1:12, r in 1:12
            push!(I_idx, dofs[r]); push!(J_idx, dofs[c]); push!(V_val, Ke[r,c])
        end
    end

    # --- CRODS ---
    prods = get(model, "PRODs", Dict())
    for (id, rod) in crods
        pid = string(rod["PID"])
        if !haskey(prods, pid); continue; end
        prop = prods[pid]
        mid = string(prop["MID"])
        if !haskey(mats, mid); continue; end
        mat = mats[mid]

        if !haskey(id_map, rod["GA"]) || !haskey(id_map, rod["GB"]); continue; end
        i1, i2 = id_map[rod["GA"]], id_map[rod["GB"]]

        p1 = SVector{3}(node_coords[i1,1], node_coords[i1,2], node_coords[i1,3])
        p2 = SVector{3}(node_coords[i2,1], node_coords[i2,2], node_coords[i2,3])
        L = norm(p2-p1)
        if L < 1e-9; continue; end
        vx = normalize(p2-p1)

        ref = abs(vx[3]) < 0.9 ? SVector(0.0,0.0,1.0) : SVector(0.0,1.0,0.0)
        vz = normalize(cross(vx, ref))
        vy = cross(vz, vx)
        Rel_t = vcat(vx', vy', vz')

        Ke_loc = zeros(12, 12)
        EA_L = mat["E"] * prop["A"] / L
        Ke_loc[1,1] = EA_L; Ke_loc[1,7] = -EA_L; Ke_loc[7,1] = -EA_L; Ke_loc[7,7] = EA_L
        GJ_L = mat["G"] * prop["J"] / L
        Ke_loc[4,4] = GJ_L; Ke_loc[4,10] = -GJ_L; Ke_loc[10,4] = -GJ_L; Ke_loc[10,10] = GJ_L

        fill!(view(T_buf, 1:12, 1:12), 0.0)
        TR1 = Rel_t * node_R[i1]; TR2 = Rel_t * node_R[i2]
        T_buf[1:3, 1:3] = TR1; T_buf[4:6, 4:6] = TR1
        T_buf[7:9, 7:9] = TR2; T_buf[10:12, 10:12] = TR2
        T_sub = view(T_buf, 1:12, 1:12)
        Ke = T_sub' * Ke_loc * T_sub

        dofs = [(i1-1)*6+k for k in 1:6]
        append!(dofs, [(i2-1)*6+k for k in 1:6])
        for c in 1:12, r in 1:12
            push!(I_idx, dofs[r]); push!(J_idx, dofs[c]); push!(V_val, Ke[r,c])
        end
    end

    # --- CONROD ---
    conrods = get(model, "CONRODs", Dict())
    for (id, rod) in conrods
        mid = string(rod["MID"])
        if !haskey(mats, mid); continue; end
        mat = mats[mid]
        if !haskey(id_map, rod["GA"]) || !haskey(id_map, rod["GB"]); continue; end
        i1, i2 = id_map[rod["GA"]], id_map[rod["GB"]]
        p1 = SVector{3}(node_coords[i1,1], node_coords[i1,2], node_coords[i1,3])
        p2 = SVector{3}(node_coords[i2,1], node_coords[i2,2], node_coords[i2,3])
        L = norm(p2-p1)
        if L < 1e-9; continue; end
        vx = normalize(p2-p1)
        ref = abs(vx[3]) < 0.9 ? SVector(0.0,0.0,1.0) : SVector(0.0,1.0,0.0)
        vz = normalize(cross(vx, ref))
        vy = cross(vz, vx)
        Rel_t = vcat(vx', vy', vz')
        Ke_loc = zeros(12, 12)
        EA_L = mat["E"] * rod["A"] / L
        Ke_loc[1,1] = EA_L; Ke_loc[1,7] = -EA_L; Ke_loc[7,1] = -EA_L; Ke_loc[7,7] = EA_L
        GJ_L = mat["G"] * rod["J"] / L
        Ke_loc[4,4] = GJ_L; Ke_loc[4,10] = -GJ_L; Ke_loc[10,4] = -GJ_L; Ke_loc[10,10] = GJ_L
        fill!(view(T_buf, 1:12, 1:12), 0.0)
        TR1 = Rel_t * node_R[i1]; TR2 = Rel_t * node_R[i2]
        T_buf[1:3, 1:3] = TR1; T_buf[4:6, 4:6] = TR1
        T_buf[7:9, 7:9] = TR2; T_buf[10:12, 10:12] = TR2
        T_sub = view(T_buf, 1:12, 1:12)
        Ke = T_sub' * Ke_loc * T_sub
        dofs = [(i1-1)*6+k for k in 1:6]
        append!(dofs, [(i2-1)*6+k for k in 1:6])
        for c in 1:12, r in 1:12
            push!(I_idx, dofs[r]); push!(J_idx, dofs[c]); push!(V_val, Ke[r,c])
        end
    end

    # --- CELAS1 ---
    celases = get(model, "CELASs", Dict())
    pelases = get(model, "PELASs", Dict())
    n_springs = 0
    for (id, spring) in celases
        pid = string(spring["PID"])
        if !haskey(pelases, pid); continue; end
        K_spring = pelases[pid]["K"]
        g1 = spring["G1"]; c1 = spring["C1"]
        g2 = spring["G2"]; c2 = spring["C2"]
        if g1 > 0 && haskey(id_map, g1) && c1 > 0
            i1 = id_map[g1]
            dof1 = (i1-1)*6 + c1
            push!(I_idx, dof1); push!(J_idx, dof1); push!(V_val, K_spring)
            if g2 > 0 && haskey(id_map, g2) && c2 > 0
                i2 = id_map[g2]
                dof2 = (i2-1)*6 + c2
                push!(I_idx, dof2); push!(J_idx, dof2); push!(V_val, K_spring)
                push!(I_idx, dof1); push!(J_idx, dof2); push!(V_val, -K_spring)
                push!(I_idx, dof2); push!(J_idx, dof1); push!(V_val, -K_spring)
            end
            n_springs += 1
        end
    end
    if n_springs > 0
        log_msg("[SOLVER] CELAS1: $n_springs spring elements assembled")
    end

    # Save max element stiffness BEFORE constraint processing
    max_elem_stiff = 0.0
    for i in 1:length(V_val)
        if abs(V_val[i]) > max_elem_stiff; max_elem_stiff = abs(V_val[i]); end
    end

    # --- Constraint assembly (RBE2, RBE3, MPC) ---
    rbe3_map, I_idx, J_idx, V_val = assemble_constraints(model, id_map, node_coords, I_idx, J_idx, V_val)

    actual_nz = length(I_idx)
    log_msg("[SOLVER] Creating Sparse Matrix (NZ: $actual_nz)...")
    K = sparse(I_idx, J_idx, V_val, ndof, ndof)
    I_idx = nothing; J_idx = nothing; V_val = nothing; GC.gc()

    return K, id_map, node_coords, ndof, node_R, max_elem_stiff, rbe3_map, snorm_normals
end
