
module FEM

using LinearAlgebra
using Statistics
using StaticArrays

# Pre-allocated workspace for QUAD4 element stiffness computation.
# Allocate once per thread, reuse across all elements to eliminate ~5M heap allocations.
struct Quad4Workspace
    # Accumulated per-element (cleared at start of each element)
    Ke::Matrix{Float64}          # 24×24 element stiffness
    K_ab::Matrix{Float64}        # 24×4  membrane incompatible mode coupling
    K_bb::Matrix{Float64}        # 4×4   membrane incompatible self-coupling
    K_ab_bend::Matrix{Float64}   # 24×4  bending incompatible mode coupling
    K_bb_bend::Matrix{Float64}   # 4×4   bending incompatible self-coupling
    # MITC4 tying point B-matrices (filled per element)
    Bs_tp::Matrix{Float64}       # 4×24  rows: [Bs_xi_A; Bs_xi_C; Bs_eta_B; Bs_eta_D]
    Bs_row::Vector{Float64}      # 24    temporary for tying point computation
    # Gauss-point matrices (reused each GP, cleared with fill!)
    Bm::Matrix{Float64}          # 3×24  membrane strain-displacement
    Bb::Matrix{Float64}          # 3×24  bending strain-displacement
    Bd::Matrix{Float64}          # 1×24  drilling B-matrix
    Bi::Matrix{Float64}          # 3×4   membrane incompatible mode B-matrix
    Bi_bend::Matrix{Float64}     # 3×4   bending incompatible mode B-matrix
    Bs_cov::Matrix{Float64}      # 2×24  covariant shear B-matrix
    # Temporaries for in-place mul!
    tmp3x24::Matrix{Float64}     # for Cm*Bm, Cb*Bb products
    tmp3x4::Matrix{Float64}      # for Cm*Bi, Cb*Bi_bend products
    tmp2x24::Matrix{Float64}     # for Cs_cov*Bs_cov products
    tmp2x2::Matrix{Float64}      # for Cs_cov = invJ'*Cs*invJ
    # B matrix coupling workspace (reused when Bmb != nothing)
    K_ab_cross::Matrix{Float64}      # 24×4
    K_ab_bend_cross::Matrix{Float64} # 24×4
    K_mb_incomp::Matrix{Float64}     # 4×4
    # Coordinate transform workspace (used in Solver.jl assembly)
    Ke_global::Matrix{Float64}   # 24×24  T'*Ke*T result
    tmp24x24::Matrix{Float64}    # 24×24  temporary for transform
    Rel_t::Matrix{Float64}       # 3×3   element-to-global rotation
    # Thread-local constitutive matrix buffers (copied from flat arrays, avoids Vector{Matrix} reads)
    Cm_buf::Matrix{Float64}      # 3×3   membrane constitutive
    Cb_buf::Matrix{Float64}      # 3×3   bending constitutive
    Cs_buf::Matrix{Float64}      # 2×2   shear constitutive
    Bmb_buf::Matrix{Float64}     # 3×3   membrane-bending coupling
    # Assembly buffers (used in Solver.jl parallel loop)
    T_buf::Matrix{Float64}       # 24×24 transformation matrix
    lc::Matrix{Float64}          # 4×2   local coordinates
    dofs::Vector{Int}            # 24    DOF indices
end

function create_quad4_workspace()
    Quad4Workspace(
        zeros(24,24), zeros(24,4), zeros(4,4), zeros(24,4), zeros(4,4),   # Ke, K_ab, K_bb, K_ab_bend, K_bb_bend
        zeros(4,24), zeros(24),                                             # Bs_tp, Bs_row
        zeros(3,24), zeros(3,24), zeros(1,24), zeros(3,4), zeros(3,4), zeros(2,24), # Bm..Bs_cov
        zeros(3,24), zeros(3,4), zeros(2,24), zeros(2,2),                  # tmp buffers
        zeros(24,4), zeros(24,4), zeros(4,4),                              # B coupling
        zeros(24,24), zeros(24,24), zeros(3,3),                            # transform
        zeros(3,3), zeros(3,3), zeros(2,2), zeros(3,3),                    # constitutive buffers
        zeros(24,24), zeros(4,2), Vector{Int}(undef, 24)                   # assembly buffers
    )
end

# Thread-safe matrix multiplication replacing BLAS mul! (which is NOT re-entrant on Windows).
# C += alpha * A * B   (A is m×k, B is k×n, C is m×n)
@inline function ts_mul_add!(C, A, B, alpha)
    m, k = size(A)
    _, n = size(B)
    @inbounds @fastmath for j in 1:n
        for l in 1:k
            val = alpha * B[l,j]
            for i in 1:m
                C[i,j] += A[i,l] * val
            end
        end
    end
end

# C += alpha * A' * B   (A is k×m, A' is m×k, B is k×n, C is m×n)
@inline function ts_mul_At_add!(C, A, B, alpha)
    k, m = size(A)
    _, n = size(B)
    @inbounds @fastmath for j in 1:n
        for l in 1:k
            val = alpha * B[l,j]
            for i in 1:m
                C[i,j] += A[l,i] * val
            end
        end
    end
end

# C = A * B   (overwrite, no accumulate)
@inline function ts_mul!(C, A, B)
    m, k = size(A)
    _, n = size(B)
    fill!(C, 0.0)
    @inbounds @fastmath for j in 1:n
        for l in 1:k
            val = B[l,j]
            for i in 1:m
                C[i,j] += A[i,l] * val
            end
        end
    end
end

function stiffness_frame3d(L, A, Iy, Iz, J, E, G; As_y=Inf, As_z=Inf, I12=0.0)
    k = zeros(12, 12)
    if L < 1e-9; return k; end

    X = E * A / L
    k[1,1] = X;  k[1,7] = -X; k[7,1] = -X; k[7,7] = X

    T = G * J / L
    k[4,4] = T;  k[4,10] = -T; k[10,4] = -T; k[10,10] = T

    # Timoshenko shear parameters (Φ=0 reduces to Euler-Bernoulli)
    # NASTRAN: K1=shear area factor for plane 1 (y-dir), K2=for plane 2 (z-dir)
    # As_y = K1*A (shear area in y), As_z = K2*A (shear area in z)
    # xz-plane bending: deflection in z → shear in z → uses As_z (K2*A)
    # xy-plane bending: deflection in y → shear in y → uses As_y (K1*A)
    Phi_y = isinf(As_z) ? 0.0 : 12*E*Iy/(G*As_z*L^2)
    Phi_z = isinf(As_y) ? 0.0 : 12*E*Iz/(G*As_y*L^2)

    # Bending in xz-plane (uses Iy, shear via As_z=K2*A)
    a_y = 12*E*Iy / (L^3*(1+Phi_y))
    b_y = 6*E*Iy / (L^2*(1+Phi_y))
    c_y = (4+Phi_y)*E*Iy / (L*(1+Phi_y))
    d_y = (2-Phi_y)*E*Iy / (L*(1+Phi_y))
    k[3,3] = a_y;  k[3,9] = -a_y; k[9,3] = -a_y; k[9,9] = a_y
    k[3,5] = -b_y; k[3,11] = -b_y; k[5,3] = -b_y; k[11,3] = -b_y
    k[9,5] = b_y;  k[9,11] = b_y;  k[5,9] = b_y;  k[11,9] = b_y
    k[5,5] = c_y;  k[5,11] = d_y;  k[11,5] = d_y;  k[11,11] = c_y

    # Bending in xy-plane (uses Iz, shear via As_y=K1*A)
    a_z = 12*E*Iz / (L^3*(1+Phi_z))
    b_z = 6*E*Iz / (L^2*(1+Phi_z))
    c_z = (4+Phi_z)*E*Iz / (L*(1+Phi_z))
    d_z = (2-Phi_z)*E*Iz / (L*(1+Phi_z))
    k[2,2] = a_z;  k[2,8] = -a_z; k[8,2] = -a_z; k[8,8] = a_z
    k[2,6] = b_z;  k[2,12] = b_z;  k[6,2] = b_z;  k[12,2] = b_z
    k[8,6] = -b_z; k[8,12] = -b_z; k[6,8] = -b_z; k[12,8] = -b_z
    k[6,6] = c_z;  k[6,12] = d_z;  k[12,6] = d_z;  k[12,12] = c_z

    # Cross-coupling from I12 (product of inertia)
    # Couples xy-plane bending {v: 2,8; θz: 6,12} with xz-plane bending {w: 3,9; θy: 5,11}
    if abs(I12) > 0.0
        a_yz = 12*E*I12 / L^3
        b_yz = 6*E*I12 / L^2
        c_yz = 4*E*I12 / L
        d_yz = 2*E*I12 / L
        # v-w coupling: DOFs (2,3), (2,9), (8,3), (8,9)
        k[2,3] += a_yz;  k[3,2] += a_yz
        k[2,9] += -a_yz; k[9,2] += -a_yz
        k[8,3] += -a_yz; k[3,8] += -a_yz
        k[8,9] += a_yz;  k[9,8] += a_yz
        # v-θy coupling: DOFs (2,5), (2,11), (8,5), (8,11)
        k[2,5] += -b_yz; k[5,2] += -b_yz
        k[2,11] += -b_yz; k[11,2] += -b_yz
        k[8,5] += b_yz;  k[5,8] += b_yz
        k[8,11] += b_yz; k[11,8] += b_yz
        # w-θz coupling: DOFs (3,6), (3,12), (9,6), (9,12)
        k[3,6] += b_yz;  k[6,3] += b_yz
        k[3,12] += b_yz; k[12,3] += b_yz
        k[9,6] += -b_yz; k[6,9] += -b_yz
        k[9,12] += -b_yz; k[12,9] += -b_yz
        # θy-θz coupling: DOFs (5,6), (5,12), (11,6), (11,12)
        k[5,6] += -c_yz; k[6,5] += -c_yz
        k[5,12] += -d_yz; k[12,5] += -d_yz
        k[11,6] += -d_yz; k[6,11] += -d_yz
        k[11,12] += -c_yz; k[12,11] += -c_yz
    end

    return k
end

function forces_frame3d(u_elem, L, A, Iy, Iz, J, E, G; As_y=Inf, As_z=Inf, I12=0.0)
    k = stiffness_frame3d(L, A, Iy, Iz, J, E, G; As_y=As_y, As_z=As_z, I12=I12)
    f_local = k * u_elem
    
    return Dict(
        "axial"     => f_local[7],
        "shear_1"   => -f_local[2],
        "shear_2"   => -f_local[3],
        "torque"    => -f_local[4],
        "moment_a1" => -f_local[6],
        "moment_a2" => f_local[5],
        "moment_b1" => f_local[12],
        "moment_b2" => -f_local[11]
    )
end

function stress_frame3d(u_elem, L, A, E)
    return E * (u_elem[7] - u_elem[1]) / L
end



@inline function shape_derivs_quad(xi, eta)
    dN_dxi  = SVector{4}(-0.25*(1-eta), 0.25*(1-eta), 0.25*(1+eta), -0.25*(1+eta))
    dN_deta = SVector{4}(-0.25*(1-xi), -0.25*(1+xi), 0.25*(1+xi), 0.25*(1-xi))
    return dN_dxi, dN_deta
end

function stiffness_quad4(coords, E, nu, h; bend_ratio=1.0, ts_t=5.0/6.0, k6rot=100.0, ws::Union{Nothing,Quad4Workspace}=nothing)
    const_mem = E * h / (1 - nu^2)
    Cm = const_mem .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]

    const_bend = bend_ratio * (E * h^3) / (12 * (1 - nu^2))
    Cb = const_bend .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]

    G = E / (2*(1+nu))
    k_shear = ts_t * G * h
    Cs = k_shear .* [1 0; 0 1]

    return stiffness_quad4_matrices(coords, Cm, Cb, Cs, h, E; bend_ratio=bend_ratio, k6rot=k6rot, ws=ws)
end

# Performance-optimized QUAD4 stiffness with optional pre-allocated workspace.
# When ws is provided, eliminates ALL heap allocations in the hot loop (~5M saved across model).
function stiffness_quad4_matrices(coords, Cm, Cb, Cs, h, E_ref; bend_ratio=1.0, k6rot=100.0, Bmb=nothing, ws::Union{Nothing,Quad4Workspace}=nothing)
    if ws === nothing; ws = create_quad4_workspace(); end

    # Clear accumulated matrices
    fill!(ws.Ke, 0.0)
    fill!(ws.K_ab, 0.0); fill!(ws.K_bb, 0.0)
    fill!(ws.K_ab_bend, 0.0); fill!(ws.K_bb_bend, 0.0)

    # B coupling accumulators (cleared even if Bmb is nothing — branch-free)
    fill!(ws.K_ab_cross, 0.0); fill!(ws.K_ab_bend_cross, 0.0); fill!(ws.K_mb_incomp, 0.0)

    # --- MITC4 transverse shear (Bathe-Dvorkin) tying points ---
    fill!(ws.Bs_tp, 0.0)
    tying_pts = (SVector(0.0, -1.0), SVector(0.0, 1.0), SVector(-1.0, 0.0), SVector(1.0, 0.0))
    for tp_idx in 1:4
        xi_tp, eta_tp = tying_pts[tp_idx][1], tying_pts[tp_idx][2]
        dNr, dNs = shape_derivs_quad(xi_tp, eta_tp)
        # Inline Jacobian 2×2
        J11 = dNr[1]*coords[1,1] + dNr[2]*coords[2,1] + dNr[3]*coords[3,1] + dNr[4]*coords[4,1]
        J12 = dNr[1]*coords[1,2] + dNr[2]*coords[2,2] + dNr[3]*coords[3,2] + dNr[4]*coords[4,2]
        J21 = dNs[1]*coords[1,1] + dNs[2]*coords[2,1] + dNs[3]*coords[3,1] + dNs[4]*coords[4,1]
        J22 = dNs[1]*coords[1,2] + dNs[2]*coords[2,2] + dNs[3]*coords[3,2] + dNs[4]*coords[4,2]

        N1 = 0.25*(1-xi_tp)*(1-eta_tp); N2 = 0.25*(1+xi_tp)*(1-eta_tp)
        N3 = 0.25*(1+xi_tp)*(1+eta_tp); N4 = 0.25*(1-xi_tp)*(1+eta_tp)
        N_tp = SVector(N1, N2, N3, N4)

        fill!(ws.Bs_row, 0.0)
        if tp_idx <= 2  # A,C: compute e_ξz (uses dNr for ∂/∂ξ)
            for k in 1:4
                idx = (k-1)*6
                ws.Bs_row[idx+3] = dNr[k]
                ws.Bs_row[idx+4] = -J12*N_tp[k]
                ws.Bs_row[idx+5] = J11*N_tp[k]
            end
        else  # B,D: compute e_ηz (uses dNs for ∂/∂η)
            for k in 1:4
                idx = (k-1)*6
                ws.Bs_row[idx+3] = dNs[k]
                ws.Bs_row[idx+4] = -J22*N_tp[k]
                ws.Bs_row[idx+5] = J21*N_tp[k]
            end
        end
        # Copy to tying point row: rows 1=A, 2=C, 3=B, 4=D
        @views copyto!(ws.Bs_tp[tp_idx, :], ws.Bs_row)
    end

    pt = 1.0/sqrt(3.0)
    gauss_pts = (SVector(-pt,-pt), SVector(pt,-pt), SVector(pt,pt), SVector(-pt,pt))

    # Hughes-Brezzi drilling
    G_drill = Cm[3,3] / h
    if G_drill < 1e-6; G_drill = E_ref / (2*3.0); end
    alpha_drill = (k6rot / 1e5) * G_drill * h

    @inbounds @fastmath for gp in 1:4
        r, s = gauss_pts[gp][1], gauss_pts[gp][2]
        dNr, dNs = shape_derivs_quad(r, s)

        # Inline Jacobian computation (eliminates [dNr';dNs']*coords allocation)
        J11 = dNr[1]*coords[1,1] + dNr[2]*coords[2,1] + dNr[3]*coords[3,1] + dNr[4]*coords[4,1]
        J12 = dNr[1]*coords[1,2] + dNr[2]*coords[2,2] + dNr[3]*coords[3,2] + dNr[4]*coords[4,2]
        J21 = dNs[1]*coords[1,1] + dNs[2]*coords[2,1] + dNs[3]*coords[3,1] + dNs[4]*coords[4,1]
        J22 = dNs[1]*coords[1,2] + dNs[2]*coords[2,2] + dNs[3]*coords[3,2] + dNs[4]*coords[4,2]
        detJ = J11*J22 - J12*J21
        abs_detJ = abs(detJ)
        if abs_detJ < 1e-12; abs_detJ = 1e-12; end
        inv_det = 1.0 / detJ
        iJ11 = J22*inv_det; iJ12 = -J12*inv_det
        iJ21 = -J21*inv_det; iJ22 = J11*inv_det

        # Fill Bm, Bb, Bd directly from inline dN/dx, dN/dy
        fill!(ws.Bm, 0.0); fill!(ws.Bb, 0.0); fill!(ws.Bd, 0.0)
        for k in 1:4
            dN_dx = iJ11*dNr[k] + iJ12*dNs[k]
            dN_dy = iJ21*dNr[k] + iJ22*dNs[k]
            idx = (k-1)*6
            ws.Bm[1, idx+1] = dN_dx;  ws.Bm[2, idx+2] = dN_dy
            ws.Bm[3, idx+1] = dN_dy;  ws.Bm[3, idx+2] = dN_dx
            ws.Bb[1, idx+5] = dN_dx;  ws.Bb[2, idx+4] = -dN_dy
            ws.Bb[3, idx+5] = dN_dy;  ws.Bb[3, idx+4] = -dN_dx
            # Drilling: Bd = [0.5*dN/dy, -0.5*dN/dx, 0, 0, 0, N_k] per node
            N_k = 0.25*(1 + (k==2||k==3 ? r : -r))*(1 + (k>=3 ? s : -s))
            ws.Bd[1, idx+1] = 0.5*dN_dy
            ws.Bd[1, idx+2] = -0.5*dN_dx
            ws.Bd[1, idx+6] = N_k
        end

        # Incompatible mode bubble function derivatives
        dphi1_dx = iJ11*(-2.0*r);  dphi1_dy = iJ21*(-2.0*r)
        dphi2_dx = iJ12*(-2.0*s);  dphi2_dy = iJ22*(-2.0*s)

        # Fill Bi (membrane incompatible)
        fill!(ws.Bi, 0.0)
        ws.Bi[1,1]=dphi1_dx; ws.Bi[3,1]=dphi1_dy; ws.Bi[2,2]=dphi1_dy; ws.Bi[3,2]=dphi1_dx
        ws.Bi[1,3]=dphi2_dx; ws.Bi[3,3]=dphi2_dy; ws.Bi[2,4]=dphi2_dy; ws.Bi[3,4]=dphi2_dx

        # Fill Bi_bend (bending incompatible)
        fill!(ws.Bi_bend, 0.0)
        ws.Bi_bend[2,1]=-dphi1_dy; ws.Bi_bend[2,2]=-dphi2_dy
        ws.Bi_bend[1,3]=dphi1_dx;  ws.Bi_bend[1,4]=dphi2_dx
        ws.Bi_bend[3,1]=-dphi1_dx; ws.Bi_bend[3,2]=-dphi2_dx
        ws.Bi_bend[3,3]=dphi1_dy;  ws.Bi_bend[3,4]=dphi2_dy

        # MITC4 interpolated covariant shear at this Gauss point
        w_A = 0.5*(1-s); w_C = 0.5*(1+s); w_B = 0.5*(1-r); w_D = 0.5*(1+r)
        for j in 1:24
            ws.Bs_cov[1,j] = w_A*ws.Bs_tp[1,j] + w_C*ws.Bs_tp[2,j]
            ws.Bs_cov[2,j] = w_B*ws.Bs_tp[3,j] + w_D*ws.Bs_tp[4,j]
        end

        # Cs_cov = invJ' * Cs * invJ (inline 2×2)
        t11 = Cs[1,1]*iJ11 + Cs[1,2]*iJ21; t12 = Cs[1,1]*iJ12 + Cs[1,2]*iJ22
        t21 = Cs[2,1]*iJ11 + Cs[2,2]*iJ21; t22 = Cs[2,1]*iJ12 + Cs[2,2]*iJ22
        ws.tmp2x2[1,1] = iJ11*t11 + iJ21*t21; ws.tmp2x2[1,2] = iJ11*t12 + iJ21*t22
        ws.tmp2x2[2,1] = iJ12*t11 + iJ22*t21; ws.tmp2x2[2,2] = iJ12*t12 + iJ22*t22

        # === In-place stiffness accumulation (thread-safe, BLAS-free) ===
        # Ke += abs_detJ * Bm' * Cm * Bm
        ts_mul!(ws.tmp3x24, Cm, ws.Bm)
        ts_mul_At_add!(ws.Ke, ws.Bm, ws.tmp3x24, abs_detJ)
        # Ke += abs_detJ * Bb' * Cb * Bb
        ts_mul!(ws.tmp3x24, Cb, ws.Bb)
        ts_mul_At_add!(ws.Ke, ws.Bb, ws.tmp3x24, abs_detJ)
        # Ke += abs_detJ * Bs_cov' * Cs_cov * Bs_cov
        ts_mul!(ws.tmp2x24, ws.tmp2x2, ws.Bs_cov)
        ts_mul_At_add!(ws.Ke, ws.Bs_cov, ws.tmp2x24, abs_detJ)
        # Ke += abs_detJ * alpha_drill * Bd' * Bd
        ts_mul_At_add!(ws.Ke, ws.Bd, ws.Bd, abs_detJ * alpha_drill)

        # B matrix coupling: Ke += abs_detJ * (Bm'*B*Bb + Bb'*B*Bm)
        if Bmb !== nothing
            ts_mul!(ws.tmp3x24, Bmb, ws.Bb)
            ts_mul_At_add!(ws.Ke, ws.Bm, ws.tmp3x24, abs_detJ)
            ts_mul!(ws.tmp3x24, Bmb, ws.Bm)
            ts_mul_At_add!(ws.Ke, ws.Bb, ws.tmp3x24, abs_detJ)
        end

        # Incompatible mode coupling: K_ab += abs_detJ * Bm' * Cm * Bi
        ts_mul!(ws.tmp3x4, Cm, ws.Bi)
        ts_mul_At_add!(ws.K_ab, ws.Bm, ws.tmp3x4, abs_detJ)
        # K_bb += abs_detJ * Bi' * Cm * Bi (reuse tmp3x4 = Cm*Bi)
        ts_mul_At_add!(ws.K_bb, ws.Bi, ws.tmp3x4, abs_detJ)
        # K_ab_bend += abs_detJ * Bb' * Cb * Bi_bend
        ts_mul!(ws.tmp3x4, Cb, ws.Bi_bend)
        ts_mul_At_add!(ws.K_ab_bend, ws.Bb, ws.tmp3x4, abs_detJ)
        # K_bb_bend += abs_detJ * Bi_bend' * Cb * Bi_bend (reuse tmp3x4)
        ts_mul_At_add!(ws.K_bb_bend, ws.Bi_bend, ws.tmp3x4, abs_detJ)

        # B coupling cross-terms for incompatible modes (accumulated during main loop)
        if Bmb !== nothing
            ts_mul!(ws.tmp3x4, Bmb, ws.Bi)
            ts_mul_At_add!(ws.K_ab_cross, ws.Bb, ws.tmp3x4, abs_detJ)
            ts_mul!(ws.tmp3x4, Bmb, ws.Bi_bend)
            ts_mul_At_add!(ws.K_ab_bend_cross, ws.Bm, ws.tmp3x4, abs_detJ)
            ts_mul!(ws.tmp3x4, Bmb, ws.Bi_bend)
            ts_mul_At_add!(ws.K_mb_incomp, ws.Bi, ws.tmp3x4, abs_detJ)
        end
    end

    # Static condensation (BLAS-free for thread safety)
    if Bmb !== nothing
        # Combined 8×8 condensation (B coupling creates cross-coupling between membrane/bending incomp modes)
        K_ab_full = hcat(ws.K_ab .+ ws.K_ab_cross, ws.K_ab_bend .+ ws.K_ab_bend_cross)
        K_bb_full = [ws.K_bb ws.K_mb_incomp; ws.K_mb_incomp' ws.K_bb_bend]
        inv_Kbb = Matrix(inv(SMatrix{8,8}(K_bb_full)))
        tmp8x24 = zeros(8, 24)
        ts_mul!(tmp8x24, inv_Kbb, K_ab_full')
        @inbounds @fastmath for j in 1:24, i in 1:24
            s = 0.0
            for l in 1:8; s += K_ab_full[i,l] * tmp8x24[l,j]; end
            ws.Ke[i,j] -= s
        end
    else
        # Standard separate 4×4 condensation (no B coupling)
        inv_Kbb_m = Matrix(inv(SMatrix{4,4}(ws.K_bb)))
        inv_Kbb_b = Matrix(inv(SMatrix{4,4}(ws.K_bb_bend)))
        # Ke -= K_ab * inv(K_bb) * K_ab'  (24×4 * 4×4 * 4×24 = 24×24)
        @inbounds @fastmath for j in 1:24, i in 1:24
            sm = 0.0; sb = 0.0
            for l in 1:4
                tmp_m = 0.0; tmp_b = 0.0
                for q in 1:4
                    tmp_m += inv_Kbb_m[l,q] * ws.K_ab[j,q]
                    tmp_b += inv_Kbb_b[l,q] * ws.K_ab_bend[j,q]
                end
                sm += ws.K_ab[i,l] * tmp_m
                sb += ws.K_ab_bend[i,l] * tmp_b
            end
            ws.Ke[i,j] -= sm + sb
        end
    end

    return ws.Ke
end

function compute_principal_2d(s11, s22, s12)
    s_avg = (s11 + s22) / 2.0
    radius = sqrt(((s11 - s22) / 2.0)^2 + s12^2)
    return s_avg + radius, s_avg - radius
end

function stress_strain_quad4(coords, u_elem, E, nu, h, t_shell; bend_ratio=1.0, Cm_override=nothing)
    const_mem = E / (1 - nu^2)
    D_mem = const_mem .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    # For PCOMP elements, use CLT Cm for incompatible mode condensation
    # (must match the Cm used in stiffness assembly for consistent strain recovery)
    Cm = isnothing(Cm_override) ? D_mem * h : Cm_override

    dNr, dNs = shape_derivs_quad(0.0, 0.0)
    J = [dNr'; dNs'] * coords
    invJ = inv(J); dN_dxy = invJ * [dNr'; dNs']

    Bm = zeros(3, 24); Bb = zeros(3, 24)

    for k in 1:4
        idx = (k-1)*6
        Bm[1, idx+1]=dN_dxy[1,k]; Bm[2, idx+2]=dN_dxy[2,k]
        Bm[3, idx+1]=dN_dxy[2,k]; Bm[3, idx+2]=dN_dxy[1,k]
        Bb[1, idx+5] = dN_dxy[1,k];
        Bb[2, idx+4] = -dN_dxy[2,k];
        Bb[3, idx+5] = dN_dxy[2,k];
        Bb[3, idx+4] = -dN_dxy[1,k];
    end

    # Recover incompatible mode amplitudes via static condensation
    # α = -K_bb^{-1} * K_ba * u  (K_ba = K_ab')
    # Recompute K_ab and K_bb (membrane incompatible coupling)
    K_ab_sr = zeros(24, 4); K_bb_sr = zeros(4, 4)
    pt = 1.0/sqrt(3.0)
    gauss_pts = [-pt -pt; pt -pt; pt pt; -pt pt]
    for i in 1:4
        r, s = gauss_pts[i,1], gauss_pts[i,2]
        dNr_g, dNs_g = shape_derivs_quad(r, s)
        J_g = [dNr_g'; dNs_g'] * coords
        detJ_g = abs(det(J_g))
        if detJ_g < 1e-12; detJ_g = 1e-12; end
        iJ = inv(J_g)
        dN_dxy_g = iJ * [dNr_g'; dNs_g']

        Bm_g = zeros(3, 24)
        for k in 1:4
            idx = (k-1)*6
            Bm_g[1, idx+1] = dN_dxy_g[1,k]; Bm_g[2, idx+2] = dN_dxy_g[2,k]
            Bm_g[3, idx+1] = dN_dxy_g[2,k]; Bm_g[3, idx+2] = dN_dxy_g[1,k]
        end

        dphi1_dx = iJ[1,1]*(-2.0*r); dphi1_dy = iJ[2,1]*(-2.0*r)
        dphi2_dx = iJ[1,2]*(-2.0*s); dphi2_dy = iJ[2,2]*(-2.0*s)

        Bi = zeros(3, 4)
        Bi[1, 1] = dphi1_dx;  Bi[3, 1] = dphi1_dy
        Bi[2, 2] = dphi1_dy;  Bi[3, 2] = dphi1_dx
        Bi[1, 3] = dphi2_dx;  Bi[3, 3] = dphi2_dy
        Bi[2, 4] = dphi2_dy;  Bi[3, 4] = dphi2_dx

        K_ab_sr .+= (Bm_g' * Cm * Bi) .* detJ_g
        K_bb_sr .+= (Bi' * Cm * Bi) .* detJ_g
    end

    alpha = -(K_bb_sr \ (K_ab_sr' * u_elem))

    # Incompatible mode B-matrix at center (ξ=η=0)
    # At center: dphi1/dξ = 0 (since -2*0=0), dphi2/dη = 0
    # So Bi at center is zero! Need to add the incompatible strain contribution
    # via the mode amplitudes integrated over the element
    # Actually, the strain at center from incompatible modes:
    iJ_c = inv(J)
    # φ1 = 1-ξ², dφ1/dξ = -2ξ = 0 at center; φ2 = 1-η², dφ2/dη = -2η = 0 at center
    # So the incompatible mode derivatives are zero at center!
    # The strain correction from incompatible modes is zero at center.
    # BUT the forces N = ∫ σ dA are affected because the incompatible modes
    # change the strain field at the Gauss points.

    # For stress recovery, compute the average strain including incompatible modes
    # by integrating over Gauss points
    eps_mem_avg = zeros(3)
    kappa_avg = zeros(3)
    total_area = 0.0
    for i in 1:4
        r, s = gauss_pts[i,1], gauss_pts[i,2]
        dNr_g, dNs_g = shape_derivs_quad(r, s)
        J_g = [dNr_g'; dNs_g'] * coords
        detJ_g = abs(det(J_g))
        if detJ_g < 1e-12; detJ_g = 1e-12; end
        iJ = inv(J_g)
        dN_dxy_g = iJ * [dNr_g'; dNs_g']

        # Standard membrane strain at this GP
        Bm_g = zeros(3, 24)
        for k in 1:4
            idx = (k-1)*6
            Bm_g[1, idx+1] = dN_dxy_g[1,k]; Bm_g[2, idx+2] = dN_dxy_g[2,k]
            Bm_g[3, idx+1] = dN_dxy_g[2,k]; Bm_g[3, idx+2] = dN_dxy_g[1,k]
        end

        # Bending strain at this GP
        Bb_g = zeros(3, 24)
        for k in 1:4
            idx = (k-1)*6
            Bb_g[1, idx+5] = dN_dxy_g[1,k]
            Bb_g[2, idx+4] = -dN_dxy_g[2,k]
            Bb_g[3, idx+5] = dN_dxy_g[2,k]
            Bb_g[3, idx+4] = -dN_dxy_g[1,k]
        end

        # Incompatible mode strain at this GP
        dphi1_dx = iJ[1,1]*(-2.0*r); dphi1_dy = iJ[2,1]*(-2.0*r)
        dphi2_dx = iJ[1,2]*(-2.0*s); dphi2_dy = iJ[2,2]*(-2.0*s)

        Bi = zeros(3, 4)
        Bi[1, 1] = dphi1_dx;  Bi[3, 1] = dphi1_dy
        Bi[2, 2] = dphi1_dy;  Bi[3, 2] = dphi1_dx
        Bi[1, 3] = dphi2_dx;  Bi[3, 3] = dphi2_dy
        Bi[2, 4] = dphi2_dy;  Bi[3, 4] = dphi2_dx

        eps_gp = Bm_g * u_elem .+ Bi * alpha
        kappa_gp = Bb_g * u_elem

        eps_mem_avg .+= eps_gp .* detJ_g
        kappa_avg .+= kappa_gp .* detJ_g
        total_area += detJ_g
    end
    eps_mem_avg ./= total_area
    kappa_avg ./= total_area

    eps_mem = eps_mem_avg
    kappa = kappa_avg

    N = (D_mem * eps_mem) * h
    M = -bend_ratio * (D_mem * kappa) * (h^3/12.0)

    # Transverse shear forces from displacement-based shear strains at center
    G = E / (2*(1+nu))
    ts_t = 5.0/6.0
    N_center = [0.25, 0.25, 0.25, 0.25]
    gamma_xz = 0.0; gamma_yz = 0.0
    for k in 1:4
        idx = (k-1)*6
        gamma_xz += dN_dxy[1,k] * u_elem[idx+3] + N_center[k] * u_elem[idx+5]
        gamma_yz += dN_dxy[2,k] * u_elem[idx+3] - N_center[k] * u_elem[idx+4]
    end
    Q = [ts_t * G * h * gamma_xz, ts_t * G * h * gamma_yz]

    z1 = -h/2.0; z2 = h/2.0

    strain_z1 = eps_mem .+ z1 .* kappa
    stress_z1 = D_mem * strain_z1
    strain_z2 = eps_mem .+ z2 .* kappa
    stress_z2 = D_mem * strain_z2

    return N, M, Q, stress_z1, stress_z2, strain_z1, strain_z2
end

function stiffness_tria3(coords, E, nu, h; bend_ratio=1.0, ts_t=5.0/6.0, k6rot=100.0)
    G = E / (2*(1+nu))
    Dm = (E * h / (1 - nu^2)) .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    Db = bend_ratio * (E * h^3 / (12*(1 - nu^2))) .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    Ds = ts_t * G * h .* [1.0 0.0; 0.0 1.0]
    return stiffness_tria3_matrices(coords, Dm, Db, Ds, h, G; bend_ratio=bend_ratio, k6rot=k6rot)
end

# Overload accepting pre-computed constitutive matrices (for orthotropic MAT8)
function stiffness_tria3_matrices(coords, Dm, Db, Ds, h, G_ref; bend_ratio=1.0, ts_t=5.0/6.0, k6rot=100.0, Bmb=nothing)
    x, y = coords[:,1], coords[:,2]
    A2 = x[1]*(y[2]-y[3]) + x[2]*(y[3]-y[1]) + x[3]*(y[1]-y[2])
    A = 0.5 * abs(A2)
    if A < 1e-12; return zeros(18,18); end

    Ke = zeros(18, 18)

    # --- Membrane (constant strain triangle, CST) ---
    bv = [y[2]-y[3], y[3]-y[1], y[1]-y[2]] ./ (2*A)
    cv = [x[3]-x[2], x[1]-x[3], x[2]-x[1]] ./ (2*A)
    Bm = zeros(3, 6)
    for i in 1:3; Bm[1, i*2-1]=bv[i]; Bm[2, i*2]=cv[i]; Bm[3, i*2-1]=cv[i]; Bm[3, i*2]=bv[i]; end
    K_mem = Bm' * Dm * Bm * A
    m_idx = [1,2, 7,8, 13,14]
    Ke[m_idx, m_idx] = K_mem

    # --- Bending: Isoparametric Mindlin-Reissner Triangle ---
    Bb = zeros(3, 9)
    for i in 1:3
        col_rx = 3*(i-1) + 2  # θx
        col_ry = 3*(i-1) + 3  # θy
        Bb[1, col_ry] = bv[i]     # κxx
        Bb[2, col_rx] = -cv[i]    # κyy
        Bb[3, col_rx] = -bv[i]    # κxy: -∂θx/∂x
        Bb[3, col_ry] = cv[i]     # κxy: +∂θy/∂y
    end
    Kb = Bb' * Db * Bb * A

    # Transverse shear at centroid (1-point averaged, Tessler-Hughes approach)
    # γxz = ∂w/∂x + θy, γyz = ∂w/∂y - θx
    Bs = zeros(2, 9)
    for i in 1:3
        col_w  = 3*(i-1) + 1
        col_rx = 3*(i-1) + 2
        col_ry = 3*(i-1) + 3
        Bs[1, col_w]  = bv[i]
        Bs[1, col_ry] = 1.0/3.0    # N[i] = 1/3 at centroid
        Bs[2, col_w]  = cv[i]
        Bs[2, col_rx] = -1.0/3.0
    end
    Ks = Bs' * Ds * Bs * A

    # Adaptive phi² correction factor to prevent shear locking
    # Based on Tessler-Hughes MIN3 / NASTRAN CoFE approach
    rot_dofs = [2,3,5,6,8,9]  # θx,θy DOFs in 9-DOF system
    tr_Kb_rot = sum(Kb[i,i] for i in rot_dofs)
    tr_Ks_rot = sum(Ks[i,i] for i in rot_dofs)
    if tr_Kb_rot > 1e-30
        alpha = (1.0/ts_t) * tr_Ks_rot / tr_Kb_rot
        phi2 = 1.0 / (1.0 + 3.0*alpha)
    else
        phi2 = 0.0
    end

    b_idx = [3,4,5, 9,10,11, 15,16,17]
    Ke[b_idx, b_idx] += Kb + phi2 * Ks

    # B matrix coupling (membrane-bending): cross-blocks between m_idx and b_idx
    if Bmb !== nothing
        K_mb = Bm' * Bmb * Bb * A  # 6x9 coupling block
        Ke[m_idx, b_idx] += K_mb
        Ke[b_idx, m_idx] += K_mb'
    end

    # --- Hughes-Brezzi drilling rotation coupling ---
    # ε_drill = θz - (1/2)(∂v/∂x - ∂u/∂y), penalized with alpha_drill * ε_drill²
    alpha_drill = (k6rot / 1e5) * G_ref * h
    Bd = zeros(1, 18)
    for i in 1:3
        idx = (i-1)*6
        Bd[1, idx+1] = 0.5 * cv[i]    # +(1/2)*∂N/∂y (from ∂u/∂y)
        Bd[1, idx+2] = -0.5 * bv[i]   # -(1/2)*∂N/∂x (from -∂v/∂x)
        Bd[1, idx+6] = 1.0/3.0        # N_i = 1/3 at centroid
    end
    Ke .+= alpha_drill .* (Bd' * Bd) .* A

    return Ke
end

function stress_strain_tria3(coords, u_elem, E, nu, h; bend_ratio=1.0)
    x, y = coords[:,1], coords[:,2]
    A = 0.5 * abs(x[1]*(y[2]-y[3]) + x[2]*(y[3]-y[1]) + x[3]*(y[1]-y[2]))
    if A < 1e-12; return zeros(3), zeros(3), zeros(2), zeros(3), zeros(3), zeros(3), zeros(3); end

    b = [y[2]-y[3], y[3]-y[1], y[1]-y[2]] ./ (2*A)
    c = [x[3]-x[2], x[1]-x[3], x[2]-x[1]] ./ (2*A)

    # Membrane strain
    Bm = zeros(3, 6)
    for i in 1:3; Bm[1, i*2-1]=b[i]; Bm[2, i*2]=c[i]; Bm[3, i*2-1]=c[i]; Bm[3, i*2]=b[i]; end
    D = (E / (1 - nu^2)) .* [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    u_mem = [u_elem[1], u_elem[2], u_elem[7], u_elem[8], u_elem[13], u_elem[14]]
    eps_mem = Bm * u_mem

    # Bending curvature
    Bb = zeros(3, 6)
    for i in 1:3
        Bb[1, i*2]   = b[i]    # dθy/dx
        Bb[2, i*2-1] = -c[i]   # -dθx/dy
        Bb[3, i*2]   = c[i]    # dθy/dy
        Bb[3, i*2-1] = -b[i]   # -dθx/dx
    end
    u_rot = [u_elem[4], u_elem[5], u_elem[10], u_elem[11], u_elem[16], u_elem[17]]
    kappa = Bb * u_rot

    # Membrane forces and bending moments
    N = (D * eps_mem) * h
    M = -bend_ratio * (D * kappa) * (h^3/12.0)

    # Transverse shear forces from displacement-based shear strains at centroid
    G = E / (2*(1+nu))
    ts_t = 5.0/6.0
    N_ctr = 1.0/3.0
    gamma_xz = 0.0; gamma_yz = 0.0
    for i in 1:3
        w_i  = u_elem[(i-1)*6+3]
        rx_i = u_elem[(i-1)*6+4]
        ry_i = u_elem[(i-1)*6+5]
        gamma_xz += b[i] * w_i + N_ctr * ry_i
        gamma_yz += c[i] * w_i - N_ctr * rx_i
    end
    Q = [ts_t * G * h * gamma_xz, ts_t * G * h * gamma_yz]

    # Stresses at top/bottom surfaces
    z1 = -h/2.0; z2 = h/2.0
    strain_z1 = eps_mem .+ z1 .* kappa
    stress_z1 = D * strain_z1
    strain_z2 = eps_mem .+ z2 .* kappa
    stress_z2 = D * strain_z2

    return N, M, Q, stress_z1, stress_z2, strain_z1, strain_z2
end

end
