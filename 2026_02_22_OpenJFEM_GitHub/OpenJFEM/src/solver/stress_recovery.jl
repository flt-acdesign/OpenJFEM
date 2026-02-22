# stress_recovery.jl — Element stress/force/strain recovery

function recover_shell_stresses!(model, id_map, X, node_R, u_global, snorm_normals, stresses, results_json)
    lc_buf = zeros(4,2)

    for (id, el) in model["CSHELLs"]
        eid = parse(Int, id)
        pid = string(el["PID"])
        if !haskey(model["PSHELLs"], pid); continue; end
        prop = model["PSHELLs"][pid]
        mid = string(prop["MID"])
        if !haskey(model["MATs"], mid); continue; end
        mat = model["MATs"][mid]

        nids = el["NODES"]; n = length(nids)
        if any(x->get(id_map,x,0)==0, nids); continue; end

        local N, M, Q, s_z1, s_z2, e_z1, e_z2, elem_key
        stress_ok = true

        if n==4
            i1, i2, i3, i4 = id_map[nids[1]], id_map[nids[2]], id_map[nids[3]], id_map[nids[4]]
            p1 = SVector{3}(X[i1,1], X[i1,2], X[i1,3])
            p2 = SVector{3}(X[i2,1], X[i2,2], X[i2,3])
            p3 = SVector{3}(X[i3,1], X[i3,2], X[i3,3])
            p4 = SVector{3}(X[i4,1], X[i4,2], X[i4,3])
            sr_indices = [i1, i2, i3, i4]

            v1, v2, v3 = shell_element_frame_fast(p1, p2, p3, p4, 4)
            v1, v2, v3 = apply_snorm_to_frame(v1, v2, v3, sr_indices, snorm_normals)
            c = (p1+p2+p3+p4)/4.0
            lc_buf[1,1]=dot(p1-c,v1); lc_buf[1,2]=dot(p1-c,v2)
            lc_buf[2,1]=dot(p2-c,v1); lc_buf[2,2]=dot(p2-c,v2)
            lc_buf[3,1]=dot(p3-c,v1); lc_buf[3,2]=dot(p3-c,v2)
            lc_buf[4,1]=dot(p4-c,v1); lc_buf[4,2]=dot(p4-c,v2)

            Rel_t = vcat(v1', v2', v3')
            u_el = zeros(24)
            for k=1:4
                idx = id_map[nids[k]]
                u_el[(k-1)*6+1:(k-1)*6+3] = Rel_t * node_R[idx] * u_global[(idx-1)*6+1:(idx-1)*6+3]
                u_el[(k-1)*6+4:(k-1)*6+6] = Rel_t * node_R[idx] * u_global[(idx-1)*6+4:(idx-1)*6+6]
            end
            br = get(prop, "BEND_RATIO", 1.0)
            clt_Cm = (get(prop, "TYPE", "") == "PCOMP_CLT" && haskey(prop, "Cm")) ? prop["Cm"] : nothing
            try
                N, M, Q, s_z1, s_z2, e_z1, e_z2 = FEM.stress_strain_quad4(view(lc_buf,1:4,:), u_el, mat["E"], mat["NU"], Float64(prop["T"]), Float64(prop["T"]); bend_ratio=br, Cm_override=clt_Cm)
            catch e
                @warn "Stress recovery failed for QUAD4 $eid: $e"
                stress_ok = false
            end
            elem_key = "quad4"
        elseif n==3
            i1, i2, i3 = id_map[nids[1]], id_map[nids[2]], id_map[nids[3]]
            p1 = SVector{3}(X[i1,1], X[i1,2], X[i1,3])
            p2 = SVector{3}(X[i2,1], X[i2,2], X[i2,3])
            p3 = SVector{3}(X[i3,1], X[i3,2], X[i3,3])
            p4=SVector(0.0,0.0,0.0)
            sr_indices = [i1, i2, i3]

            v1, v2, v3 = shell_element_frame_fast(p1, p2, p3, p4, 3)
            v1, v2, v3 = apply_snorm_to_frame(v1, v2, v3, sr_indices, snorm_normals)
            c = (p1+p2+p3)/3.0
            lc_buf[1,1]=dot(p1-c,v1); lc_buf[1,2]=dot(p1-c,v2)
            lc_buf[2,1]=dot(p2-c,v1); lc_buf[2,2]=dot(p2-c,v2)
            lc_buf[3,1]=dot(p3-c,v1); lc_buf[3,2]=dot(p3-c,v2)
            Rel_t = vcat(v1', v2', v3')
            u_el = zeros(18)
            for k=1:3
                idx = id_map[nids[k]]
                u_el[(k-1)*6+1:(k-1)*6+3] = Rel_t * node_R[idx] * u_global[(idx-1)*6+1:(idx-1)*6+3]
                u_el[(k-1)*6+4:(k-1)*6+6] = Rel_t * node_R[idx] * u_global[(idx-1)*6+4:(idx-1)*6+6]
            end
            br = get(prop, "BEND_RATIO", 1.0)
            try
                N, M, Q, s_z1, s_z2, e_z1, e_z2 = FEM.stress_strain_tria3(view(lc_buf,1:3,:), u_el, mat["E"], mat["NU"], Float64(prop["T"]); bend_ratio=br)
            catch e
                @warn "Stress recovery failed for TRIA3 $eid: $e"
                stress_ok = false
            end
            elem_key = "tria3"
        else
            continue
        end

        if !stress_ok; continue; end

        eps_mem_out = (e_z1 .+ e_z2) ./ 2.0
        kappa_nast_out = (e_z1 .- e_z2) ./ prop["T"]

        is_pcomp = get(prop, "TYPE", "") == "PCOMP_CLT" && haskey(prop, "PLY_DATA")
        if is_pcomp
            t_total = prop["T"]
            eps_mem = (e_z1 .+ e_z2) ./ 2.0
            kappa = (e_z2 .- e_z1) ./ t_total
            N = prop["Cm"] * eps_mem
            M = -prop["Cb"] * kappa

            ply_data = prop["PLY_DATA"]
            vm_max = 0.0
            s_z1_out = zeros(3)
            s_z2_out = zeros(3)
            e_z1_out = zeros(3)
            e_z2_out = zeros(3)
            for (ip, pd) in enumerate(ply_data)
                Qbar = pd["Qbar"]
                z_mid = (pd["z_bot"] + pd["z_top"]) / 2.0
                strain_ply = eps_mem .+ z_mid .* kappa
                stress_ply = Qbar * strain_ply
                vm_ply = sqrt(stress_ply[1]^2 - stress_ply[1]*stress_ply[2] + stress_ply[2]^2 + 3*stress_ply[3]^2)
                if vm_ply > vm_max; vm_max = vm_ply; end
                if ip == 1; s_z1_out .= stress_ply; e_z1_out .= strain_ply; end
                if ip == length(ply_data); s_z2_out .= stress_ply; e_z2_out .= strain_ply; end
            end
            s_z1 = s_z1_out
            s_z2 = s_z2_out
            e_z1 = e_z1_out
            e_z2 = e_z2_out
            stresses[eid] = vm_max
        else
            stresses[eid] = FEM.compute_principal_2d(s_z1[1], s_z1[2], s_z1[3])[1]
        end

        push!(results_json["forces"][elem_key], Dict("eid" => eid, "fx" => N[1], "fy" => N[2], "fxy" => N[3], "mx" => M[1], "my" => M[2], "mxy" => M[3], "qx" => Q[1], "qy" => Q[2]))

        make_stress_entry(s, t) = Dict("fiber_dist" => t, "normal_x" => s[1], "normal_y" => s[2], "shear_xy" => s[3], "von_mises" => sqrt(s[1]^2-s[1]*s[2]+s[2]^2+3*s[3]^2), "major" => 0.0, "minor" => 0.0)
        make_strain_entry(e, t) = Dict("fiber_dist" => t, "normal_x" => e[1], "normal_y" => e[2], "shear_xy" => e[3], "major" => 0.0, "minor" => 0.0)

        push!(results_json["stresses"][elem_key], Dict("eid" => eid, "z1" => make_stress_entry(s_z1, -prop["T"]/2), "z2" => make_stress_entry(s_z2, prop["T"]/2)))
        push!(results_json["strains"][elem_key], Dict("eid" => eid, "z1" => make_strain_entry(eps_mem_out, 0.0), "z2" => make_strain_entry(kappa_nast_out, -1.0)))
    end
end

function recover_bar_stresses!(model, id_map, X, node_R, u_global, stresses, results_json)
    for (id, bar) in model["CBARs"]
        eid = parse(Int, id)
        pid = string(bar["PID"])
        if !haskey(model["PBARLs"], pid); continue; end
        prop = model["PBARLs"][pid]
        mid = string(prop["MID"])
        if !haskey(model["MATs"], mid); continue; end
        mat = model["MATs"][mid]

        if !haskey(id_map, bar["GA"]) || !haskey(id_map, bar["GB"]); continue; end
        i1, i2 = id_map[bar["GA"]], id_map[bar["GB"]]

        p1 = SVector{3}(X[i1,1], X[i1,2], X[i1,3])
        p2 = SVector{3}(X[i2,1], X[i2,2], X[i2,3])

        wa, wb, has_offset, p1_eff, p2_eff = bar_offsets_and_endpoints(bar, p1, p2)
        L = norm(p2_eff - p1_eff)
        if L < 1e-9; continue; end
        vx = normalize(p2_eff - p1_eff)
        v_ref = resolve_bar_vref(bar, p1, id_map, X)

        if norm(v_ref) < 1e-6
             v_ref = SVector(0.0,0.0,1.0)
             if abs(dot(vx, v_ref)) > 0.9; v_ref = SVector(0.0,1.0,0.0); end
        end
        vz = normalize(cross(vx, v_ref))
        vy = cross(vz, vx)
        Rel_t = vcat(vx', vy', vz')

        TR1 = Rel_t * node_R[i1]
        TR2 = Rel_t * node_R[i2]
        u_el = zeros(12)
        u_el[1:3] = TR1 * u_global[(i1-1)*6+1:(i1-1)*6+3]
        u_el[4:6] = TR1 * u_global[(i1-1)*6+4:(i1-1)*6+6]
        u_el[7:9] = TR2 * u_global[(i2-1)*6+1:(i2-1)*6+3]
        u_el[10:12] = TR2 * u_global[(i2-1)*6+4:(i2-1)*6+6]
        if has_offset
            S_wa = skew3(wa); S_wb = skew3(wb)
            θ_glob_A = node_R[i1] * u_global[(i1-1)*6+4:(i1-1)*6+6]
            θ_glob_B = node_R[i2] * u_global[(i2-1)*6+4:(i2-1)*6+6]
            u_el[1:3] -= Rel_t * S_wa * θ_glob_A
            u_el[7:9] -= Rel_t * S_wb * θ_glob_B
        end

        Iy = get(prop, "I2", get(prop, "I", 0.0))
        Iz = get(prop, "I1", get(prop, "I", 0.0))
        if Iy == 0.0; Iy = Iz; end
        if Iz == 0.0; Iz = Iy; end
        Iyz = Float64(get(prop, "I12", 0.0))
        K1 = get(prop, "K1", 0.0); K2 = get(prop, "K2", 0.0)
        As_y = (K1 > 0.0) ? K1 * prop["A"] : Inf
        As_z = (K2 > 0.0) ? K2 * prop["A"] : Inf
        forces = FEM.forces_frame3d(u_el, L, prop["A"], Iy, Iz, prop["J"], mat["E"], mat["G"]; As_y=As_y, As_z=As_z, I12=Iyz)
        A_bar = Float64(prop["A"])
        sig_axial = (abs(A_bar) > 1e-30) ? forces["axial"]/A_bar : 0.0
        stresses[eid] = abs(sig_axial)
        push!(results_json["forces"]["cbar"], Dict("eid" => eid, "axial" => forces["axial"], "shear_1" => forces["shear_1"], "shear_2" => forces["shear_2"], "torque" => forces["torque"], "moment_a1" => forces["moment_a1"], "moment_a2" => forces["moment_a2"], "moment_b1" => forces["moment_b1"], "moment_b2" => forces["moment_b2"]))

        I12_sr = Float64(get(prop, "I12", 0.0))
        det_I = Iz*Iy - I12_sr^2
        sr_pts = [(get(prop,"C1",0.0), get(prop,"C2",0.0)),
                  (get(prop,"D1",0.0), get(prop,"D2",0.0)),
                  (get(prop,"E1",0.0), get(prop,"E2",0.0)),
                  (get(prop,"F1",0.0), get(prop,"F2",0.0))]
        end_a = Dict{String,Float64}()
        end_b = Dict{String,Float64}()
        for (j, (yj, zj)) in enumerate(sr_pts)
            if abs(det_I) > 1e-30
                end_a["p$j"] = sig_axial - (forces["moment_a1"]*Iy - forces["moment_a2"]*I12_sr)*yj/det_I - (forces["moment_a2"]*Iz - forces["moment_a1"]*I12_sr)*zj/det_I
                end_b["p$j"] = sig_axial - (forces["moment_b1"]*Iy - forces["moment_b2"]*I12_sr)*yj/det_I - (forces["moment_b2"]*Iz - forces["moment_b1"]*I12_sr)*zj/det_I
            else
                end_a["p$j"] = sig_axial - forces["moment_a1"]*yj/max(Iz,1e-30) - forces["moment_a2"]*zj/max(Iy,1e-30)
                end_b["p$j"] = sig_axial - forces["moment_b1"]*yj/max(Iz,1e-30) - forces["moment_b2"]*zj/max(Iy,1e-30)
            end
        end
        push!(results_json["stresses"]["cbar"], Dict("eid"=>eid, "end_a"=>end_a, "end_b"=>end_b, "axial"=>sig_axial))
    end
end

function recover_rod_stresses!(model, id_map, X, node_R, u_global, stresses, results_json)
    # CROD recovery
    crods = get(model, "CRODs", Dict())
    prods = get(model, "PRODs", Dict())
    for (id, rod) in crods
        eid = parse(Int, id)
        pid = string(rod["PID"])
        if !haskey(prods, pid); continue; end
        prop = prods[pid]
        mid = string(prop["MID"])
        if !haskey(model["MATs"], mid); continue; end
        mat = model["MATs"][mid]

        if !haskey(id_map, rod["GA"]) || !haskey(id_map, rod["GB"]); continue; end
        i1, i2 = id_map[rod["GA"]], id_map[rod["GB"]]

        p1 = SVector{3}(X[i1,1], X[i1,2], X[i1,3])
        p2 = SVector{3}(X[i2,1], X[i2,2], X[i2,3])
        L = norm(p2-p1)
        if L < 1e-9; continue; end
        vx = normalize(p2-p1)
        ref = abs(vx[3]) < 0.9 ? SVector(0.0,0.0,1.0) : SVector(0.0,1.0,0.0)
        vz = normalize(cross(vx, ref))
        vy = cross(vz, vx)
        Rel_t = vcat(vx', vy', vz')

        u_el = zeros(12)
        u_el[1:3] = Rel_t * node_R[i1] * u_global[(i1-1)*6+1:(i1-1)*6+3]
        u_el[4:6] = Rel_t * node_R[i1] * u_global[(i1-1)*6+4:(i1-1)*6+6]
        u_el[7:9] = Rel_t * node_R[i2] * u_global[(i2-1)*6+1:(i2-1)*6+3]
        u_el[10:12] = Rel_t * node_R[i2] * u_global[(i2-1)*6+4:(i2-1)*6+6]

        axial_force = mat["E"] * prop["A"] / L * (u_el[7] - u_el[1])
        torque = mat["G"] * prop["J"] / L * (u_el[10] - u_el[4])
        axial_stress = prop["A"] > 0 ? axial_force / prop["A"] : 0.0
        torsional_stress = prop["J"] > 0 && haskey(prop, "C") && prop["C"] > 0 ? torque * prop["C"] / prop["J"] : 0.0
        stresses[eid] = abs(axial_stress)
        axial_strain = mat["E"] > 0 ? axial_stress / mat["E"] : 0.0
        push!(results_json["forces"]["crod"], Dict("eid" => eid, "axial" => axial_force, "torque" => torque))
        push!(results_json["stresses"]["crod"], Dict("eid" => eid, "axial" => axial_stress, "torsional" => torsional_stress))
        push!(results_json["strains"]["crod"], Dict("eid" => eid, "axial" => axial_strain, "torsional" => mat["G"] > 0 ? torsional_stress / mat["G"] : 0.0))
    end

    # CONROD recovery
    conrods = get(model, "CONRODs", Dict())
    for (id, rod) in conrods
        eid = parse(Int, id)
        mid = string(rod["MID"])
        if !haskey(model["MATs"], mid); continue; end
        mat = model["MATs"][mid]
        if !haskey(id_map, rod["GA"]) || !haskey(id_map, rod["GB"]); continue; end
        i1, i2 = id_map[rod["GA"]], id_map[rod["GB"]]
        p1 = SVector{3}(X[i1,1], X[i1,2], X[i1,3])
        p2 = SVector{3}(X[i2,1], X[i2,2], X[i2,3])
        L = norm(p2-p1)
        if L < 1e-9; continue; end
        vx = normalize(p2-p1)
        ref = abs(vx[3]) < 0.9 ? SVector(0.0,0.0,1.0) : SVector(0.0,1.0,0.0)
        vz = normalize(cross(vx, ref))
        vy = cross(vz, vx)
        Rel_t = vcat(vx', vy', vz')
        u_el = zeros(12)
        u_el[1:3] = Rel_t * node_R[i1] * u_global[(i1-1)*6+1:(i1-1)*6+3]
        u_el[4:6] = Rel_t * node_R[i1] * u_global[(i1-1)*6+4:(i1-1)*6+6]
        u_el[7:9] = Rel_t * node_R[i2] * u_global[(i2-1)*6+1:(i2-1)*6+3]
        u_el[10:12] = Rel_t * node_R[i2] * u_global[(i2-1)*6+4:(i2-1)*6+6]
        axial_force = mat["E"] * rod["A"] / L * (u_el[7] - u_el[1])
        torque = mat["G"] * rod["J"] / L * (u_el[10] - u_el[4])
        axial_stress = rod["A"] > 0 ? axial_force / rod["A"] : 0.0
        torsional_stress = rod["J"] > 0 && rod["C"] > 0 ? torque * rod["C"] / rod["J"] : 0.0
        stresses[eid] = abs(axial_stress)
        axial_strain = mat["E"] > 0 ? axial_stress / mat["E"] : 0.0
        push!(results_json["forces"]["conrod"], Dict("eid" => eid, "axial" => axial_force, "torque" => torque))
        push!(results_json["stresses"]["conrod"], Dict("eid" => eid, "axial" => axial_stress, "torsional" => torsional_stress))
        push!(results_json["strains"]["conrod"], Dict("eid" => eid, "axial" => axial_strain, "torsional" => mat["G"] > 0 ? torsional_stress / mat["G"] : 0.0))
    end
end

function recover_spring_forces!(model, id_map, u_global, stresses, results_json)
    celases = get(model, "CELASs", Dict())
    pelases = get(model, "PELASs", Dict())
    for (id, spring) in celases
        eid = parse(Int, id)
        pid = string(spring["PID"])
        if !haskey(pelases, pid); continue; end
        K_spring = pelases[pid]["K"]
        g1 = spring["G1"]; c1 = spring["C1"]
        g2 = spring["G2"]; c2 = spring["C2"]
        u1 = 0.0; u2 = 0.0
        if g1 > 0 && haskey(id_map, g1) && c1 > 0
            u1 = u_global[(id_map[g1]-1)*6 + c1]
        end
        if g2 > 0 && haskey(id_map, g2) && c2 > 0
            u2 = u_global[(id_map[g2]-1)*6 + c2]
        end
        spring_force = K_spring * (u1 - u2)
        stresses[eid] = abs(spring_force)
        push!(results_json["forces"]["celas1"], Dict("eid" => eid, "force" => spring_force))
        push!(results_json["stresses"]["celas1"], Dict("eid" => eid, "force" => spring_force))
        push!(results_json["strains"]["celas1"], Dict("eid" => eid, "deformation" => u1 - u2))
    end
end
