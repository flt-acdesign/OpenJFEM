# ModelBuilder.jl
# Functions for building the finite element model from parsed Nastran cards.
# This file is included at the top level and has access to the NastranParser module.
#
# Contains:
#   transform_geometry!(model) - Transforms grid coordinates from local coordinate systems
#   build_model(cards, cc)     - Constructs the full model dict from parsed cards and case control

function transform_geometry!(model)
    grids = model["GRIDs"]
    cords = model["CORDs"]

    for (sid, g) in grids
        if g["CP"] != 0 && haskey(cords, string(g["CP"]))
            c = cords[string(g["CP"])]
            R = hcat(c["U"], c["V"], c["W"])
            ctype = get(c, "TYPE", "RECTANGULAR")

            if ctype == "CYLINDRICAL"
                r, theta_deg, z = g["X"][1], g["X"][2], g["X"][3]
                theta = deg2rad(theta_deg)
                local_xyz = [r * cos(theta), r * sin(theta), z]
                g["X"] = c["Origin"] + R * local_xyz
            elseif ctype == "SPHERICAL"
                r, theta_deg, phi_deg = g["X"][1], g["X"][2], g["X"][3]
                theta = deg2rad(theta_deg)
                phi = deg2rad(phi_deg)
                local_xyz = [r * sin(phi) * cos(theta), r * sin(phi) * sin(theta), r * cos(phi)]
                g["X"] = c["Origin"] + R * local_xyz
            else
                g["X"] = c["Origin"] + R * g["X"]
            end
        end
    end
end

function build_model(cards, cc)
    println(">>> Constructing Model Data...")

    # Parse PBARL and PBAR, then merge into unified PBARLs dict
    pbarls = haskey(cards,"PBARL") ? NastranParser.extract_pbarl(cards["PBARL"]) : Dict()
    pbars  = haskey(cards,"PBAR")  ? NastranParser.extract_pbar(cards["PBAR"])   : Dict()
    # Also handle PBAR* (long-field PBAR) - process_cards stores these under "PBAR*"
    if haskey(cards, "PBAR*")
        pbars_long = NastranParser.extract_pbar(cards["PBAR*"])
        merge!(pbars, pbars_long)
    end
    # Merge PBAR into PBARLs (PBAR takes precedence if same PID)
    merged_bars = merge(pbarls, pbars)

    # Parse PLOAD2 and merge with PLOAD4
    pload4s = haskey(cards,"PLOAD4") ? NastranParser.extract_pload4(cards["PLOAD4"]) : []
    if haskey(cards, "PLOAD2")
        pload2s = NastranParser.extract_pload2(cards["PLOAD2"])
        append!(pload4s, pload2s)
    end

    # Parse SPC (non-SPC1) and merge with SPC1
    spc1s = haskey(cards,"SPC1") ? NastranParser.extract_spc1(cards["SPC1"]) : []
    if haskey(cards, "SPC")
        spcs = NastranParser.extract_spc(cards["SPC"])
        append!(spc1s, spcs)
    end

    # Parse PROD properties
    prods = haskey(cards,"PROD") ? NastranParser.extract_prod(cards["PROD"]) : Dict()

    # Parse MAT8, MAT2 and merge with MAT1
    mats = haskey(cards,"MAT1") ? NastranParser.extract_mats(cards["MAT1"]) : Dict()
    if haskey(cards, "MAT8")
        mat8s = NastranParser.extract_mat8(cards["MAT8"])
        merge!(mats, mat8s)
    end
    if haskey(cards, "MAT2")
        mat2s = NastranParser.extract_mat2(cards["MAT2"])
        merge!(mats, mat2s)
    end

    # Parse PELAS properties
    pelases = haskey(cards,"PELAS") ? NastranParser.extract_pelas(cards["PELAS"]) : Dict()

    # Parse PCOMP and merge with PSHELL
    # PCOMP uses proper CLT (Classical Laminate Theory) to compute full ABD matrices
    pshells = haskey(cards,"PSHELL") ? NastranParser.extract_props_shell(cards["PSHELL"]) : Dict()
    if haskey(cards, "PCOMP")
        pcomps = NastranParser.extract_pcomp(cards["PCOMP"])
        for (pid, pc) in pcomps
            if isempty(pc["PLIES"]); continue; end
            total_t = pc["T"]
            if total_t <= 0; continue; end

            # CLT: compute full ABD matrices with z-coordinates
            A = zeros(3, 3)  # membrane stiffness
            B = zeros(3, 3)  # membrane-bending coupling
            D = zeros(3, 3)  # bending stiffness
            z0 = get(pc, "Z0", -total_t / 2.0)  # bottom of laminate
            z_bot = z0
            G12_ref = 0.0  # for transverse shear
            E_max = 0.0  # for drill stiffness reference

            ply_data = []  # store per-ply Qbar and z for stress recovery
            for ply in pc["PLIES"]
                pmid = string(ply["MID"])
                if !haskey(mats, pmid); continue; end
                pm = mats[pmid]
                t = ply["T"]; theta = deg2rad(ply["THETA"])
                z_top = z_bot + t

                if haskey(pm, "E1")  # MAT8
                    E1 = pm["E1"]; E2 = pm["E2"]; nu12 = pm["NU12"]; G12 = pm["G12"]
                else  # MAT1 (isotropic)
                    E1 = pm["E"]; E2 = pm["E"]; nu12 = pm["NU"]; G12 = pm["G"]
                end
                if G12 > G12_ref; G12_ref = G12; end
                if max(E1, E2) > E_max; E_max = max(E1, E2); end

                nu21 = nu12 * E2 / max(E1, 1e-30)
                denom = 1.0 - nu12 * nu21
                Q11 = E1/denom; Q22 = E2/denom; Q12 = nu12*E2/denom; Q66 = G12

                # Full rotated Q-bar (including Q16, Q26)
                c = cos(theta); s = sin(theta)
                c2 = c^2; s2 = s^2; cs = c*s
                Qb = zeros(3, 3)
                Qb[1,1] = Q11*c2^2 + 2*(Q12+2*Q66)*c2*s2 + Q22*s2^2
                Qb[2,2] = Q11*s2^2 + 2*(Q12+2*Q66)*c2*s2 + Q22*c2^2
                Qb[1,2] = (Q11+Q22-4*Q66)*c2*s2 + Q12*(c2^2+s2^2)
                Qb[2,1] = Qb[1,2]
                Qb[1,3] = (Q11-Q12-2*Q66)*c*s*c2 + (Q12-Q22+2*Q66)*c*s*s2
                Qb[3,1] = Qb[1,3]
                Qb[2,3] = (Q11-Q12-2*Q66)*c*s*s2 + (Q12-Q22+2*Q66)*c*s*c2
                Qb[3,2] = Qb[2,3]
                Qb[3,3] = (Q11+Q22-2*Q12-2*Q66)*c2*s2 + Q66*(c2^2+s2^2)

                # CLT integration: A += Qb*(z_top-z_bot), B += Qb*(z_top^2-z_bot^2)/2, D += Qb*(z_top^3-z_bot^3)/3
                A .+= Qb .* (z_top - z_bot)
                B .+= Qb .* (z_top^2 - z_bot^2) / 2.0
                D .+= Qb .* (z_top^3 - z_bot^3) / 3.0

                # Store ply data for stress recovery
                push!(ply_data, Dict("Qbar"=>copy(Qb), "z_bot"=>z_bot, "z_top"=>z_top,
                                     "theta"=>ply["THETA"], "sout"=>get(ply, "SOUT", "")))

                z_bot = z_top
            end

            # Transverse shear: use weighted average of ply G values with shear correction 5/6
            ts_t = 5.0/6.0
            Cs_lam = ts_t * total_t .* [G12_ref 0.0; 0.0 G12_ref]

            # Derive equivalent E and nu for stress recovery and drill stiffness
            nu_eq = A[1,1] > 0 ? clamp(A[1,2] / A[1,1], 0.0, 0.49) : 0.3
            E_eq = A[1,1] * (1 - nu_eq^2) / total_t
            G_eq = G12_ref

            pid_int = pc["PID"]
            synth_mid = 900000 + pid_int
            mats[string(synth_mid)] = Dict("MID"=>synth_mid, "E"=>E_eq, "G"=>G_eq, "NU"=>nu_eq, "TYPE"=>"MAT1_EQUIV")
            # Check if B is effectively zero (symmetric laminate)
            B_max = maximum(abs.(B))
            Bmb = B_max > 1e-10 * maximum(abs.(A)) ? B : nothing

            pshells[pid] = Dict("PID"=>pid_int, "MID"=>synth_mid, "T"=>total_t,
                                "TYPE"=>"PCOMP_CLT",
                                "Cm" => A, "Bmb" => Bmb, "Cb" => D, "Cs" => Cs_lam, "E_ref" => E_max,
                                "PLY_DATA" => ply_data)
        end
    end

    model = Dict(
        "CASE_CONTROL" => cc,
        "GRIDs"       => haskey(cards,"GRID")   ? NastranParser.extract_grid(cards["GRID"]) : Dict(),
        "CORDs"       => merge(haskey(cards,"CORD2R") ? NastranParser.extract_coords(cards["CORD2R"]; coord_type="RECTANGULAR") : Dict(),
                               haskey(cards,"CORD1R") ? NastranParser.extract_coords(cards["CORD1R"]; coord_type="RECTANGULAR") : Dict(),
                               haskey(cards,"CORD2C") ? NastranParser.extract_coords(cards["CORD2C"]; coord_type="CYLINDRICAL") : Dict(),
                               haskey(cards,"CORD2S") ? NastranParser.extract_coords(cards["CORD2S"]; coord_type="SPHERICAL") : Dict()),
        "CSHELLs"     => merge(haskey(cards,"CTRIA3") ? NastranParser.extract_shells(cards["CTRIA3"]) : Dict(),
                               haskey(cards,"CQUAD4") ? NastranParser.extract_shells(cards["CQUAD4"]) : Dict()),
        "CBARs"       => haskey(cards,"CBAR")   ? NastranParser.extract_cbar(cards["CBAR"]) : Dict(),
        "CRODs"       => haskey(cards,"CROD")   ? NastranParser.extract_crod(cards["CROD"]) : Dict(),
        "CONRODs"     => haskey(cards,"CONROD") ? NastranParser.extract_conrod(cards["CONROD"]) : Dict(),
        "CELASs"      => haskey(cards,"CELAS1") ? NastranParser.extract_celas1(cards["CELAS1"]) : Dict(),
        "PELASs"      => pelases,
        "CONM2s"      => haskey(cards,"CONM2")  ? NastranParser.extract_conm2(cards["CONM2"]) : Dict(),
        "RBE2s"       => haskey(cards,"RBE2")   ? NastranParser.extract_rbe2(cards["RBE2"]) : Dict(),
        "RBE3s"       => haskey(cards,"RBE3")   ? NastranParser.extract_rbe3(cards["RBE3"]) : Dict(),
        "MPCs"        => haskey(cards,"MPC")    ? NastranParser.extract_mpc(cards["MPC"]) : [],
        "MPCADDs"     => haskey(cards,"MPCADD") ? NastranParser.extract_mpcadd(cards["MPCADD"]) : Dict(),
        "PSHELLs"     => pshells,
        "PBARLs"      => merged_bars,
        "PRODs"       => prods,
        "MATs"        => mats,
        "FORCEs"      => haskey(cards,"FORCE")  ? NastranParser.extract_loads(cards["FORCE"]) : [],
        "MOMENTs"     => haskey(cards,"MOMENT") ? NastranParser.extract_moments(cards["MOMENT"]) : [],
        "PLOAD4s"     => pload4s,
        "PLOAD1s"     => haskey(cards,"PLOAD1") ? NastranParser.extract_pload1(cards["PLOAD1"]) : [],
        "GRAVs"       => haskey(cards,"GRAV")   ? NastranParser.extract_grav(cards["GRAV"]) : [],
        "LOAD_COMBOS" => haskey(cards,"LOAD")   ? NastranParser.extract_load_combos(cards["LOAD"]) : [],
        "SPC1s"       => spc1s,
        "SPCADDs"     => haskey(cards,"SPCADD") ? NastranParser.extract_spcadd(cards["SPCADD"]) : Dict()
    )

    # Extract PARAM cards
    if haskey(cards, "PARAM")
        for c in cards["PARAM"]
            pname = uppercase(strip(string(NastranParser.safe_get(c, 3))))
            pval = NastranParser.parse_nastran_number(NastranParser.safe_get(c, 4), 0.0)
            model["PARAM_$pname"] = pval
        end
    end

    return model
end
