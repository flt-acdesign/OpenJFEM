# extract_properties.jl — PSHELL, PBARL, PBAR, PCOMP, PROD, PELAS

function extract_props_shell(cards)
    d = Dict()
    for c in cards
        pid = to_id(parse_nastran_number(safe_get(c, 3)))
        mid = to_id(parse_nastran_number(safe_get(c, 4)))
        t = parse_nastran_number(safe_get(c, 5), 0.0)
        # MID2 = field 5 (c[6]), 12I/T^3 = field 6 (c[7]), MID3 = field 7 (c[8]), TS/T = field 8 (c[9])
        bend_ratio = parse_nastran_number(safe_get(c, 7), 1.0)  # 12I/T^3, default=1.0
        ts_t = parse_nastran_number(safe_get(c, 9), 5.0/6.0)    # TS/T, default=5/6
        d[string(pid)] = Dict("PID"=>pid, "MID"=>mid, "T"=>t, "BEND_RATIO"=>bend_ratio, "TS_T"=>ts_t)
    end
    return d
end

function extract_pbarl(cards)
    d = Dict()
    for c in cards
        pid = to_id(parse_nastran_number(safe_get(c, 3)))
        mid = to_id(parse_nastran_number(safe_get(c, 4)))
        type = "ROD"
        dim_start_idx = 7
        for k in 5:min(12, length(c))
            val = strip(string(safe_get(c, k, "")))
            if !isempty(val) && occursin(r"^[A-Za-z]", val)
                type = uppercase(val); dim_start_idx = k + 1; break
            end
        end
        # Collect all dimension values
        dims = Float64[]
        for k in dim_start_idx:length(c)
            val = parse_nastran_number(safe_get(c, k), nothing)
            if isa(val, Number)
                push!(dims, Float64(val))
            end
        end

        A, I1, I2, J = 1.0, 1.0, 1.0, 1.0
        K1, K2 = 0.0, 0.0
        nu_approx = 0.3  # approximate NU for shear factor computation

        if type == "ROD" && length(dims) >= 1
            R = dims[1]; A = pi*R^2; I1 = pi*R^4/4; I2 = I1; J = pi*R^4/2
            K1 = 6*(1+nu_approx)/(7+6*nu_approx); K2 = K1
        elseif (type == "TUBE" || type == "TUBE2") && length(dims) >= 2
            R_out = dims[1]; R_in = dims[2]
            if R_in < 0; R_in = 0.0; end
            A = pi*(R_out^2 - R_in^2)
            I1 = pi*(R_out^4 - R_in^4)/4
            I2 = I1
            J = pi*(R_out^4 - R_in^4)/2
            # Cowper formula for hollow circular section
            m = (R_out > 0) ? R_in/R_out : 0.0
            m2 = m^2
            K1 = 6*(1+nu_approx)*(1+m2)^2 / ((7+6*nu_approx)*(1+m2)^2 + (20+12*nu_approx)*m2)
            K2 = K1
        elseif type == "BAR" && length(dims) >= 2
            w = dims[1]; h = dims[2]
            A = w * h
            # NASTRAN BAR: DIM1=z-extent, DIM2=y-extent
            # I1 = Iz (Plane 1) = DIM1*DIM2^3/12;  I2 = Iy (Plane 2) = DIM2*DIM1^3/12
            I1 = w * h^3 / 12
            I2 = h * w^3 / 12
            if w >= h && w > 0
                J = w * h^3 / 3 * (1 - 0.63 * h / w)
            elseif h > 0
                J = h * w^3 / 3 * (1 - 0.63 * w / h)
            else
                J = I1 + I2
            end
            K1 = 10*(1+nu_approx)/(12+11*nu_approx); K2 = K1
        elseif type == "BOX" && length(dims) >= 4
            w = dims[1]; h = dims[2]; tw = dims[3]; th = dims[4]
            w_in = w - 2*tw; h_in = h - 2*th
            if w_in < 0; w_in = 0; end
            if h_in < 0; h_in = 0; end
            A = w*h - w_in*h_in
            I1 = (w*h^3 - w_in*h_in^3)/12
            I2 = (h*w^3 - h_in*w_in^3)/12
            J = 2*tw*th*(w-tw)^2*(h-th)^2 / (tw*(w-tw) + th*(h-th))
            K1 = 5.0/6.0; K2 = 5.0/6.0
        elseif type == "I" && length(dims) >= 6
            H = dims[1]; Bb = dims[2]; Bt = dims[3]; tw = dims[4]; tfb = dims[5]; tft = dims[6]
            hw = H - tfb - tft
            A = Bb*tfb + hw*tw + Bt*tft
            yb = (Bb*tfb*tfb/2 + hw*tw*(tfb+hw/2) + Bt*tft*(H-tft/2)) / A
            I1 = Bb*tfb^3/12 + Bb*tfb*(yb-tfb/2)^2 + tw*hw^3/12 + tw*hw*(yb-tfb-hw/2)^2 + Bt*tft^3/12 + Bt*tft*(yb-(H-tft/2))^2
            I2 = tfb*Bb^3/12 + hw*tw^3/12 + tft*Bt^3/12
            J = (Bb*tfb^3 + hw*tw^3 + Bt*tft^3) / 3
            K1 = hw*tw/A; K2 = 5.0/6.0
        elseif type == "CHAN" && length(dims) >= 4
            bf = dims[1]; h = dims[2]; tw = dims[3]; tf = dims[4]
            hw = h - 2*tf
            A = 2*bf*tf + hw*tw
            yc = (2*bf*tf*(bf/2) + hw*tw*0) / A
            I1 = (tw*hw^3)/12 + 2*(bf*tf^3/12 + bf*tf*((h-tf)/2)^2)
            I2 = (hw*tw^3)/12 + 2*(tf*bf^3/12 + bf*tf*(bf/2 - yc)^2) + hw*tw*yc^2
            J = (2*bf*tf^3 + hw*tw^3) / 3
            K1 = hw*tw/A; K2 = 5.0/6.0
        elseif type == "HAT" && length(dims) >= 4
            w_top = dims[1]; t = dims[2]; w_bot = dims[3]; h_hat = dims[4]
            h_web = h_hat - t
            A = w_top*t + 2*h_web*t + w_bot*t
            yc = (w_bot*t*t/2 + 2*h_web*t*(t + h_web/2) + w_top*t*(h_hat - t/2)) / A
            I1 = w_bot*t^3/12 + w_bot*t*(yc - t/2)^2 + 2*(t*h_web^3/12 + t*h_web*(yc - t - h_web/2)^2) + w_top*t^3/12 + w_top*t*(yc - h_hat + t/2)^2
            I2 = t*w_bot^3/12 + 2*h_web*t^3/12 + t*w_top^3/12
            J = (w_bot*t^3 + 2*h_web*t^3 + w_top*t^3) / 3
            K1 = 5.0/6.0; K2 = 5.0/6.0
        elseif (type == "T" || type == "T1" || type == "T2") && length(dims) >= 4
            h_t = dims[1]; bf = dims[2]; tw = dims[3]; tf = dims[4]
            hw = h_t - tf
            A = bf*tf + hw*tw
            yc = (hw*tw*hw/2 + bf*tf*(hw + tf/2)) / A
            I1 = tw*hw^3/12 + tw*hw*(yc - hw/2)^2 + bf*tf^3/12 + bf*tf*(hw + tf/2 - yc)^2
            I2 = hw*tw^3/12 + tf*bf^3/12
            J = (hw*tw^3 + bf*tf^3) / 3
            K1 = hw*tw/A; K2 = 5.0/6.0
        elseif type == "Z" && length(dims) >= 4
            bf = dims[1]; tf = dims[2]; hw = dims[3]; h_z = dims[4]
            tw = tf
            A = 2*bf*tf + hw*tw
            I1 = tw*hw^3/12 + 2*(bf*tf^3/12 + bf*tf*(hw/2 + tf/2)^2)
            I2 = 2*(tf*bf^3/12) + hw*tw^3/12
            J = (2*bf*tf^3 + hw*tw^3) / 3
            K1 = hw*tw/A; K2 = 5.0/6.0
        else
            if !isempty(dims)
                R = dims[1]; A = pi*R^2; I1 = pi*R^4/4; I2 = I1; J = pi*R^4/2
            end
            println("[WARNING] PBARL type '$type' not fully supported. Using approximate properties.")
        end
        # Stress recovery coefficients (y,z coordinates of 4 corner points)
        C1, C2, D1, D2, E1, E2, F1, F2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        if type == "ROD" && length(dims) >= 1
            R = dims[1]
            C1=-R; C2=0.0; D1=0.0; D2=-R; E1=R; E2=0.0; F1=0.0; F2=R
        elseif (type == "TUBE" || type == "TUBE2") && length(dims) >= 1
            R_o = dims[1]
            C1=-R_o; C2=0.0; D1=0.0; D2=-R_o; E1=R_o; E2=0.0; F1=0.0; F2=R_o
        elseif type == "BAR" && length(dims) >= 2
            w2 = dims[1]/2; h2 = dims[2]/2
            C1=h2; C2=w2; D1=-h2; D2=w2; E1=-h2; E2=-w2; F1=h2; F2=-w2
        elseif type == "BOX" && length(dims) >= 2
            w2 = dims[1]/2; h2 = dims[2]/2
            C1=h2; C2=w2; D1=-h2; D2=w2; E1=-h2; E2=-w2; F1=h2; F2=-w2
        end
        d[string(pid)] = Dict("PID"=>pid, "MID"=>mid, "A"=>A, "I"=>I1, "I1"=>I1, "I2"=>I2, "J"=>J, "TYPE"=>type, "K1"=>K1, "K2"=>K2,
            "C1"=>C1, "C2"=>C2, "D1"=>D1, "D2"=>D2, "E1"=>E1, "E2"=>E2, "F1"=>F1, "F2"=>F2)
    end
    return d
end

function extract_pbar(cards)
    d = Dict()
    for c in cards
        pid = to_id(parse_nastran_number(safe_get(c, 3)))
        mid = to_id(parse_nastran_number(safe_get(c, 4)))
        A   = parse_nastran_number(safe_get(c, 5), 0.0)
        I1  = parse_nastran_number(safe_get(c, 6), 0.0)
        I2  = parse_nastran_number(safe_get(c, 7), 0.0)
        J   = parse_nastran_number(safe_get(c, 8), 0.0)
        # Stress recovery points on continuation line 1 (indices 11-18)
        C1 = parse_nastran_number(safe_get(c, 11), 0.0)
        C2 = parse_nastran_number(safe_get(c, 12), 0.0)
        D1 = parse_nastran_number(safe_get(c, 13), 0.0)
        D2 = parse_nastran_number(safe_get(c, 14), 0.0)
        E1 = parse_nastran_number(safe_get(c, 15), 0.0)
        E2 = parse_nastran_number(safe_get(c, 16), 0.0)
        F1 = parse_nastran_number(safe_get(c, 17), 0.0)
        F2 = parse_nastran_number(safe_get(c, 18), 0.0)
        # Second continuation: K1, K2, I12
        K1 = parse_nastran_number(safe_get(c, 19), 0.0)
        K2 = parse_nastran_number(safe_get(c, 20), 0.0)
        I12 = parse_nastran_number(safe_get(c, 21), 0.0)
        if pid > 0
            d[string(pid)] = Dict(
                "PID"=>pid, "MID"=>mid, "A"=>Float64(A),
                "I1"=>Float64(I1), "I2"=>Float64(I2), "I12"=>Float64(I12), "J"=>Float64(J),
                "I"=>Float64(I1),
                "TYPE"=>"PBAR", "K1"=>Float64(K1), "K2"=>Float64(K2),
                "C1"=>Float64(C1), "C2"=>Float64(C2), "D1"=>Float64(D1), "D2"=>Float64(D2),
                "E1"=>Float64(E1), "E2"=>Float64(E2), "F1"=>Float64(F1), "F2"=>Float64(F2)
            )
        end
    end
    return d
end

function extract_pcomp(cards)
    d = Dict()
    for c in cards
        pid  = to_id(parse_nastran_number(safe_get(c, 3)))
        z0   = parse_nastran_number(safe_get(c, 4), nothing)
        nsm  = parse_nastran_number(safe_get(c, 5), 0.0)
        # Field 10 (c[10]) is LAM (SYM, MEM, BEND, SMEAR, SMCORE)
        lam_field = strip(string(safe_get(c, 10, "")))
        is_sym = uppercase(lam_field) == "SYM"

        # Parse plies starting from field 11 (continuation)
        plies = []
        k = 11
        while k + 3 <= length(c)
            ply_mid   = to_id(parse_nastran_number(safe_get(c, k), 0))
            ply_t     = parse_nastran_number(safe_get(c, k+1), 0.0)
            ply_theta = parse_nastran_number(safe_get(c, k+2), 0.0)
            ply_sout  = strip(string(safe_get(c, k+3, "")))
            if ply_mid > 0 && ply_t > 0
                push!(plies, Dict("MID"=>ply_mid, "T"=>Float64(ply_t),
                                  "THETA"=>Float64(ply_theta), "SOUT"=>ply_sout))
            end
            k += 4
        end

        # Handle SYM: mirror plies
        if is_sym && !isempty(plies)
            full_plies = copy(plies)
            for i in length(plies):-1:1
                push!(full_plies, plies[i])
            end
            plies = full_plies
        end

        # Compute total thickness
        total_t = sum(p["T"] for p in plies; init=0.0)

        if pid > 0
            d[string(pid)] = Dict("PID"=>pid, "Z0"=>isnothing(z0) ? -total_t/2 : Float64(z0),
                "NSM"=>Float64(nsm), "LAM"=>lam_field, "PLIES"=>plies,
                "T"=>Float64(total_t), "TYPE"=>"PCOMP")
        end
    end
    return d
end

function extract_prod(cards)
    d = Dict()
    for c in cards
        pid = to_id(parse_nastran_number(safe_get(c, 3)))
        mid = to_id(parse_nastran_number(safe_get(c, 4)))
        A   = parse_nastran_number(safe_get(c, 5), 0.0)
        J   = parse_nastran_number(safe_get(c, 6), 0.0)
        C   = parse_nastran_number(safe_get(c, 7), 0.0)
        NSM = parse_nastran_number(safe_get(c, 8), 0.0)
        if pid > 0
            d[string(pid)] = Dict("PID"=>pid, "MID"=>mid, "A"=>Float64(A), "J"=>Float64(J),
                                  "C"=>Float64(C), "NSM"=>Float64(NSM), "TYPE"=>"PROD")
        end
    end
    return d
end

function extract_pelas(cards)
    d = Dict()
    for c in cards
        pid = to_id(parse_nastran_number(safe_get(c, 3)))
        K   = Float64(parse_nastran_number(safe_get(c, 4), 0.0))
        GE  = Float64(parse_nastran_number(safe_get(c, 5), 0.0))
        if pid > 0
            d[string(pid)] = Dict("PID"=>pid, "K"=>K, "GE"=>GE)
        end
    end
    return d
end
