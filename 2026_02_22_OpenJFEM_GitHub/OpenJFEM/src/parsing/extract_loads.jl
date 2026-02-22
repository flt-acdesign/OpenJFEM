# extract_loads.jl — FORCE, MOMENT, PLOAD4, PLOAD2, PLOAD1, GRAV, LOAD combos

function extract_loads(cards)
    f = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        gid = to_id(parse_nastran_number(safe_get(c, 4)))
        cid = to_id(parse_nastran_number(safe_get(c, 5), 0))
        mag = parse_nastran_number(safe_get(c, 6), 0.0)
        dir = [parse_nastran_number(safe_get(c, 7),0.0), parse_nastran_number(safe_get(c, 8),0.0), parse_nastran_number(safe_get(c, 9),0.0)]
        push!(f, Dict("TYPE"=>"FORCE", "SID"=>sid, "GID"=>gid, "CID"=>cid, "Mag"=>mag, "Dir"=>dir))
    end
    return f
end

function extract_moments(cards)
    m = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        gid = to_id(parse_nastran_number(safe_get(c, 4)))
        cid = to_id(parse_nastran_number(safe_get(c, 5), 0))
        mag = parse_nastran_number(safe_get(c, 6), 0.0)
        dir = [parse_nastran_number(safe_get(c, 7),0.0), parse_nastran_number(safe_get(c, 8),0.0), parse_nastran_number(safe_get(c, 9),0.0)]
        push!(m, Dict("TYPE"=>"MOMENT", "SID"=>sid, "GID"=>gid, "CID"=>cid, "Mag"=>mag, "Dir"=>dir))
    end
    return m
end

function extract_pload4(cards)
    p = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        eid = to_id(parse_nastran_number(safe_get(c, 4)))
        press = [parse_nastran_number(safe_get(c, 5), 0.0)]
        push!(p, Dict("TYPE"=>"PLOAD4", "SID"=>sid, "EID"=>eid, "P"=>press[1]))
    end
    return p
end

function extract_pload2(cards)
    p = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        press = parse_nastran_number(safe_get(c, 4), 0.0)
        for k in 5:length(c)
            eid = to_id(parse_nastran_number(safe_get(c, k), 0))
            if eid > 0
                push!(p, Dict("TYPE"=>"PLOAD4", "SID"=>sid, "EID"=>eid, "P"=>press))
            end
        end
    end
    return p
end

function extract_pload1(cards)
    p = []
    for c in cards
        sid   = to_id(parse_nastran_number(safe_get(c, 3)))
        eid   = to_id(parse_nastran_number(safe_get(c, 4)))
        ltype = to_id(parse_nastran_number(safe_get(c, 5), 0))
        scale_str = strip(string(safe_get(c, 6, "")))
        x1    = Float64(parse_nastran_number(safe_get(c, 7), 0.0))
        p1    = Float64(parse_nastran_number(safe_get(c, 8), 0.0))
        x2    = Float64(parse_nastran_number(safe_get(c, 9), 1.0))
        p2    = Float64(parse_nastran_number(safe_get(c, 10), 0.0))
        if sid > 0 && eid > 0
            push!(p, Dict("TYPE"=>"PLOAD1", "SID"=>sid, "EID"=>eid,
                          "LOAD_TYPE"=>ltype, "SCALE"=>scale_str,
                          "X1"=>x1, "P1"=>p1, "X2"=>x2, "P2"=>p2))
        end
    end
    return p
end

function extract_grav(cards)
    g = []
    for c in cards
        sid   = to_id(parse_nastran_number(safe_get(c, 3)))
        cid   = to_id(parse_nastran_number(safe_get(c, 4), 0))
        accel = Float64(parse_nastran_number(safe_get(c, 5), 0.0))
        n1    = Float64(parse_nastran_number(safe_get(c, 6), 0.0))
        n2    = Float64(parse_nastran_number(safe_get(c, 7), 0.0))
        n3    = Float64(parse_nastran_number(safe_get(c, 8), 0.0))
        if sid > 0
            push!(g, Dict("TYPE"=>"GRAV", "SID"=>sid, "CID"=>cid,
                           "A"=>accel, "N"=>[n1, n2, n3]))
        end
    end
    return g
end

function extract_load_combos(cards)
    combos = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        s = parse_nastran_number(safe_get(c, 4), 1.0)
        comps = []
        for i in 5:2:length(c)-1
            s_i = parse_nastran_number(safe_get(c, i), nothing)
            l_i = to_id(parse_nastran_number(safe_get(c, i+1), nothing))
            if !isnothing(s_i) && l_i > 0
                push!(comps, Dict("S"=>s_i, "LID"=>l_i))
            end
        end
        push!(combos, Dict("SID"=>sid, "S"=>s, "COMPS"=>comps))
    end
    return combos
end
