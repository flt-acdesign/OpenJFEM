# extract_elements.jl — CROD, CELAS1, CONROD, CONM2

function extract_crod(cards)
    d = Dict()
    for c in cards
        id  = to_id(parse_nastran_number(safe_get(c, 3)))
        pid = to_id(parse_nastran_number(safe_get(c, 4)))
        ga  = to_id(parse_nastran_number(safe_get(c, 5)))
        gb  = to_id(parse_nastran_number(safe_get(c, 6)))
        if id > 0
            d[string(id)] = Dict("ID"=>id, "PID"=>pid, "GA"=>ga, "GB"=>gb, "TYPE"=>"CROD")
        end
    end
    return d
end

function extract_celas1(cards)
    d = Dict()
    for c in cards
        eid = to_id(parse_nastran_number(safe_get(c, 3)))
        pid = to_id(parse_nastran_number(safe_get(c, 4)))
        g1  = to_id(parse_nastran_number(safe_get(c, 5), 0))
        c1  = to_id(parse_nastran_number(safe_get(c, 6), 0))
        g2  = to_id(parse_nastran_number(safe_get(c, 7), 0))
        c2  = to_id(parse_nastran_number(safe_get(c, 8), 0))
        if eid > 0
            d[string(eid)] = Dict("ID"=>eid, "PID"=>pid,
                                  "G1"=>g1, "C1"=>c1, "G2"=>g2, "C2"=>c2, "TYPE"=>"CELAS1")
        end
    end
    return d
end

function extract_conrod(cards)
    d = Dict()
    for c in cards
        eid = to_id(parse_nastran_number(safe_get(c, 3)))
        g1  = to_id(parse_nastran_number(safe_get(c, 4)))
        g2  = to_id(parse_nastran_number(safe_get(c, 5)))
        mid = to_id(parse_nastran_number(safe_get(c, 6)))
        A   = Float64(parse_nastran_number(safe_get(c, 7), 0.0))
        J   = Float64(parse_nastran_number(safe_get(c, 8), 0.0))
        C_  = Float64(parse_nastran_number(safe_get(c, 9), 0.0))
        NSM = Float64(parse_nastran_number(safe_get(c, 10), 0.0))
        if eid > 0
            d[string(eid)] = Dict("ID"=>eid, "GA"=>g1, "GB"=>g2,
                "MID"=>mid, "A"=>A, "J"=>J, "C"=>C_, "NSM"=>NSM, "TYPE"=>"CONROD")
        end
    end
    return d
end

function extract_conm2(cards)
    d = Dict()
    for c in cards
        eid  = to_id(parse_nastran_number(safe_get(c, 3)))
        gid  = to_id(parse_nastran_number(safe_get(c, 4)))
        cid  = to_id(parse_nastran_number(safe_get(c, 5), 0))
        mass = Float64(parse_nastran_number(safe_get(c, 6), 0.0))
        x1   = Float64(parse_nastran_number(safe_get(c, 7), 0.0))
        x2   = Float64(parse_nastran_number(safe_get(c, 8), 0.0))
        x3   = Float64(parse_nastran_number(safe_get(c, 9), 0.0))
        I11  = Float64(parse_nastran_number(safe_get(c, 11), 0.0))
        I21  = Float64(parse_nastran_number(safe_get(c, 12), 0.0))
        I22  = Float64(parse_nastran_number(safe_get(c, 13), 0.0))
        I31  = Float64(parse_nastran_number(safe_get(c, 14), 0.0))
        I32  = Float64(parse_nastran_number(safe_get(c, 15), 0.0))
        I33  = Float64(parse_nastran_number(safe_get(c, 16), 0.0))
        if eid > 0
            d[string(eid)] = Dict("ID"=>eid, "GID"=>gid, "CID"=>cid,
                "M"=>mass, "X"=>[x1, x2, x3],
                "I"=>[I11, I21, I22, I31, I32, I33], "TYPE"=>"CONM2")
        end
    end
    return d
end
