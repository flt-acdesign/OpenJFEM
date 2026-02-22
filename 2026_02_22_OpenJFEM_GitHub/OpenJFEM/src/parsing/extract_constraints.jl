# extract_constraints.jl — SPC1, SPC, SPCADD, RBE2, RBE3, MPC, MPCADD

function extract_spc1(cards)
    spcs = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        comp = string(parse_nastran_number(safe_get(c, 4), ""))
        raw_nodes = c[5:end]
        nodes = expand_nastran_list(raw_nodes)
        push!(spcs, Dict("SID"=>sid, "C"=>comp, "NODES"=>nodes))
    end
    return spcs
end

function extract_spcadd(cards)
    d = Dict()
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        raw_sets = c[4:end]
        sets = expand_nastran_list(raw_sets)
        d[sid] = sets
    end
    return d
end

function extract_spc(cards)
    spcs = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        gid = to_id(parse_nastran_number(safe_get(c, 4)))
        comp = string(parse_nastran_number(safe_get(c, 5), ""))
        if gid > 0
            push!(spcs, Dict("SID"=>sid, "C"=>comp, "NODES"=>[gid]))
        end
        gid2 = to_id(parse_nastran_number(safe_get(c, 7), 0))
        if gid2 > 0
            comp2 = string(parse_nastran_number(safe_get(c, 8), ""))
            push!(spcs, Dict("SID"=>sid, "C"=>comp2, "NODES"=>[gid2]))
        end
    end
    return spcs
end

function extract_rbe2(cards)
    d = Dict()
    for c in cards
        eid = to_id(parse_nastran_number(safe_get(c, 3)))
        gn  = to_id(parse_nastran_number(safe_get(c, 4)))  # master node
        cm  = to_id(parse_nastran_number(safe_get(c, 5)))   # constrained components (e.g. 123456)
        # Remaining fields are slave nodes
        slave_grids = Int[]
        for k in 6:length(c)
            g = to_id(parse_nastran_number(safe_get(c, k), 0))
            if g > 0; push!(slave_grids, g); end
        end
        if eid > 0 && gn > 0
            d[string(eid)] = Dict("ID"=>eid, "GN"=>gn, "CM"=>cm, "GM"=>slave_grids)
        end
    end
    return d
end

function extract_rbe3(cards)
    d = Dict()
    for c in cards
        eid = to_id(parse_nastran_number(safe_get(c, 3)))
        refgrid = to_id(parse_nastran_number(safe_get(c, 5)))
        refc = to_id(parse_nastran_number(safe_get(c, 6)))
        wt = parse_nastran_number(safe_get(c, 7), 1.0)
        comps = to_id(parse_nastran_number(safe_get(c, 8)))
        raw_grids = c[9:end]
        dep_grids = expand_nastran_list(raw_grids)
        d[string(eid)] = Dict("ID"=>eid, "REFGRID"=>refgrid, "REFC"=>refc, "WT"=>wt, "COMPS"=>comps, "DEP_GRIDS"=>dep_grids)
    end
    return d
end

function extract_mpc(cards)
    mpcs = []
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        terms = []
        k = 4
        while k + 2 <= length(c)
            gid   = to_id(parse_nastran_number(safe_get(c, k), 0))
            comp  = to_id(parse_nastran_number(safe_get(c, k+1), 0))
            coeff = Float64(parse_nastran_number(safe_get(c, k+2), 0.0))
            if gid > 0 && comp > 0
                push!(terms, Dict("G"=>gid, "C"=>comp, "A"=>coeff))
            end
            k += 3
        end
        if length(terms) >= 2
            push!(mpcs, Dict("SID"=>sid, "TERMS"=>terms))
        end
    end
    return mpcs
end

function extract_mpcadd(cards)
    d = Dict()
    for c in cards
        sid = to_id(parse_nastran_number(safe_get(c, 3)))
        subs = Int[]
        for k in 4:length(c)
            s = to_id(parse_nastran_number(safe_get(c, k), 0))
            if s > 0; push!(subs, s); end
        end
        if sid > 0 && !isempty(subs)
            d[sid] = subs
        end
    end
    return d
end
