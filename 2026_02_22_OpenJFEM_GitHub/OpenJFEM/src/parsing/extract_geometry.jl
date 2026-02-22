# extract_geometry.jl — GRID, coordinate systems, shell elements, CBAR

function extract_grid(cards)
    d = Dict()
    for c in cards
        id = to_id(parse_nastran_number(safe_get(c, 3)))
        cp = to_id(parse_nastran_number(safe_get(c, 4), 0))
        x = [parse_nastran_number(safe_get(c, 5), 0.0),
             parse_nastran_number(safe_get(c, 6), 0.0),
             parse_nastran_number(safe_get(c, 7), 0.0)]
        cd = to_id(parse_nastran_number(safe_get(c, 8), 0))
        if id > 0; d[string(id)] = Dict("ID"=>id, "CP"=>cp, "CD"=>cd, "X"=>x); end
    end
    return d
end

function extract_coords(cards; coord_type::String="RECTANGULAR")
    d = Dict()
    for c in cards
        id = to_id(parse_nastran_number(safe_get(c, 3)))
        A = [parse_nastran_number(safe_get(c, 5),0.0), parse_nastran_number(safe_get(c, 6),0.0), parse_nastran_number(safe_get(c, 7),0.0)]
        B = [parse_nastran_number(safe_get(c, 8),0.0), parse_nastran_number(safe_get(c, 9),0.0), parse_nastran_number(safe_get(c, 10),0.0)]
        C = [parse_nastran_number(safe_get(c, 11),0.0), parse_nastran_number(safe_get(c, 12),0.0), parse_nastran_number(safe_get(c, 13),0.0)]
        w = B - A
        if norm(w) < 1e-9; w=[0.0,0.0,1.0]; else; w=normalize(w); end
        v_t = C - A
        v = cross(w, v_t)
        if norm(v) < 1e-9; v=[0.0,1.0,0.0]; else; v=normalize(v); end
        u = normalize(cross(v, w))
        d[string(id)] = Dict("Origin"=>A, "U"=>u, "V"=>v, "W"=>w, "TYPE"=>coord_type)
    end
    return d
end

function extract_shells(cards)
    d = Dict()
    for c in cards
        id = to_id(parse_nastran_number(safe_get(c, 3)))
        pid = to_id(parse_nastran_number(safe_get(c, 4)))
        # Determine node count from card name: CTRIA3 has 3, CQUAD4 has 4
        card_name = uppercase(strip(string(safe_get(c, 2, ""))))
        n_nodes = startswith(card_name, "CTRIA") ? 3 : 4
        nodes = []
        for k in 5:(4+n_nodes)
            val = to_id(parse_nastran_number(safe_get(c, k), 0))
            if val > 0; push!(nodes, val); end
        end
        # THETA/MCID field: immediately after node fields
        theta = Float64(parse_nastran_number(safe_get(c, 5+n_nodes), 0.0))
        if id > 0; d[string(id)] = Dict("ID"=>id, "PID"=>pid, "NODES"=>nodes, "THETA"=>theta); end
    end
    return d
end

function extract_cbar(cards)
    d = Dict()
    for c in cards
        id = to_id(parse_nastran_number(safe_get(c, 3)))
        pid = to_id(parse_nastran_number(safe_get(c, 4)))
        ga = to_id(parse_nastran_number(safe_get(c, 5)))
        gb = to_id(parse_nastran_number(safe_get(c, 6)))

        # Detect G0 (grid point) vs X1,X2,X3 (vector) format
        x1_raw = strip(string(safe_get(c, 7, "")))
        x2_raw = strip(string(safe_get(c, 8, "")))
        x3_raw = strip(string(safe_get(c, 9, "")))

        g0 = 0
        v = [0.0, 0.0, 0.0]
        if !isempty(x1_raw) && !occursin(".", x1_raw) && (isempty(x2_raw) || x2_raw == "0")
            # G0 format: integer grid ID, X2/X3 blank
            g0 = to_id(parse_nastran_number(x1_raw, 0))
        else
            v = [parse_nastran_number(safe_get(c, 7),0.0), parse_nastran_number(safe_get(c, 8),0.0), parse_nastran_number(safe_get(c, 9),0.0)]
        end
        if g0 == 0 && norm(v) < 1e-6; v = [0.0, 0.0, 1.0]; end

        # Parse WA/WB offset vectors from continuation line (fields 13-18)
        wa = [parse_nastran_number(safe_get(c, 13), 0.0), parse_nastran_number(safe_get(c, 14), 0.0), parse_nastran_number(safe_get(c, 15), 0.0)]
        wb = [parse_nastran_number(safe_get(c, 16), 0.0), parse_nastran_number(safe_get(c, 17), 0.0), parse_nastran_number(safe_get(c, 18), 0.0)]

        d[string(id)] = Dict("ID"=>id, "PID"=>pid, "GA"=>ga, "GB"=>gb, "V"=>v, "G0"=>g0, "WA"=>wa, "WB"=>wb, "TYPE"=>"CBAR")
    end
    return d
end
