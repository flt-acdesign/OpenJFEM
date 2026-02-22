# card_processor.jl — Card processing and case control / bulk data splitting

function process_cards(lines)
    processed = Dict{String, Vector{Any}}()
    i = 1
    while i <= length(lines)
        line = lines[i]
        clean_line = strip(line)
        if startswith(clean_line, '$') || isempty(clean_line)
            i += 1; continue
        end

        name = get_nastran_card_name(line)
        if startswith(name, "+") || startswith(name, "*")
             i += 1; continue
        end

        if !isnothing(name) && !isempty(name)
            fields = Any["SMALL"; String(name)]
            append!(fields, get_nastran_fields_from_line(line))

            steps = 1
            while i + steps <= length(lines)
                next_line = lines[i+steps]
                next_clean = strip(next_line)
                if startswith(next_clean, '$')
                    steps += 1; continue
                end

                is_cont = false
                if startswith(next_line, " ") || startswith(next_clean, "+") || startswith(next_clean, "*")
                    is_cont = true
                elseif occursin(",", next_line) && startswith(next_clean, ",")
                    is_cont = true
                end

                if is_cont
                    append!(fields, get_nastran_fields_from_line(next_line))
                    steps += 1
                else
                    break
                end
            end

            if !haskey(processed, name); processed[name] = []; end
            push!(processed[name], fields)
            i += steps
        else
            i += 1
        end
    end
    return processed
end

function read_bulk_and_case(lines::Vector{String})
    case_control = Dict("SUBCASES" => Dict{Int, Dict{String, Any}}())
    bulk_lines = String[]
    in_bulk = false
    global_load, global_spc, global_mpc = nothing, nothing, nothing
    current_sub = 0

    for line in lines
        cl = uppercase(split(line, '$')[1])
        if occursin("BEGIN BULK", cl); in_bulk=true; continue; end
        if occursin("ENDDATA", cl); break; end

        if !in_bulk
            stripped_cl = strip(cl)
            # Skip MYSTRAN-specific case control keywords
            if startswith(stripped_cl, "ELDATA"); continue; end
            if startswith(stripped_cl, "LABEL"); continue; end
            if startswith(stripped_cl, "ECHO"); continue; end
            if startswith(stripped_cl, "SET "); continue; end
            if startswith(stripped_cl, "GPFORCE"); continue; end
            if startswith(stripped_cl, "MPCFORCE"); continue; end
            if startswith(stripped_cl, "OLOAD"); continue; end
            if startswith(stripped_cl, "TITLE"); continue; end
            if startswith(stripped_cl, "SUBTI"); continue; end

            if startswith(stripped_cl, "SUBCASE")
                parts = split(cl)
                if length(parts) >= 2
                    val = try parse(Int, parts[2]) catch; 0 end
                    if val > 0
                        current_sub = val
                        case_control["SUBCASES"][current_sub] = Dict{String, Any}("LOAD"=>global_load, "SPC"=>global_spc, "MPC"=>global_mpc)
                    end
                end
            elseif occursin("=", cl)
                # Parse key = value, handling parenthetical modifiers like DISP(PRINT,PLOT) = ALL
                eq_parts = split(cl, "="; limit=2)
                k_raw = strip(eq_parts[1])
                v_raw = strip(eq_parts[2])
                # Strip parenthetical modifiers from key: "DISPLACEMENT(SORT1,REAL)" -> "DISPLACEMENT"
                k = replace(k_raw, r"\(.*\)" => "")
                k = strip(k)
                val = try parse(Int, v_raw) catch; v_raw end
                if current_sub > 0
                    case_control["SUBCASES"][current_sub][k] = val
                else
                    if k == "LOAD"; global_load = val; end
                    if k == "SPC"; global_spc = val; end
                    if k == "MPC"; global_mpc = val; end
                end
            end
        else
            if length(strip(cl)) > 1
                push!(bulk_lines, String(rstrip(cl)))
            end
        end
    end

    # Create default subcase if none defined but global load/spc exists
    if isempty(case_control["SUBCASES"]) && (!isnothing(global_load) || !isnothing(global_spc))
        case_control["SUBCASES"][1] = Dict{String, Any}("LOAD"=>global_load, "SPC"=>global_spc, "MPC"=>global_mpc)
        println("[INFO] No SUBCASE defined. Created default SUBCASE 1 with LOAD=$global_load, SPC=$global_spc")
    end

    return case_control, bulk_lines
end
