# mystran_converter.jl — Convert MYSTRAN-format BDF files to NASTRAN format

function convert_mystran_to_nastran(lines::Vector{String})
    is_mystran = false
    for line in lines
        cl = strip(uppercase(line))
        if startswith(cl, "SOL ")
            m = match(r"SOL\s+(\d+)", cl)
            if m !== nothing
                sol_num = parse(Int, m.captures[1])
                if sol_num < 100; is_mystran = true; end
            end
            break
        end
    end

    if is_mystran
        println("[INFO] MYSTRAN format detected. Converting to NASTRAN format...")
    end

    # Phase 1: Collect GRID PS constraints (always, for both NASTRAN and MYSTRAN)
    in_bulk = false
    grid_ps_entries = Dict{Int, String}()
    existing_spc_sid = nothing

    for line in lines
        cl = strip(uppercase(split(line, '$')[1]))
        if occursin("BEGIN BULK", cl); in_bulk = true; continue; end
        if occursin("ENDDATA", cl); break; end
        if !in_bulk
            if occursin("SPC", cl) && occursin("=", cl)
                parts = split(cl, "=")
                k = strip(replace(parts[1], r"\(.*\)" => ""))
                if k == "SPC"
                    existing_spc_sid = try parse(Int, strip(parts[2])) catch; nothing end
                end
            end
        else
            name = ""
            if occursin(",", cl)
                name = uppercase(strip(split(cl, ",")[1]))
            elseif length(cl) >= 4
                name = uppercase(strip(cl[1:min(8, length(cl))]))
            end
            if name == "GRID"
                if occursin(",", cl)
                    parts = split(cl, ",")
                    # Free-field GRID: GRID,ID,CP,X1,X2,X3,CD,PS,SEID
                    # PS is field index 8 (1-based after GRID)
                    if length(parts) >= 8
                        ps_field = strip(parts[8])
                        if !isempty(ps_field) && ps_field != "0"
                            gid = try parse(Int, strip(parts[2])) catch; 0 end
                            if gid > 0; grid_ps_entries[gid] = ps_field; end
                        end
                    end
                else
                    padded = rpad(cl, 80, ' ')
                    # Fixed-field GRID: cols 1-8=GRID, 9-16=ID, 17-24=CP, 25-32=X1, 33-40=X2, 41-48=X3, 49-56=CD, 57-64=PS
                    if length(padded) >= 64
                        ps_str = strip(padded[57:64])
                        if !isempty(ps_str) && ps_str != "0"
                            id_str = strip(padded[9:16])
                            gid = try parse(Int, id_str) catch; 0 end
                            if gid > 0; grid_ps_entries[gid] = ps_str; end
                        end
                    end
                end
            end
        end
    end

    # If not MYSTRAN and no PS fields, nothing to convert
    if !is_mystran && isempty(grid_ps_entries); return lines; end

    ps_spc_sid = 999999

    # Phase 2: Transform lines
    result = String[]
    in_bulk = false

    for line in lines
        cl = strip(uppercase(line))

        if is_mystran && startswith(cl, "ID ") && !in_bulk; continue; end

        if is_mystran && startswith(cl, "SOL ") && !in_bulk
            m = match(r"SOL\s+(\d+)", cl)
            if m !== nothing
                sol_num = parse(Int, m.captures[1])
                if sol_num < 100
                    push!(result, "SOL $(sol_num + 100)")
                    println("[INFO]   SOL $sol_num -> SOL $(sol_num + 100)")
                    continue
                end
            end
        end

        if occursin("BEGIN BULK", cl)
            in_bulk = true
            if !isempty(grid_ps_entries)
                if isnothing(existing_spc_sid)
                    push!(result, "SPC = $ps_spc_sid")
                else
                    spcadd_sid = ps_spc_sid - 1
                    push!(result, "SPC = $spcadd_sid")
                end
            end
            push!(result, line)
            if !isempty(grid_ps_entries)
                for (gid, ps) in grid_ps_entries
                    spc1_line = rpad("SPC1", 8) * rpad(string(ps_spc_sid), 8) * rpad(ps, 8) * rpad(string(gid), 8)
                    push!(result, spc1_line)
                end
                println("[INFO]   Generated $(length(grid_ps_entries)) SPC1 cards from GRID PS fields")
                if !isnothing(existing_spc_sid)
                    spcadd_sid = ps_spc_sid - 1
                    spcadd_line = rpad("SPCADD", 8) * rpad(string(spcadd_sid), 8) * rpad(string(existing_spc_sid), 8) * rpad(string(ps_spc_sid), 8)
                    push!(result, spcadd_line)
                end
            end
            continue
        end

        if in_bulk
            card_name = ""
            if occursin(",", cl)
                card_name = uppercase(strip(split(cl, ",")[1]))
            elseif length(cl) >= 4
                card_name = uppercase(strip(cl[1:min(8, length(cl))]))
            end

            if is_mystran && card_name == "DEBUG"; continue; end
            if is_mystran && card_name == "PARAM" && occursin("SOLLIB", cl); continue; end

            if card_name == "GRID"
                if occursin(",", line)
                    parts = split(line, ",")
                    if length(parts) >= 8; parts[8] = ""; end
                    push!(result, join(parts, ","))
                    continue
                else
                    padded = rpad(line, 80, ' ')
                    clean = padded[1:min(56, length(padded))] * "                "
                    push!(result, rstrip(clean))
                    continue
                end
            end
        else
            if is_mystran && startswith(cl, "ELDATA"); continue; end
        end

        push!(result, line)
    end

    return result
end
