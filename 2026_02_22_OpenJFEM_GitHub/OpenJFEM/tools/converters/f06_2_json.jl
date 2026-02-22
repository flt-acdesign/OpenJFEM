
using JSON

"""
Parses Nastran numbers. Handles "1.0-5" shorthand and standard "1.0E-5".
"""
function parse_nastran_float(s::AbstractString)
    s = strip(s)
    if isempty(s); return 0.0; end
    
    
    val = tryparse(Float64, s)
    if val !== nothing; return val; end

    # Handle "1.234-5" -> "1.234E-5" (sign without E)
    try
        # Loop backwards to find a +/- that isn't E+ or E-
        for i in length(s):-1:2
            c = s[i]
            prev = s[i-1]
            if (c == '+' || c == '-') && (prev != 'E' && prev != 'e')
                new_s = s[1:i-1] * "E" * s[i:end]
                return parse(Float64, new_s)
            end
        end
    catch; end
    return 0.0
end

"""
Checks header by removing spaces. "D I S P" -> "DISP"
"""
function is_header(clean_line::AbstractString, keyword::String)
    # Remove all spaces to match "D I S P L A C E M E N T"
    normalized = replace(clean_line, " " => "")
    return occursin(keyword, normalized)
end

function parse_f06_to_json(file_content::String)
    
    results = Dict(
        "displacements" => [],
        "spc_forces" => [],
        "forces" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => []),
        "stresses" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => []),
        "strains" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => [])
    )

    current_section = :none
    cbar_buffer = Dict() 

    
    lines = split(file_content, "\n")
    println("... Scanning $(length(lines)) lines...")

    for line in lines
        # 1. Handle Page Breaks: NASTRAN F06 carriage control character '1' at position 1
        #    Page breaks: "1    TITLE..." (no leading whitespace)
        #    Data lines:  "             1      G      0.0..." (leading whitespace before ID)
        if !isempty(line) && line[1] == '1' && (length(line) == 1 || !isdigit(line[min(2,end)]))
            continue
        end

        
        raw_line = string(line)
        line_stripped = strip(raw_line)
        
        if isempty(line_stripped); continue; end

        
        if is_header(line_stripped, "DISPLACEMENTVECTOR")
            current_section = :displacement; continue
        elseif is_header(line_stripped, "LOADVECTOR")
            current_section = :skip; continue
        elseif is_header(line_stripped, "FORCESOFSINGLE-POINTCONSTRAINT")
            current_section = :spc_force; continue
        elseif is_header(line_stripped, "FORCESINBARELEMENTS")
            current_section = :force_cbar; continue
        elseif is_header(line_stripped, "FORCESINQUADRILATERAL")
            current_section = :force_quad4; continue
        elseif is_header(line_stripped, "FORCESINTRIANGULAR")
            current_section = :force_tria3; continue
        elseif is_header(line_stripped, "STRESSESINBARELEMENTS")
            current_section = :stress_cbar; continue
        elseif is_header(line_stripped, "STRAINSINBARELEMENTS")
            current_section = :strain_cbar; continue
        elseif is_header(line_stripped, "STRESSESINQUADRILATERAL")
            current_section = :stress_quad4; continue
        elseif is_header(line_stripped, "STRAINSINQUADRILATERAL")
             current_section = :strain_quad4; continue
        elseif is_header(line_stripped, "STRESSESINTRIANGULAR")
            current_section = :stress_tria3; continue
        elseif is_header(line_stripped, "STRAINSINTRIANGULAR")
            current_section = :strain_tria3; continue
        elseif is_header(line_stripped, "FORCESINRODELEMENTS")
            current_section = :force_crod; continue
        elseif is_header(line_stripped, "STRESSESINRODELEMENTS")
            current_section = :stress_crod; continue
        elseif is_header(line_stripped, "STRAINSINRODELEMENTS")
            current_section = :strain_crod; continue
        elseif is_header(line_stripped, "OLOADRESULTANT")
            current_section = :skip; continue
        elseif is_header(line_stripped, "ANALYSISSUMMARYTABLE")
            current_section = :skip; continue
        elseif is_header(line_stripped, "DBDICTPRINT")
            current_section = :skip; continue
        end

        
        if is_header(line_stripped, "POINTID") || is_header(line_stripped, "ELEMENTID"); continue; end
        if is_header(line_stripped, "SUBCASE") || is_header(line_stripped, "DAREA"); continue; end
        if is_header(line_stripped, "BEND-MOMENT") || occursin("MSC Nastran", line_stripped); continue; end
        # Skip ROD element table headers
        if is_header(line_stripped, "AXIALFORCE") && is_header(line_stripped, "TORQUE"); continue; end
        if is_header(line_stripped, "AXIALSTRESS") && is_header(line_stripped, "MARGIN"); continue; end
        if is_header(line_stripped, "AXIALSTRAIN") && is_header(line_stripped, "MARGIN"); continue; end
        # Skip CBAR stress/strain table headers (SA1/SB1 etc.)
        if occursin("SA1", line_stripped) && occursin("SA2", line_stripped); continue; end
        if occursin("SB1", line_stripped) && occursin("SB2", line_stripped); continue; end
        # Skip any line with "----" separators
        if occursin("----", line_stripped) && !occursin("E", line_stripped) && !occursin("e", line_stripped); continue; end

        
        try
            parts = split(line_stripped)
            if isempty(parts); continue; end

            # --- HANDLE PRINTER CONTROL CHAR '0' (Double Space) ---
            # Case A: Detached "0" -> ["0", "123", ...]
            if parts[1] == "0"
                popfirst!(parts)
            # Case B: Attached "0" -> ["0123", ...] (Only if 0 is followed by digits)
            elseif startswith(parts[1], "0") && length(parts[1]) > 1 && isdigit(parts[1][2])
                
                parts[1] = parts[1][2:end]
            end
            
            if isempty(parts); continue; end

            
            
            if current_section == :displacement
                id = tryparse(Int, parts[1])
                if id !== nothing
                    # Check if 'G' or 'S' type is present
                    idx = (length(parts) == 8) ? 3 : 2 
                    if length(parts) >= idx+5
                         push!(results["displacements"], Dict(
                            "grid_id" => id,
                            "t1" => parse_nastran_float(parts[idx]),
                            "t2" => parse_nastran_float(parts[idx+1]),
                            "t3" => parse_nastran_float(parts[idx+2]),
                            "r1" => parse_nastran_float(parts[idx+3]),
                            "r2" => parse_nastran_float(parts[idx+4]),
                            "r3" => parse_nastran_float(parts[idx+5])
                        ))
                    end
                end

            
            elseif current_section == :spc_force
                id = tryparse(Int, parts[1])
                if id !== nothing
                    idx = (length(parts) == 8) ? 3 : 2
                    if length(parts) >= idx+5
                        push!(results["spc_forces"], Dict(
                            "grid_id" => id,
                            "t1" => parse_nastran_float(parts[idx]),
                            "t2" => parse_nastran_float(parts[idx+1]),
                            "t3" => parse_nastran_float(parts[idx+2]),
                            "r1" => parse_nastran_float(parts[idx+3]),
                            "r2" => parse_nastran_float(parts[idx+4]),
                            "r3" => parse_nastran_float(parts[idx+5])
                        ))
                    end
                end

            
            elseif current_section == :force_cbar
                id = tryparse(Int, parts[1])
                if id !== nothing && length(parts) >= 9
                    push!(results["forces"]["cbar"], Dict(
                        "eid" => id,
                        "moment_a1" => parse_nastran_float(parts[2]),
                        "moment_a2" => parse_nastran_float(parts[3]),
                        "moment_b1" => parse_nastran_float(parts[4]),
                        "moment_b2" => parse_nastran_float(parts[5]),
                        "shear_1"   => parse_nastran_float(parts[6]),
                        "shear_2"   => parse_nastran_float(parts[7]),
                        "axial"     => parse_nastran_float(parts[8]),
                        "torque"    => parse_nastran_float(parts[9])
                    ))
                end

            
            elseif (current_section == :force_quad4 || current_section == :force_tria3)
                key = (current_section == :force_quad4) ? "quad4" : "tria3"
                id = tryparse(Int, parts[1])
                if id !== nothing && length(parts) >= 9
                    push!(results["forces"][key], Dict(
                        "eid" => id,
                        "fx" => parse_nastran_float(parts[2]), "fy" => parse_nastran_float(parts[3]),
                        "fxy" => parse_nastran_float(parts[4]), "mx" => parse_nastran_float(parts[5]),
                        "my" => parse_nastran_float(parts[6]), "mxy" => parse_nastran_float(parts[7]),
                        "qx" => parse_nastran_float(parts[8]), "qy" => parse_nastran_float(parts[9])
                    ))
                end

            
            elseif current_section in [:stress_cbar, :strain_cbar]
                cat_key = (current_section == :stress_cbar) ? "stresses" : "strains"
                
                
                id = tryparse(Int, parts[1])
                
                if id !== nothing && length(parts) >= 6
                    
                    cbar_buffer = Dict(
                        "eid" => id,
                        "end_a" => Dict(
                            "p1" => parse_nastran_float(parts[2]), "p2" => parse_nastran_float(parts[3]), 
                            "p3" => parse_nastran_float(parts[4]), "p4" => parse_nastran_float(parts[5])
                        ),
                        "axial" => parse_nastran_float(parts[6])
                    )
                elseif !isempty(cbar_buffer) && length(parts) >= 4
                    
                    cbar_buffer["end_b"] = Dict(
                        "p1" => parse_nastran_float(parts[1]), "p2" => parse_nastran_float(parts[2]), 
                        "p3" => parse_nastran_float(parts[3]), "p4" => parse_nastran_float(parts[4])
                    )
                    push!(results[cat_key]["cbar"], deepcopy(cbar_buffer))
                    empty!(cbar_buffer)
                end

            
            # --- CROD Forces ---
            elseif current_section == :force_crod
                id = tryparse(Int, parts[1])
                if id !== nothing && length(parts) >= 3
                    push!(results["forces"]["crod"], Dict(
                        "eid" => id,
                        "axial" => parse_nastran_float(parts[2]),
                        "torque" => parse_nastran_float(parts[3])
                    ))
                    # F06 can have two elements per line side by side
                    if length(parts) >= 6
                        id2 = tryparse(Int, parts[4])
                        if id2 !== nothing
                            push!(results["forces"]["crod"], Dict(
                                "eid" => id2,
                                "axial" => parse_nastran_float(parts[5]),
                                "torque" => parse_nastran_float(parts[6])
                            ))
                        end
                    end
                end

            # --- CROD Stresses/Strains ---
            # F06 format: two elements per line, safety margins may be blank
            # With margins:    EID1  AXIAL1  SM_T1  TORSIONAL1  SM_C1    EID2  AXIAL2  SM_T2  TORSIONAL2  SM_C2  (10 fields)
            # Without margins: EID1  AXIAL1  TORSIONAL1    EID2  AXIAL2  TORSIONAL2  (6 fields)
            # Detect by checking if parts[4] looks like an integer (EID)
            elseif current_section in [:stress_crod, :strain_crod]
                cat_key = (current_section == :stress_crod) ? "stresses" : "strains"
                id = tryparse(Int, parts[1])
                if id !== nothing && length(parts) >= 2
                    # Determine if safety margins are present
                    # If parts[4] exists and parses as an integer, margins are absent (parts[4] is EID2)
                    has_margins = true
                    if length(parts) >= 4
                        id_test = tryparse(Int, parts[4])
                        if id_test !== nothing && !occursin(".", parts[4]) && !occursin("E", uppercase(parts[4]))
                            has_margins = false
                        end
                    end

                    if has_margins
                        # 5 fields per element: EID AXIAL SM_T TORSIONAL SM_C
                        entry = Dict(
                            "eid" => id,
                            "axial" => parse_nastran_float(parts[2]),
                            "safety_margin_t" => length(parts) >= 3 ? parse_nastran_float(parts[3]) : 0.0,
                            "torsional" => length(parts) >= 4 ? parse_nastran_float(parts[4]) : 0.0,
                            "safety_margin_c" => length(parts) >= 5 ? parse_nastran_float(parts[5]) : 0.0
                        )
                        push!(results[cat_key]["crod"], entry)
                        if length(parts) >= 10
                            id2 = tryparse(Int, parts[6])
                            if id2 !== nothing
                                push!(results[cat_key]["crod"], Dict(
                                    "eid" => id2,
                                    "axial" => parse_nastran_float(parts[7]),
                                    "safety_margin_t" => parse_nastran_float(parts[8]),
                                    "torsional" => parse_nastran_float(parts[9]),
                                    "safety_margin_c" => parse_nastran_float(parts[10])
                                ))
                            end
                        end
                    else
                        # 3 fields per element: EID AXIAL TORSIONAL (margins blank)
                        entry = Dict(
                            "eid" => id,
                            "axial" => parse_nastran_float(parts[2]),
                            "torsional" => length(parts) >= 3 ? parse_nastran_float(parts[3]) : 0.0,
                            "safety_margin_t" => 0.0,
                            "safety_margin_c" => 0.0
                        )
                        push!(results[cat_key]["crod"], entry)
                        # Two elements per line
                        if length(parts) >= 6
                            id2 = tryparse(Int, parts[4])
                            if id2 !== nothing
                                push!(results[cat_key]["crod"], Dict(
                                    "eid" => id2,
                                    "axial" => parse_nastran_float(parts[5]),
                                    "torsional" => length(parts) >= 6 ? parse_nastran_float(parts[6]) : 0.0,
                                    "safety_margin_t" => 0.0,
                                    "safety_margin_c" => 0.0
                                ))
                            end
                        end
                    end
                end

            elseif current_section in [:stress_quad4, :stress_tria3, :strain_quad4, :strain_tria3]
                is_stress = occursin("stress", string(current_section))
                elem_type = occursin("quad", string(current_section)) ? "quad4" : "tria3"
                root_cat = is_stress ? "stresses" : "strains"
                
                
                
                
                
                
                id_val = tryparse(Int, parts[1])
                
                
                
                # Usually IDs are integers > 0. Fiber distances contain "." or "E".
                is_id = (id_val !== nothing) && !occursin(".", parts[1])
                
                if is_id && length(parts) >= 9
                    
                    
                    entry = Dict(
                        "fiber_dist" => parse_nastran_float(parts[2]),
                        "normal_x"   => parse_nastran_float(parts[3]),
                        "normal_y"   => parse_nastran_float(parts[4]),
                        "shear_xy"   => parse_nastran_float(parts[5]),
                        "angle"      => parse_nastran_float(parts[6]),
                        "major"      => parse_nastran_float(parts[7]),
                        "minor"      => parse_nastran_float(parts[8]),
                        "von_mises"  => parse_nastran_float(parts[9])
                    )
                    push!(results[root_cat][elem_type], Dict("eid" => id_val, "z1" => entry))

                elseif !is_id && !isempty(results[root_cat][elem_type]) && length(parts) >= 8
                    
                    
                    entry = Dict(
                        "fiber_dist" => parse_nastran_float(parts[1]),
                        "normal_x"   => parse_nastran_float(parts[2]),
                        "normal_y"   => parse_nastran_float(parts[3]),
                        "shear_xy"   => parse_nastran_float(parts[4]),
                        "angle"      => parse_nastran_float(parts[5]),
                        "major"      => parse_nastran_float(parts[6]),
                        "minor"      => parse_nastran_float(parts[7]),
                        "von_mises"  => parse_nastran_float(parts[8])
                    )
                    results[root_cat][elem_type][end]["z2"] = entry
                end
            end

        catch e
            
        end
    end

    return results
end

function convert_f06(f06_path::String, output_path::String)
    if !isfile(f06_path); println("ERROR: File not found: $f06_path"); return nothing; end

    println(">>> Reading F06: $f06_path")
    file_content = read(f06_path, String)

    println("... Parsing content")
    json_data = parse_f06_to_json(file_content)

    println("... Extracted $(length(json_data["displacements"])) displacement entries")
    println("... Extracted $(length(json_data["spc_forces"])) spc force entries")
    println("... Extracted $(length(json_data["forces"]["cbar"])) cbar force entries")
    println("... Extracted $(length(json_data["forces"]["crod"])) crod force entries")
    println("... Extracted $(length(json_data["stresses"]["cbar"])) cbar stress entries")
    println("... Extracted $(length(json_data["stresses"]["crod"])) crod stress entries")
    println("... Extracted $(length(json_data["stresses"]["quad4"])) quad4 stress entries")
    println("... Extracted $(length(json_data["strains"]["cbar"])) cbar strain entries")
    println("... Extracted $(length(json_data["strains"]["crod"])) crod strain entries")

    out_dir = dirname(output_path)
    if !isdir(out_dir); mkpath(out_dir); end

    println("... Writing JSON to: $output_path")
    open(output_path, "w") do f
        JSON.print(f, json_data, 4)
    end
    println(">>> Success!")
    return json_data
end

if abspath(PROGRAM_FILE) == @__FILE__
    model_name = length(ARGS) >= 1 ? ARGS[1] : "ala3"
    script_dir = @__DIR__
    f06_path = normpath(joinpath(script_dir, "..", "..", "references", "$(model_name).f06"))
    if !isfile(f06_path)
        f06_path = normpath(joinpath(script_dir, "..", "..", "$(model_name).f06"))
    end
    out_path = normpath(joinpath(script_dir, "..", "..", "references", "$(model_name).NAST.json"))
    convert_f06(f06_path, out_path)
end
