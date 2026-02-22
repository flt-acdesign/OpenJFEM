
using JSON
using Printf
using LinearAlgebra
using Statistics


const REF_DIR = joinpath(@__DIR__, "..", "..", "references")
const OUTPUT_DIR = joinpath(@__DIR__, "..", "..", "output")
const SIGNIFICANCE_THRESHOLD = 0.001 

struct FieldStats
    max_ref::Float64
    max_diff::Float64      
    mean_diff_pct::Float64 
    median_ratio::Float64  
    l2_ratio::Float64      
    sign_err_pct::Float64  
    corr_coeff::Float64    
end

function safe_get(d, k, default)
    return haskey(d, k) ? d[k] : default
end

function calculate_stats(ref_vals, test_vals)
    n = length(ref_vals)
    if n == 0; return FieldStats(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); end

    max_val = maximum(abs.(ref_vals))
    threshold = max_val * 0.001 
    if threshold < SIGNIFICANCE_THRESHOLD; threshold = SIGNIFICANCE_THRESHOLD; end
    
    sig_idx = findall(x -> abs(x) > threshold, ref_vals)
    
    if isempty(sig_idx)
        max_test = maximum(abs.(test_vals))
        if max_test < threshold
            return FieldStats(max_val, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0)
        else
            return FieldStats(max_val, max_test, 100.0, 0.0, 0.0, 0.0, 0.0)
        end
    end

    r = ref_vals[sig_idx]
    t = test_vals[sig_idx]

    abs_diffs = abs.(t .- r)
    max_diff = maximum(abs_diffs)
    errs = abs_diffs ./ abs.(r) .* 100.0
    mean_err = mean(errs)

    ratios = t ./ r
    valid_ratios = filter(x -> !isnan(x) && !isinf(x), ratios)
    med_ratio = isempty(valid_ratios) ? 0.0 : median(valid_ratios)
    
    norm_r = norm(r)
    norm_t = norm(t)
    l2_ratio = (norm_r > 1e-12) ? norm_t / norm_r : 0.0

    sign_mismatches = count(x -> x < 0, r .* t)
    sign_err_pct = (length(r) > 0) ? (sign_mismatches / length(r)) * 100.0 : 0.0

    if std(r) == 0 || std(t) == 0
        corr = 0.0
    else
        corr = cor(r, t)
    end

    return FieldStats(max_val, max_diff, mean_err, med_ratio, l2_ratio, sign_err_pct, corr)
end

function split_data_into_chunks(full_list, id_key)
    if isempty(full_list); return []; end
    
    chunks = []
    current_chunk = Dict{Int, Any}()
    
    for item in full_list
        id = item[id_key]
        if haskey(current_chunk, id)
            push!(chunks, current_chunk)
            current_chunk = Dict{Int, Any}()
        end
        current_chunk[id] = item
    end
    push!(chunks, current_chunk)
    return chunks
end

function analyze_chunk_pair(ref_chunk, test_chunk, fields, label)
    if isempty(ref_chunk) || isempty(test_chunk)
        return
    end
    
    common_ids = intersect(keys(ref_chunk), keys(test_chunk))
    if isempty(common_ids); return; end

    println("\n... $label ...")
    @printf("%-10s | %-9s | %-9s | %-7s | %-7s | %-7s | %-7s | %-6s\n", 
            "Field", "MaxRef", "MaxDiff", "MeanErr", "MedRat", "L2 Rat", "Sgn%", "Corr")
    println("-"^95)

    for f in fields
        r_vals = Float64[]
        t_vals = Float64[]
        for id in common_ids
            push!(r_vals, Float64(safe_get(ref_chunk[id], f, 0.0)))
            push!(t_vals, Float64(safe_get(test_chunk[id], f, 0.0)))
        end

        s = calculate_stats(r_vals, t_vals)
        
        is_bad = (s.mean_diff_pct > 10.0) || (s.corr_coeff < 0.9)
        is_zero = (abs(s.l2_ratio) < 1e-4) && (s.max_ref > 1e-4)
        status = is_zero ? "ZERO" : (is_bad ? "FAIL" : "OK")

        @printf("%-10s | %9.2E | %9.2E | %6.0f%% | %7.3f | %7.3f | %6.0f%% | %6.3f [%s]\n", 
            f, s.max_ref, s.max_diff, s.mean_diff_pct, s.median_ratio, s.l2_ratio, s.sign_err_pct, s.corr_coeff, status)
    end
end

function flatten_shell_stress(chunk)
    """Flatten z1/z2 stress/strain into separate entries for QUAD4/TRIA3."""
    new_chunk = Dict{String, Any}()
    for (id, item) in chunk
        for lay in ["z1", "z2"]
            if haskey(item, lay)
                sub = copy(item[lay])
                sub["_id"] = "$(id)_$lay"
                new_chunk["$(id)_$lay"] = sub
            end
        end
    end
    return new_chunk
end

function flatten_cbar_stress(chunk)
    """Flatten end_a/end_b stress into separate entries for CBAR, keeping axial."""
    new_chunk = Dict{String, Any}()
    for (id, item) in chunk
        axial_val = get(item, "axial", 0.0)
        for endname in ["end_a", "end_b"]
            if haskey(item, endname)
                sub = copy(item[endname])
                sub["axial"] = axial_val
                sub["_id"] = "$(id)_$endname"
                new_chunk["$(id)_$endname"] = sub
            end
        end
    end
    return new_chunk
end

function compare_category(ref_data, test_data, top_key, sub_key, id_key, fields, label; flatten_mode=:auto)
    r_list = isnothing(sub_key) ? get(ref_data, top_key, []) : get(get(ref_data, top_key, Dict()), sub_key, [])
    t_list = isnothing(sub_key) ? get(test_data, top_key, []) : get(get(test_data, top_key, Dict()), sub_key, [])

    if isempty(r_list) || isempty(t_list); return; end

    r_chunks = split_data_into_chunks(r_list, id_key)
    t_chunks = split_data_into_chunks(t_list, id_key)

    println("\n===============================================================================================")
    println(" ANALYSIS: $label")
    println(" Ref entries: $(length(r_list)) | Test entries: $(length(t_list))")
    println(" Ref Subcases: $(length(r_chunks)) | Test Subcases: $(length(t_chunks))")
    println("===============================================================================================")

    n_compare = min(length(r_chunks), length(t_chunks))

    for i in 1:n_compare
        println("\n>>> SUBCASE $i Comparison:")

        r_c = r_chunks[i]
        t_c = t_chunks[i]

        if flatten_mode == :shell_stress
            analyze_chunk_pair(flatten_shell_stress(r_c), flatten_shell_stress(t_c), fields, label)
        elseif flatten_mode == :cbar_stress
            analyze_chunk_pair(flatten_cbar_stress(r_c), flatten_cbar_stress(t_c), fields, label)
        else
            analyze_chunk_pair(r_c, t_c, fields, label)
        end
    end
end

function run_full_comparison(model_name::String)
    
    ref_file = joinpath(REF_DIR, "$(model_name).NAST.json")
    test_file = joinpath(OUTPUT_DIR, "$(model_name).JU.JSON")

    if !isfile(ref_file) || !isfile(test_file)
        println("ERROR: Files not found for model '$model_name'.")
        println("  Ref:  $ref_file")
        println("  Test: $test_file")
        return
    end

    println("Comparing Model: $model_name")
    println("  Ref:  $ref_file")
    println("  Test: $test_file")

    ref = JSON.parsefile(ref_file)
    test = JSON.parsefile(test_file)

    compare_category(ref, test, "displacements", nothing, "grid_id", 
        ["t1", "t2", "t3", "r1", "r2", "r3"], "DISPLACEMENTS")

    compare_category(ref, test, "spc_forces", nothing, "grid_id", 
        ["t1", "t2", "t3", "r1", "r2", "r3"], "SPC FORCES")

    compare_category(ref, test, "forces", "quad4", "eid",
        ["fx", "fy", "fxy", "mx", "my", "mxy", "qx", "qy"], "QUAD4 FORCES")

    compare_category(ref, test, "forces", "tria3", "eid",
        ["fx", "fy", "fxy", "mx", "my", "mxy", "qx", "qy"], "TRIA3 FORCES")

    compare_category(ref, test, "forces", "cbar", "eid",
        ["axial", "shear_1", "shear_2", "torque", "moment_a1", "moment_a2", "moment_b1", "moment_b2"], "CBAR FORCES")

    compare_category(ref, test, "forces", "crod", "eid",
        ["axial", "torque"], "CROD FORCES")

    compare_category(ref, test, "stresses", "quad4", "eid",
        ["normal_x", "normal_y", "shear_xy", "von_mises"], "QUAD4 STRESSES";
        flatten_mode=:shell_stress)

    compare_category(ref, test, "stresses", "tria3", "eid",
        ["normal_x", "normal_y", "shear_xy", "von_mises"], "TRIA3 STRESSES";
        flatten_mode=:shell_stress)

    compare_category(ref, test, "stresses", "cbar", "eid",
        ["axial", "p1", "p2", "p3", "p4"], "CBAR STRESSES";
        flatten_mode=:cbar_stress)

    compare_category(ref, test, "stresses", "crod", "eid",
        ["axial", "torsional"], "CROD STRESSES")

    compare_category(ref, test, "strains", "quad4", "eid",
        ["normal_x", "normal_y", "shear_xy"], "QUAD4 STRAINS";
        flatten_mode=:shell_stress)
end

model_name = length(ARGS) >= 1 ? ARGS[1] : "OpenJFEM"
run_full_comparison(model_name)
