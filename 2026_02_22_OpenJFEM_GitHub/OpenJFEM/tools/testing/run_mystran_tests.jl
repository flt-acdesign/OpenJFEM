
# Automated test runner: runs Julia solver on MYSTRAN test cases and compares with NASTRAN reference f06

using JSON
using Printf

# Include the f06 parser
include(joinpath(@__DIR__, "..", "converters", "f06_2_json.jl"))

"""
Run a single test case:
1. Run Julia solver on .bdf -> produces .JU.JSON
2. Parse NASTRAN reference .f06 -> produces .NAST.json
3. Compare displacements, forces, stresses
"""
function run_single_test(bdf_path::String, f06_path::String, output_dir::String)
    model_name = replace(basename(bdf_path), ".bdf" => "")

    println("\n" * "="^80)
    println("TEST: $model_name")
    println("="^80)

    # Step 1: Run Julia solver
    julia_json = joinpath(output_dir, "$(model_name).JU.JSON")
    println(">>> Running Julia solver on: $bdf_path")

    solver_script = joinpath(@__DIR__, "..", "..", "src", "main.jl")
    cmd = `julia --project=$(joinpath(@__DIR__, "..", "..")) $solver_script $bdf_path $output_dir`
    try
        run(cmd)
    catch e
        println("  ERROR running solver: $e")
        return Dict("status" => "SOLVER_ERROR", "error" => string(e))
    end

    if !isfile(julia_json)
        println("  ERROR: Julia output not found: $julia_json")
        return Dict("status" => "NO_OUTPUT")
    end

    # Step 2: Parse NASTRAN f06
    nast_json_path = joinpath(output_dir, "$(model_name).NAST.json")
    println(">>> Parsing NASTRAN reference: $f06_path")
    nast_data = convert_f06(f06_path, nast_json_path)
    if isnothing(nast_data)
        println("  ERROR: Could not parse NASTRAN f06")
        return Dict("status" => "F06_PARSE_ERROR")
    end

    # Step 3: Load Julia results
    julia_data = Dict{String,Any}(JSON.parsefile(julia_json))

    # Step 4: Compare
    results = compare_results(julia_data, nast_data, model_name)
    return results
end

"""
Compare Julia solver output with NASTRAN reference, return metrics dict.
"""
function compare_results(julia::Dict, nastran::Dict, model_name::String)
    results = Dict(
        "model" => model_name,
        "status" => "OK",
        "displacements" => Dict(),
        "spc_forces" => Dict(),
        "cbar_forces" => Dict(),
        "cbar_stresses" => Dict(),
        "quad4_stresses" => Dict(),
        "tria3_stresses" => Dict()
    )

    # Compare displacements
    println("\n--- DISPLACEMENT COMPARISON ---")
    nast_disp = nastran["displacements"]
    julia_disp = julia["displacements"]

    # Build dicts by grid_id
    nd = Dict{Int, Dict}()
    for d in nast_disp; nd[d["grid_id"]] = d; end
    jd = Dict{Int, Dict}()
    for d in julia_disp; jd[d["grid_id"]] = d; end

    common_grids = intersect(keys(nd), keys(jd))
    println("  NASTRAN grids: $(length(nd)), Julia grids: $(length(jd)), Common: $(length(common_grids))")

    if !isempty(common_grids)
        max_err = 0.0
        max_err_grid = 0
        max_err_dof = ""
        total_err = 0.0
        n_compared = 0
        errors_by_grid = Dict{Int, Float64}()

        for gid in common_grids
            n_d = nd[gid]; j_d = jd[gid]
            for dof in ["t1", "t2", "t3", "r1", "r2", "r3"]
                nv = Float64(get(n_d, dof, 0.0))
                jv = Float64(get(j_d, dof, 0.0))
                abs_err = abs(nv - jv)
                ref = max(abs(nv), 1e-30)
                rel_err = abs_err / ref

                if abs(nv) > 1e-15  # only count meaningful values
                    n_compared += 1
                    total_err += rel_err
                    if rel_err > max_err
                        max_err = rel_err
                        max_err_grid = gid
                        max_err_dof = dof
                    end
                    errors_by_grid[gid] = max(get(errors_by_grid, gid, 0.0), rel_err)
                end
            end
        end

        avg_err = n_compared > 0 ? total_err / n_compared : 0.0

        results["displacements"] = Dict(
            "n_nastran" => length(nd),
            "n_julia" => length(jd),
            "n_common" => length(common_grids),
            "n_compared" => n_compared,
            "max_rel_error" => max_err,
            "max_err_grid" => max_err_grid,
            "max_err_dof" => max_err_dof,
            "avg_rel_error" => avg_err
        )

        if max_err < 0.01
            println("  PASS - Max relative error: $(@sprintf("%.6e", max_err)) at grid $max_err_grid ($max_err_dof)")
        elseif max_err < 0.05
            println("  WARN - Max relative error: $(@sprintf("%.6e", max_err)) at grid $max_err_grid ($max_err_dof)")
        else
            println("  FAIL - Max relative error: $(@sprintf("%.6e", max_err)) at grid $max_err_grid ($max_err_dof)")
        end
        println("  Avg relative error: $(@sprintf("%.6e", avg_err)) over $n_compared DOF comparisons")

        # Show specific values for worst grid
        if max_err_grid > 0 && haskey(nd, max_err_grid) && haskey(jd, max_err_grid)
            println("  Worst grid $max_err_grid:")
            for dof in ["t1", "t2", "t3", "r1", "r2", "r3"]
                nv = Float64(get(nd[max_err_grid], dof, 0.0))
                jv = Float64(get(jd[max_err_grid], dof, 0.0))
                if abs(nv) > 1e-15 || abs(jv) > 1e-15
                    println("    $dof: NASTRAN=$(@sprintf("%.6e", nv))  Julia=$(@sprintf("%.6e", jv))  err=$(@sprintf("%.2f%%", abs(nv-jv)/max(abs(nv),1e-30)*100))")
                end
            end
        end
    end

    # Compare CBAR forces
    println("\n--- CBAR FORCE COMPARISON ---")
    nf = nastran["forces"]["cbar"]
    jf = julia["forces"]["cbar"]
    println("  NASTRAN: $(length(nf)), Julia: $(length(jf))")

    if !isempty(nf) && !isempty(jf)
        nf_d = Dict{Int, Dict}()
        for f in nf; nf_d[f["eid"]] = f; end
        jf_d = Dict{Int, Dict}()
        for f in jf; jf_d[f["eid"]] = f; end

        common_eids = intersect(keys(nf_d), keys(jf_d))
        max_force_err = 0.0

        for eid in common_eids
            for key in ["moment_a1", "moment_a2", "moment_b1", "moment_b2", "shear_1", "shear_2", "axial", "torque"]
                nv = Float64(get(nf_d[eid], key, 0.0))
                jv = Float64(get(jf_d[eid], key, 0.0))
                if abs(nv) > 1e-10
                    rel = abs(nv - jv) / abs(nv)
                    max_force_err = max(max_force_err, rel)
                end
            end
        end

        results["cbar_forces"]["max_rel_error"] = max_force_err
        results["cbar_forces"]["n_common"] = length(common_eids)

        if max_force_err < 0.01
            println("  PASS - Max relative error: $(@sprintf("%.6e", max_force_err))")
        else
            println("  FAIL - Max relative error: $(@sprintf("%.6e", max_force_err))")
        end
    end

    # Compare CROD forces
    println("\n--- CROD FORCE COMPARISON ---")
    nf_rod = get(get(nastran, "forces", Dict()), "crod", [])
    jf_rod = get(get(julia, "forces", Dict()), "crod", [])
    println("  NASTRAN: $(length(nf_rod)), Julia: $(length(jf_rod))")

    if !isempty(nf_rod) && !isempty(jf_rod)
        nfr_d = Dict{Int, Dict}()
        for f in nf_rod; nfr_d[f["eid"]] = f; end
        jfr_d = Dict{Int, Dict}()
        for f in jf_rod; jfr_d[f["eid"]] = f; end

        common_eids = intersect(keys(nfr_d), keys(jfr_d))
        max_rod_err = 0.0

        for eid in common_eids
            for key in ["axial", "torque"]
                nv = Float64(get(nfr_d[eid], key, 0.0))
                jv = Float64(get(jfr_d[eid], key, 0.0))
                if abs(nv) > 1e-10
                    rel = abs(nv - jv) / abs(nv)
                    max_rod_err = max(max_rod_err, rel)
                end
            end
        end

        results["crod_forces"] = Dict("max_rel_error" => max_rod_err, "n_common" => length(common_eids))
        if max_rod_err < 0.01
            println("  PASS - Max relative error: $(@sprintf("%.6e", max_rod_err))")
        else
            println("  FAIL - Max relative error: $(@sprintf("%.6e", max_rod_err))")
        end
    end

    # Compare CBAR stresses
    println("\n--- CBAR STRESS COMPARISON ---")
    ns = nastran["stresses"]["cbar"]
    js = julia["stresses"]["cbar"]
    println("  NASTRAN: $(length(ns)), Julia: $(length(js))")

    if !isempty(ns) && !isempty(js)
        ns_d = Dict{Int, Dict}()
        for s in ns; ns_d[s["eid"]] = s; end
        js_d = Dict{Int, Dict}()
        for s in js; js_d[s["eid"]] = s; end

        common_eids = intersect(keys(ns_d), keys(js_d))
        max_stress_err = 0.0

        for eid in common_eids
            for end_key in ["end_a", "end_b"]
                if haskey(ns_d[eid], end_key) && haskey(js_d[eid], end_key)
                    for pt in ["p1", "p2", "p3", "p4"]
                        nv = Float64(get(ns_d[eid][end_key], pt, 0.0))
                        jv = Float64(get(js_d[eid][end_key], pt, 0.0))
                        if abs(nv) > 1e-10
                            rel = abs(nv - jv) / abs(nv)
                            max_stress_err = max(max_stress_err, rel)
                        end
                    end
                end
            end
        end

        results["cbar_stresses"]["max_rel_error"] = max_stress_err
        results["cbar_stresses"]["n_common"] = length(common_eids)

        if max_stress_err < 0.01
            println("  PASS - Max relative error: $(@sprintf("%.6e", max_stress_err))")
        else
            println("  FAIL - Max relative error: $(@sprintf("%.6e", max_stress_err))")
        end
    end

    # Compare QUAD4 stresses
    println("\n--- QUAD4 STRESS COMPARISON ---")
    nsq = nastran["stresses"]["quad4"]
    jsq = julia["stresses"]["quad4"]
    println("  NASTRAN: $(length(nsq)), Julia: $(length(jsq))")

    if !isempty(nsq) && !isempty(jsq)
        nsq_d = Dict{Int, Dict}()
        for s in nsq; nsq_d[s["eid"]] = s; end
        jsq_d = Dict{Int, Dict}()
        for s in jsq; jsq_d[s["eid"]] = s; end

        common_eids = intersect(keys(nsq_d), keys(jsq_d))
        max_stress_err = 0.0

        for eid in common_eids
            for z_key in ["z1", "z2"]
                if haskey(nsq_d[eid], z_key) && haskey(jsq_d[eid], z_key)
                    nv = Float64(get(nsq_d[eid][z_key], "von_mises", 0.0))
                    jv = Float64(get(jsq_d[eid][z_key], "von_mises", 0.0))
                    if abs(nv) > 1e-10
                        rel = abs(nv - jv) / abs(nv)
                        max_stress_err = max(max_stress_err, rel)
                    end
                end
            end
        end

        results["quad4_stresses"]["max_rel_error"] = max_stress_err
        results["quad4_stresses"]["n_common"] = length(common_eids)

        if max_stress_err < 0.05
            println("  PASS - Max relative error: $(@sprintf("%.6e", max_stress_err))")
        else
            println("  FAIL - Max relative error: $(@sprintf("%.6e", max_stress_err))")
        end
    end

    # Compare SPC forces
    println("\n--- SPC FORCE COMPARISON ---")
    nspc = nastran["spc_forces"]
    jspc = julia["spc_forces"]
    println("  NASTRAN: $(length(nspc)), Julia: $(length(jspc))")

    if !isempty(nspc) && !isempty(jspc)
        nspc_d = Dict{Int, Dict}()
        for s in nspc; nspc_d[s["grid_id"]] = s; end
        jspc_d = Dict{Int, Dict}()
        for s in jspc; jspc_d[s["grid_id"]] = s; end

        common = intersect(keys(nspc_d), keys(jspc_d))
        max_spc_err = 0.0

        for gid in common
            for dof in ["t1", "t2", "t3", "r1", "r2", "r3"]
                nv = Float64(get(nspc_d[gid], dof, 0.0))
                jv = Float64(get(jspc_d[gid], dof, 0.0))
                if abs(nv) > 1e-10
                    rel = abs(nv - jv) / abs(nv)
                    max_spc_err = max(max_spc_err, rel)
                end
            end
        end

        results["spc_forces"]["max_rel_error"] = max_spc_err
        results["spc_forces"]["n_common"] = length(common)

        if max_spc_err < 0.01
            println("  PASS - Max relative error: $(@sprintf("%.6e", max_spc_err))")
        else
            println("  FAIL - Max relative error: $(@sprintf("%.6e", max_spc_err))")
            # Show details
            for gid in common
                for dof in ["t1", "t2", "t3", "r1", "r2", "r3"]
                    nv = Float64(get(nspc_d[gid], dof, 0.0))
                    jv = Float64(get(jspc_d[gid], dof, 0.0))
                    if abs(nv) > 1e-10
                        println("    Grid $gid $dof: NASTRAN=$(@sprintf("%.6e", nv))  Julia=$(@sprintf("%.6e", jv))")
                    end
                end
            end
        end
    end

    return results
end

function main()
    script_dir = @__DIR__
    mystran_dir = normpath(joinpath(script_dir, "..", "..", "MYSTRAN"))
    msc_dir = joinpath(mystran_dir, "MSC")
    fem_input_dir = normpath(joinpath(script_dir, "..", "..", "models"))
    fem_output_dir = normpath(joinpath(script_dir, "..", "..", "output"))
    ref_dir = normpath(joinpath(script_dir, "..", "..", "references"))
    output_dir = normpath(joinpath(script_dir, "..", "..", "output", "MYSTRAN_tests"))

    if !isdir(output_dir); mkpath(output_dir); end

    # Build test list: (bdf_path, f06_path, output_dir)
    test_list = Tuple{String, String, String}[]

    if !isempty(ARGS) && ARGS[1] == "all"
        # Run ALL tests from both directories
        # 1. MYSTRAN tests
        for tc in ["bar_tube", "bar_tube2", "BAR-I12", "bar_static_large", "bar_tube_dollar",
                    "cquad4_pshell_center", "ctria3_pshell_center", "cquad4_bad_quality", "cquad4_pcomp"]
            bdf = joinpath(mystran_dir, "$(tc).bdf")
            f06_name = uppercase(tc)
            f06 = joinpath(msc_dir, "$(f06_name).f06")
            if !isfile(f06); f06 = joinpath(msc_dir, "$(tc).f06"); end
            if isfile(bdf) && isfile(f06)
                push!(test_list, (bdf, f06, output_dir))
            end
        end
        # 2. FEM_input tests
        for bdf_file in readdir(fem_input_dir; join=true)
            if !endswith(bdf_file, ".bdf"); continue; end
            model_name = replace(basename(bdf_file), ".bdf" => "")
            f06 = joinpath(ref_dir, "$(model_name).f06")
            if isfile(f06)
                push!(test_list, (bdf_file, f06, joinpath(fem_output_dir, "test_results")))
            end
        end
    elseif !isempty(ARGS)
        # Single test by name
        tc = ARGS[1]
        bdf = joinpath(mystran_dir, "$(tc).bdf")
        if !isfile(bdf); bdf = joinpath(fem_input_dir, "$(tc).bdf"); end
        f06_name = uppercase(tc)
        f06 = joinpath(msc_dir, "$(f06_name).f06")
        if !isfile(f06); f06 = joinpath(msc_dir, "$(tc).f06"); end
        if !isfile(f06); f06 = joinpath(ref_dir, "$(tc).f06"); end
        if isfile(bdf) && isfile(f06)
            out = startswith(bdf, fem_input_dir) ? joinpath(fem_output_dir, "test_results") : output_dir
            push!(test_list, (bdf, f06, out))
        else
            println("ERROR: Could not find BDF ($bdf) or F06 ($f06)")
            return
        end
    else
        # Default: run MYSTRAN tests only
        for tc in ["bar_tube", "bar_tube2", "BAR-I12", "bar_static_large", "bar_tube_dollar",
                    "cquad4_pshell_center", "ctria3_pshell_center", "cquad4_bad_quality", "cquad4_pcomp"]
            bdf = joinpath(mystran_dir, "$(tc).bdf")
            f06_name = uppercase(tc)
            f06 = joinpath(msc_dir, "$(f06_name).f06")
            if !isfile(f06); f06 = joinpath(msc_dir, "$(tc).f06"); end
            if isfile(bdf) && isfile(f06)
                push!(test_list, (bdf, f06, output_dir))
            end
        end
    end

    all_results = Dict[]

    for (bdf_path, f06_path, out_dir) in test_list
        if !isdir(out_dir); mkpath(out_dir); end
        result = run_single_test(bdf_path, f06_path, out_dir)
        push!(all_results, result)
    end

    # Print summary
    println("\n" * "="^80)
    println("SUMMARY")
    println("="^80)
    @printf("%-35s %-12s %-20s %-20s\n", "Test Case", "Status", "Disp Max Err", "Force Max Err")
    println("-"^90)

    n_pass = 0; n_fail = 0; n_na = 0
    for r in all_results
        name = get(r, "model", "?")
        status = get(r, "status", "?")
        disp_err = get(get(r, "displacements", Dict()), "max_rel_error", NaN)
        force_err = get(get(r, "cbar_forces", Dict()), "max_rel_error", NaN)

        disp_str = isnan(disp_err) ? "N/A" : @sprintf("%.2e (%.1f%%)", disp_err, disp_err*100)
        force_str = isnan(force_err) ? "N/A" : @sprintf("%.2e (%.1f%%)", force_err, force_err*100)

        pass = !isnan(disp_err) && disp_err < 0.05
        if !isnan(disp_err)
            if pass; n_pass += 1; else; n_fail += 1; end
        else
            n_na += 1
        end

        @printf("%-35s %-12s %-20s %-20s\n", name, pass ? "PASS" : (isnan(disp_err) ? "N/A" : "FAIL"), disp_str, force_str)
    end

    println("-"^90)
    println("Total: $(length(all_results)) tests | PASS: $n_pass | FAIL: $n_fail | N/A: $n_na")
end

main()
