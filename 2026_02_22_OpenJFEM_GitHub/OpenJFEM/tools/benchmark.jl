"""
    benchmark.jl

Benchmarks MYSTRAN vs Julia solver on ala3, fwing, and OpenJFEM models.
Measures wall-clock time for each solver run.

Usage:
    julia benchmark.jl [mystran|julia|both]
"""

using Dates

const JFEM_ROOT = normpath(joinpath(@__DIR__, ".."))
const MODELS_DIR = joinpath(JFEM_ROOT, "models")
const CONVERTER = joinpath(@__DIR__, "converters", "nastran_to_mystran.jl")
const JULIA_SOLVER = joinpath(JFEM_ROOT, "src", "main.jl")
const MYSTRAN_EXE = normpath(joinpath(JFEM_ROOT, "..", "MYSTRAN", "mystran-17.0.0-windows-x86_64.exe"))

const MODELS = ["ala3", "fwing", "OpenJFEM"]

function run_mystran_timed(model_name::String)
    bdf_path = joinpath(MODELS_DIR, model_name, "$model_name.bdf")
    if !isfile(bdf_path)
        println("  ERROR: BDF not found: $bdf_path")
        return NaN
    end

    println("\n  ▶ MYSTRAN: $model_name")
    println("    Input: $bdf_path")

    t_start = time()
    try
        cmd = `julia --project=$(JFEM_ROOT) $(CONVERTER) $(bdf_path) $(MYSTRAN_EXE)`
        run(pipeline(cmd, stdout=stdout, stderr=stderr))
    catch e
        println("    WARNING: MYSTRAN run error: $e")
    end
    t_elapsed = time() - t_start

    println("    ⏱  MYSTRAN $model_name: $(round(t_elapsed, digits=2)) seconds")
    return t_elapsed
end

function run_julia_timed(model_name::String)
    bdf_path = joinpath(MODELS_DIR, model_name, "$model_name.bdf")
    out_dir = joinpath(MODELS_DIR, model_name)
    if !isfile(bdf_path)
        println("  ERROR: BDF not found: $bdf_path")
        return NaN
    end

    println("\n  ▶ Julia solver: $model_name")
    println("    Input: $bdf_path")

    t_start = time()
    try
        cmd = `julia --project=$(JFEM_ROOT) $(JULIA_SOLVER) $(bdf_path) $(out_dir)`
        run(pipeline(cmd, stdout=stdout, stderr=stderr))
    catch e
        println("    WARNING: Julia solver error: $e")
    end
    t_elapsed = time() - t_start

    println("    ⏱  Julia $model_name: $(round(t_elapsed, digits=2)) seconds")
    return t_elapsed
end

function main()
    mode = length(ARGS) >= 1 ? lowercase(ARGS[1]) : "both"

    println("╔══════════════════════════════════════════════════╗")
    println("║       Solver Benchmark: MYSTRAN vs Julia         ║")
    println("╚══════════════════════════════════════════════════╝")
    println()
    println("  Models:  $(join(MODELS, ", "))")
    println("  Mode:    $mode")
    println("  Time:    $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
    println()

    mystran_times = Dict{String, Float64}()
    julia_times = Dict{String, Float64}()

    if mode in ["mystran", "both"]
        println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        println("  MYSTRAN RUNS")
        println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        for model in MODELS
            mystran_times[model] = run_mystran_timed(model)
        end
    end

    if mode in ["julia", "both"]
        println()
        println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        println("  JULIA SOLVER RUNS")
        println("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        for model in MODELS
            julia_times[model] = run_julia_timed(model)
        end
    end

    # ── Summary ──
    println()
    println("╔══════════════════════════════════════════════════════════════════╗")
    println("║                    BENCHMARK RESULTS                            ║")
    println("╠════════════════╦══════════════╦══════════════╦══════════════════╣")
    println("║ Model          ║ MYSTRAN (s)  ║ Julia (s)    ║ Ratio (M/J)      ║")
    println("╠════════════════╬══════════════╬══════════════╬══════════════════╣")

    for model in MODELS
        mt = get(mystran_times, model, NaN)
        jt = get(julia_times, model, NaN)
        ratio = (isnan(mt) || isnan(jt) || jt == 0) ? "N/A" : string(round(mt / jt, digits=2)) * "x"
        mt_s = isnan(mt) ? "N/A" : string(round(mt, digits=2))
        jt_s = isnan(jt) ? "N/A" : string(round(jt, digits=2))
        println("║ $(rpad(model, 14)) ║ $(rpad(mt_s, 12)) ║ $(rpad(jt_s, 12)) ║ $(rpad(ratio, 16)) ║")
    end

    println("╚════════════════╩══════════════╩══════════════╩══════════════════╝")

    # Also write results to a file
    results_path = joinpath(JFEM_ROOT, "tools", "benchmark_results.txt")
    open(results_path, "w") do io
        println(io, "Solver Benchmark Results")
        println(io, "Date: $(Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "")
        println(io, "Model           | MYSTRAN (s) | Julia (s)   | Ratio (M/J)")
        println(io, "----------------|-------------|-------------|-------------")
        for model in MODELS
            mt = get(mystran_times, model, NaN)
            jt = get(julia_times, model, NaN)
            ratio = (isnan(mt) || isnan(jt) || jt == 0) ? "N/A" : string(round(mt / jt, digits=2)) * "x"
            mt_s = isnan(mt) ? "N/A" : string(round(mt, digits=2))
            jt_s = isnan(jt) ? "N/A" : string(round(jt, digits=2))
            println(io, "$(rpad(model, 16))| $(rpad(mt_s, 12))| $(rpad(jt_s, 12))| $(ratio)")
        end
    end
    println("\nResults saved to: $results_path")
end

main()
