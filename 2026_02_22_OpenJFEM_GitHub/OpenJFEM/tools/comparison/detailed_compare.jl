import JSON
using Statistics
using LinearAlgebra
using Printf

function detailed_compare(model_name)
    ref_file = joinpath(@__DIR__, "..", "..", "references", "$(model_name).NAST.json")
    test_file = joinpath(@__DIR__, "..", "..", "output", "$(model_name).JU.JSON")

    nast = JSON.parsefile(ref_file)
    ju   = JSON.parsefile(test_file)

    nd = Dict(d["grid_id"]=>d for d in nast["displacements"])
    jd = Dict(d["grid_id"]=>d for d in ju["displacements"])
    common = sort(collect(intersect(keys(nd), keys(jd))))

    println("Model: $model_name")
    println("  Common grids: $(length(common))")

    # Find max displacement grids in NASTRAN
    println("\n=== TOP 10 NASTRAN DISPLACEMENTS (by |T3|) ===")
    sorted_by_t3 = sort(common, by=g->abs(Float64(get(nd[g],"t3",0.0))), rev=true)
    for g in sorted_by_t3[1:min(10, length(sorted_by_t3))]
        nt = [Float64(get(nd[g],"t$i",0.0)) for i in 1:3]
        jt = [Float64(get(jd[g],"t$i",0.0)) for i in 1:3]
        mag_n = norm(nt); mag_j = norm(jt)
        @printf("  Grid %10d: NAST(%10.3f,%10.3f,%10.3f) JU(%10.3f,%10.3f,%10.3f) |N|=%10.3f |J|=%10.3f ratio=%.4f\n",
                g, nt[1], nt[2], nt[3], jt[1], jt[2], jt[3], mag_n, mag_j, mag_n > 0 ? mag_j/mag_n : 0.0)
    end

    # Find max displacement grids in Julia
    println("\n=== TOP 10 JULIA DISPLACEMENTS (by |T3|) ===")
    sorted_by_jt3 = sort(common, by=g->abs(Float64(get(jd[g],"t3",0.0))), rev=true)
    for g in sorted_by_jt3[1:min(10, length(sorted_by_jt3))]
        nt = [Float64(get(nd[g],"t$i",0.0)) for i in 1:3]
        jt = [Float64(get(jd[g],"t$i",0.0)) for i in 1:3]
        mag_n = norm(nt); mag_j = norm(jt)
        @printf("  Grid %10d: NAST(%10.3f,%10.3f,%10.3f) JU(%10.3f,%10.3f,%10.3f) |N|=%10.3f |J|=%10.3f ratio=%.4f\n",
                g, nt[1], nt[2], nt[3], jt[1], jt[2], jt[3], mag_n, mag_j, mag_n > 0 ? mag_j/mag_n : 0.0)
    end

    # Statistics per DOF
    println("\n=== DISPLACEMENT STATISTICS ===")
    for dof in ["t1","t2","t3"]
        nvals = Float64[Float64(get(nd[g], dof, 0.0)) for g in common]
        jvals = Float64[Float64(get(jd[g], dof, 0.0)) for g in common]
        @printf("  %s: NAST range [%10.3f, %10.3f] mean=%10.3f std=%10.3f\n", dof, minimum(nvals), maximum(nvals), mean(nvals), std(nvals))
        @printf("       JU   range [%10.3f, %10.3f] mean=%10.3f std=%10.3f\n", minimum(jvals), maximum(jvals), mean(jvals), std(jvals))
    end

    # Check RBE3 reference grids specifically
    println("\n=== FORCE APPLICATION GRIDS ===")
    for g in [7922000, 48000774]
        if haskey(nd, g) && haskey(jd, g)
            nt = [Float64(get(nd[g],"t$i",0.0)) for i in 1:3]
            jt = [Float64(get(jd[g],"t$i",0.0)) for i in 1:3]
            @printf("  Grid %10d: NAST(%10.3f,%10.3f,%10.3f) JU(%10.3f,%10.3f,%10.3f)\n",
                    g, nt[1], nt[2], nt[3], jt[1], jt[2], jt[3])
        end
    end

    # Check constrained grid
    println("\n=== CONSTRAINED GRID ===")
    if haskey(nd, 99002281) && haskey(jd, 99002281)
        nt = [Float64(get(nd[99002281],"t$i",0.0)) for i in 1:3]
        jt = [Float64(get(jd[99002281],"t$i",0.0)) for i in 1:3]
        @printf("  Grid 99002281: NAST(%10.6f,%10.6f,%10.6f) JU(%10.6f,%10.6f,%10.6f)\n",
                nt[1], nt[2], nt[3], jt[1], jt[2], jt[3])
    end
end

model = length(ARGS) >= 1 ? ARGS[1] : "OpenJFEM"
detailed_compare(model)
