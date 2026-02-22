import JSON
using Statistics
using LinearAlgebra

function compare(model_name)
    ref_file = joinpath(@__DIR__, "..", "..", "references", "$(model_name).NAST.json")
    test_file = joinpath(@__DIR__, "..", "..", "output", "$(model_name).JU.JSON")

    if !isfile(ref_file)
        # Try parsing the f06
        include(joinpath(@__DIR__, "..", "converters", "f06_2_json.jl"))
        f06_file = joinpath(@__DIR__, "..", "..", "references", "$(model_name).f06")
        if isfile(f06_file)
            convert_f06(f06_file, ref_file)
        else
            println("ERROR: No NAST.json or f06 for $model_name")
            return
        end
    end

    nast = JSON.parsefile(ref_file)
    ju   = JSON.parsefile(test_file)

    nd = Dict(d["grid_id"]=>d for d in nast["displacements"])
    jd = Dict(d["grid_id"]=>d for d in ju["displacements"])
    common = intersect(keys(nd), keys(jd))
    println("Model: $model_name")
    println("  NAST: $(length(nd)) grids, JU: $(length(jd)) grids, Common: $(length(common))")

    for dof in ["t1","t2","t3"]
        nvals = Float64[Float64(get(nd[g], dof, 0.0)) for g in common]
        jvals = Float64[Float64(get(jd[g], dof, 0.0)) for g in common]
        l2n = norm(nvals); l2j = norm(jvals)
        corr = (std(nvals) > 0 && std(jvals) > 0) ? cor(nvals, jvals) : 0.0
        ratio = l2n > 1e-30 ? l2j/l2n : 0.0
        println("  $dof: L2_ratio=$(round(ratio, digits=4)) corr=$(round(corr, digits=4))")
    end

    # Force comparison for quad4
    for ftype in ["quad4", "tria3", "cbar"]
        nf = get(get(nast, "forces", Dict()), ftype, [])
        jf = get(get(ju, "forces", Dict()), ftype, [])
        if isempty(nf) || isempty(jf); continue; end
        nfd = Dict(f["eid"]=>f for f in nf)
        jfd = Dict(f["eid"]=>f for f in jf)
        ce = intersect(keys(nfd), keys(jfd))
        if isempty(ce); continue; end

        fields = ftype == "cbar" ? ["axial"] : ["fx", "fy", "fxy"]
        for field in fields
            nvals = Float64[Float64(get(nfd[e], field, 0.0)) for e in ce]
            jvals = Float64[Float64(get(jfd[e], field, 0.0)) for e in ce]
            l2n = norm(nvals); l2j = norm(jvals)
            ratio = l2n > 1e-30 ? l2j/l2n : 0.0
            corr = (std(nvals) > 0 && std(jvals) > 0) ? cor(nvals, jvals) : 0.0
            println("  $ftype.$field: L2_ratio=$(round(ratio, digits=4)) corr=$(round(corr, digits=4))")
        end
    end
end

model = length(ARGS) >= 1 ? ARGS[1] : "fwing"
compare(model)
