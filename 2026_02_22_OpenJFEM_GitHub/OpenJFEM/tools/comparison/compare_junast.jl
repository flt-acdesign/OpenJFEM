#!/usr/bin/env julia
# Compare OpenJFEM Julia vs MSC NASTRAN displacements (standalone)

using JSON

function parse_f06_singularity(f06_path)
    """Parse GRID POINT SINGULARITY TABLE from MSC NASTRAN F06.
    These are DOFs moved to the SB set (constrained due to singularity)."""
    singular = Dict{Int,Set{Int}}()
    lines = readlines(f06_path)
    local in_sing = false
    for line in lines
        if occursin("S I N G U L A R I T Y   T A B L E", line)
            in_sing = true
            continue
        end
        # Skip sub-headers
        if in_sing && (occursin("POINT", line) && occursin("TYPE", line))
            continue
        end
        if in_sing && (occursin("ID", line) && occursin("DIRECTION", line))
            continue
        end
        # Stop at page breaks or new sections
        if in_sing && (occursin("D I S P L A C E", line) || occursin("OLOAD", line) ||
                       occursin("S T R E S S", line) || occursin("F O R C E", line) ||
                       (startswith(strip(line), "1") && length(strip(line)) > 2 && !occursin("G", line)))
            in_sing = false
            continue
        end
        if in_sing && occursin("SB", line)
            parts = split(strip(line))
            if length(parts) >= 3
                try
                    gid = parse(Int, parts[1])
                    dof = parse(Int, parts[3])
                    if !haskey(singular, gid); singular[gid] = Set{Int}(); end
                    push!(singular[gid], dof)
                catch; end
            end
        end
    end
    return singular
end

function parse_f06_displacements(f06_path)
    disp = Dict{Int, Vector{Float64}}()
    lines = readlines(f06_path)
    local in_disp = false
    for line in lines
        if occursin("D I S P L A C E M E N T", line)
            in_disp = true; continue
        end
        if in_disp && (occursin("---", line) || occursin("MAX*", line) ||
                       occursin("OUTPUT FOR", line) || occursin("S P C", line) ||
                       occursin("E L E M E N T", line) || occursin("OLOAD", line) ||
                       occursin("RESULTANT", line) || occursin("SUBCASE", line))
            in_disp = false
            continue
        end
        if in_disp
            parts = split(strip(line))
            if length(parts) >= 8
                try
                    gid = parse(Int, parts[1])
                    vals = [parse(Float64, parts[k]) for k in 3:8]
                    disp[gid] = vals
                catch; end
            end
        end
    end
    return disp
end

function main()
    f06_path = joinpath(@__DIR__, "..", "..", "references", "OpenJFEM.f06")
    json_path = joinpath(@__DIR__, "..", "..", "output", "OpenJFEM.JU.JSON")

    println("Parsing MSC NASTRAN F06...")
    msc_disp = parse_f06_displacements(f06_path)
    autospc = parse_f06_singularity(f06_path)
    println("  MSC grids: $(length(msc_disp))")

    # Count AUTOSPC DOFs
    total_autospc = sum(length(v) for (_, v) in autospc; init=0)
    println("  Singular (SB) grids: $(length(autospc)), total DOFs: $(total_autospc)")
    dof_names = ["T1","T2","T3","R1","R2","R3"]
    dof_counts = zeros(Int, 6)
    for (_, dofs) in autospc
        for d in dofs
            if d >= 1 && d <= 6; dof_counts[d] += 1; end
        end
    end
    for k in 1:6
        if dof_counts[k] > 0; println("    $(dof_names[k]): $(dof_counts[k])"); end
    end

    println("Parsing Julia JSON...")
    jul_data = JSON.parsefile(json_path)
    jul_disp = Dict{Int, Vector{Float64}}()
    for d in jul_data["displacements"]
        jul_disp[d["grid_id"]] = [d["t1"], d["t2"], d["t3"], d["r1"], d["r2"], d["r3"]]
    end
    println("  Julia grids: $(length(jul_disp))")

    common = sort(collect(intersect(keys(msc_disp), keys(jul_disp))))
    println("  Common grids: $(length(common))")

    # Build set of AUTOSPC'd (gid, dof) pairs
    auto_set = Set{Tuple{Int,Int}}()
    for (gid, dofs) in autospc
        for d in dofs; push!(auto_set, (gid, d)); end
    end

    # IN-PLANE (T1, T2) comparison
    println("\n" * "="^70)
    println("IN-PLANE (T1, T2) COMPARISON")
    println("="^70)

    errors_noauto = Float64[]
    errors_all = Float64[]
    big_errors = Tuple{Int,Int,Float64,Float64,Float64}[]

    for gid in common
        for k in 1:2
            mv, jv = msc_disp[gid][k], jul_disp[gid][k]
            if abs(mv) < 1e-10 && abs(jv) < 1e-10; continue; end
            ref = max(abs(mv), 1e-10)
            err = abs(mv - jv) / ref * 100
            push!(errors_all, err)

            is_auto = (gid, k) in auto_set
            if !is_auto
                push!(errors_noauto, err)
                if err > 20.0
                    push!(big_errors, (gid, k, mv, jv, err))
                end
            end
        end
    end

    sort!(errors_all); sort!(errors_noauto)

    n = length(errors_all)
    println("\nIncluding AUTOSPC'd DOFs ($n values):")
    println("  Median: $(round(errors_all[n÷2], digits=2))%")
    println("  90th:   $(round(errors_all[Int(ceil(0.9*n))], digits=2))%")
    println("  Max:    $(round(errors_all[end], digits=2))%")

    n2 = length(errors_noauto)
    println("\nExcluding AUTOSPC'd DOFs ($n2 values):")
    println("  Median: $(round(errors_noauto[n2÷2], digits=2))%")
    println("  90th:   $(round(errors_noauto[Int(ceil(0.9*n2))], digits=2))%")
    println("  95th:   $(round(errors_noauto[Int(ceil(0.95*n2))], digits=2))%")
    println("  99th:   $(round(errors_noauto[Int(ceil(0.99*n2))], digits=2))%")
    println("  Max:    $(round(errors_noauto[end], digits=2))%")
    println("  <1%:  $(count(e->e<1, errors_noauto)) / $n2 ($(round(100*count(e->e<1, errors_noauto)/n2, digits=1))%)")
    println("  <5%:  $(count(e->e<5, errors_noauto)) / $n2 ($(round(100*count(e->e<5, errors_noauto)/n2, digits=1))%)")
    println("  <10%: $(count(e->e<10, errors_noauto)) / $n2 ($(round(100*count(e->e<10, errors_noauto)/n2, digits=1))%)")
    println("  <20%: $(count(e->e<20, errors_noauto)) / $n2 ($(round(100*count(e->e<20, errors_noauto)/n2, digits=1))%)")

    if !isempty(big_errors)
        sort!(big_errors, by=x->x[5], rev=true)
        println("\nTop 20 non-AUTOSPC errors >20%:")
        for i in 1:min(20, length(big_errors))
            gid, k, mv, jv, err = big_errors[i]
            println("  Grid $(gid) $(dof_names[k]): MSC=$(round(mv,sigdigits=5))  JUL=$(round(jv,sigdigits=5))  err=$(round(err,digits=1))%")
        end
    end

    # DISPLACEMENT MAGNITUDE comparison
    println("\n" * "="^70)
    println("DISPLACEMENT MAGNITUDE (T1+T2+T3) - excl. AUTOSPC grids")
    println("="^70)

    local msc_max_mag = 0.0
    local msc_max_gid = 0
    local jul_max_mag = 0.0
    local jul_max_gid = 0
    mag_errors = Float64[]

    for gid in common
        # Skip grids with any translation AUTOSPC
        has_auto = false
        for k in 1:3
            if (gid, k) in auto_set; has_auto = true; break; end
        end
        if has_auto; continue; end

        mv = sqrt(sum(msc_disp[gid][k]^2 for k in 1:3))
        jv = sqrt(sum(jul_disp[gid][k]^2 for k in 1:3))
        if mv > msc_max_mag; msc_max_mag = mv; msc_max_gid = gid; end
        if jv > jul_max_mag; jul_max_mag = jv; jul_max_gid = gid; end
        if mv > 1e-8
            push!(mag_errors, abs(mv-jv)/mv * 100)
        end
    end
    sort!(mag_errors)
    nm = length(mag_errors)
    println("MSC max disp: $(round(msc_max_mag, sigdigits=6)) at Grid $(msc_max_gid)")
    println("JUL max disp: $(round(jul_max_mag, sigdigits=6)) at Grid $(jul_max_gid)")
    if nm > 0
        println("Mag error median: $(round(mag_errors[nm÷2], digits=2))%")
        println("Mag error 90th:   $(round(mag_errors[Int(ceil(0.9*nm))], digits=2))%")
        println("Mag error 95th:   $(round(mag_errors[Int(ceil(0.95*nm))], digits=2))%")
        println("Mag error 99th:   $(round(mag_errors[Int(ceil(0.99*nm))], digits=2))%")
        println("Mag error max:    $(round(mag_errors[end], digits=2))%")
    end

    # Details at max displacement grid
    if msc_max_gid > 0 && haskey(jul_disp, msc_max_gid)
        m = msc_disp[msc_max_gid]; j = jul_disp[msc_max_gid]
        println("\nAt MSC max-displacement Grid $(msc_max_gid):")
        for k in 1:6
            if abs(m[k]) > 1e-10 || abs(j[k]) > 1e-10
                ref = max(abs(m[k]), 1e-10)
                err = abs(m[k]-j[k])/ref * 100
                println("  $(dof_names[k]): MSC=$(round(m[k],sigdigits=6))  JUL=$(round(j[k],sigdigits=6))  err=$(round(err,digits=2))%")
            end
        end
    end
end

main()
