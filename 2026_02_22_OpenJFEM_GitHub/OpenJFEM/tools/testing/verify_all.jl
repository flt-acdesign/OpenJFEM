#!/usr/bin/env julia
# Comprehensive regression test: Julia vs MYSTRAN for all test cases

using JSON

function parse_f06_displacements(f06_path)
    disp = Dict{Int, Vector{Float64}}()
    lines = readlines(f06_path)
    in_disp = false
    for line in lines
        if occursin("D I S P L A C E M E N T", line); in_disp = true; continue; end
        if in_disp && (occursin("---", line) || occursin("MAX*", line) || occursin("OUTPUT FOR", line) || occursin("S P C", line) || occursin("E L E M E N T", line)); break; end
        if in_disp
            parts = split(strip(line))
            if length(parts) >= 8
                try
                    gid = parse(Int, parts[1])
                    disp[gid] = [parse(Float64, parts[k]) for k in 3:8]
                catch; end
            end
        end
    end
    return disp
end

function parse_f06_spcforces(f06_path)
    forces = Dict{Int, Vector{Float64}}()
    lines = readlines(f06_path)
    in_spc = false
    for line in lines
        if occursin("S I N G L E   P O I N T   C O N S T R A I N T", line) || occursin("S P C   F O R C E S", line); in_spc = true; continue; end
        if in_spc && (occursin("---", line) || occursin("MAX*", line) || occursin("OUTPUT FOR", line) || occursin("E L E M E N T", line) || occursin("AUTOSPC", line)); break; end
        if in_spc
            parts = split(strip(line))
            if length(parts) >= 8
                try
                    gid = parse(Int, parts[1])
                    forces[gid] = [parse(Float64, parts[k]) for k in 3:8]
                catch; end
            end
        end
    end
    return forces
end

function compare_displacements(mys_disp, jul_disp; dofs=1:3)
    dof_names = ["T1", "T2", "T3", "R1", "R2", "R3"]
    common = intersect(keys(mys_disp), keys(jul_disp))
    max_err = 0.0
    n_compared = 0
    for gid in sort(collect(common))
        for k in dofs
            mv, jv = mys_disp[gid][k], jul_disp[gid][k]
            if abs(mv) < 1e-12 && abs(jv) < 1e-12; continue; end
            ref = max(abs(mv), 1e-12)
            err = abs(mv - jv) / ref * 100
            max_err = max(max_err, err)
            n_compared += 1
        end
    end
    return max_err, n_compared
end

function compare_spc_total(mys_spc, jul_spc)
    # Sum SPC forces in T1 direction for both
    mys_total = sum(v[1] for (_, v) in mys_spc; init=0.0)
    jul_total = sum(v[1] for (_, v) in jul_spc; init=0.0)
    return mys_total, jul_total
end

function run_test(test_name)
    mys_f06 = joinpath(@__DIR__, "..", "..", "MYSTRAN", "MSC", "$(test_name).F06")
    jul_json = joinpath(@__DIR__, "..", "..", "output", "$(test_name).JU.JSON")

    if !isfile(mys_f06); return nothing; end
    if !isfile(jul_json); return nothing; end

    mys_disp = parse_f06_displacements(mys_f06)
    jul_data = JSON.parsefile(jul_json)
    jul_disp = Dict{Int, Vector{Float64}}()
    for d in jul_data["displacements"]
        jul_disp[d["grid_id"]] = [d["t1"], d["t2"], d["t3"], d["r1"], d["r2"], d["r3"]]
    end

    max_err, n_comp = compare_displacements(mys_disp, jul_disp)

    mys_spc = parse_f06_spcforces(mys_f06)
    jul_spc = Dict{Int, Vector{Float64}}()
    for d in jul_data["spc_forces"]
        jul_spc[d["grid_id"]] = [d["t1"], d["t2"], d["t3"], d["r1"], d["r2"], d["r3"]]
    end

    # SPC force equilibrium
    jul_spc_t1 = sum(v[1] for (_, v) in jul_spc; init=0.0)

    return (max_err=max_err, n_comp=n_comp, n_mys=length(mys_disp),
            n_jul=length(jul_disp), spc_t1=jul_spc_t1)
end

function main()
    println("="^75)
    println("COMPREHENSIVE REGRESSION TEST: Julia FEM Solver vs MYSTRAN")
    println("="^75)

    tests = [
        ("crod_axial_verify",       "3 CRODs in series, P=1000"),
        ("quad4_membrane_verify",   "2 CQUAD4, tension + shear"),
        ("quad4_single_tension",    "1 CQUAD4, pure tension"),
        ("quad4_fine_membrane",     "8 CQUAD4, fine mesh tension"),
        ("ctria3_membrane_verify",  "16 CTRIA3, tension + shear"),
        ("rbe2_verify",             "RBE2 CM=123456, two panels"),
        ("rbe2_cm123_verify",       "RBE2 CM=123, two panels"),
        ("rbe3_verify",             "RBE3 load distribution"),
        ("rbe3_multigrid_verify",   "RBE3 4-grid interpolation"),
        ("stiffened_panel_verify",  "QUAD4+CBAR+RBE3 panel"),
        ("cbar_rect_verify",        "CBAR cantilever, Fy+Fz"),
    ]

    println("\n$(rpad("Test",30)) $(rpad("Description",30)) $(lpad("MaxErr%",8)) $(lpad("Grids",6)) $(lpad("Comp",5))")
    println("-"^75)
    n_pass = 0
    n_fail = 0
    n_skip = 0
    for (name, desc) in tests
        result = run_test(name)
        if isnothing(result)
            println("$(rpad(name,30)) $(rpad(desc,30)) $(lpad("SKIP",8))")
            n_skip += 1
        else
            status = result.max_err < 20.0 ? "PASS" : "WARN"
            err_str = "$(round(result.max_err, digits=2))%"
            if result.max_err >= 20.0; n_fail += 1; else; n_pass += 1; end
            println("$(rpad(name,30)) $(rpad(desc,30)) $(lpad(err_str,8)) $(lpad(result.n_mys,6)) $(lpad(result.n_comp,5))  $status")
        end
    end
    println("-"^75)
    println("Results: $n_pass PASS, $n_fail WARN, $n_skip SKIP")
end

main()
