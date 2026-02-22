using JSON

function compare_disp(ju, mys, label)
    println("\n=== $label: DISPLACEMENTS ===")
    # Parse MYSTRAN F06 displacements
    mys_disp = Dict{Int,Dict{String,Float64}}()
    lines = split(mys, "\n")
    in_disp = false
    for line in lines
        if occursin("D I S P L A C E M E N T S", line)
            in_disp = true; continue
        end
        if in_disp && occursin("SYS", line); continue; end
        if in_disp && occursin("-----", line); in_disp = false; continue; end
        if in_disp && occursin("MAX*", line); in_disp = false; continue; end
        if in_disp
            parts = split(strip(line))
            if length(parts) >= 7
                gid = tryparse(Int, parts[1])
                if gid !== nothing && tryparse(Int, parts[2]) !== nothing
                    mys_disp[gid] = Dict(
                        "t1" => parse(Float64, parts[3]),
                        "t2" => parse(Float64, parts[4]),
                        "t3" => parse(Float64, parts[5]),
                        "r1" => parse(Float64, parts[6]),
                        "r2" => parse(Float64, parts[7]),
                        "r3" => parse(Float64, parts[8])
                    )
                end
            end
        end
    end

    ju_disp = Dict{Int,Any}()
    for d in ju["displacements"]; ju_disp[d["grid_id"]] = d; end

    max_err = Dict("t1"=>0.0, "t2"=>0.0, "t3"=>0.0, "r1"=>0.0, "r2"=>0.0, "r3"=>0.0)
    for (gid, md) in mys_disp
        if !haskey(ju_disp, gid); continue; end
        jd = ju_disp[gid]
        for k in ["t1","t2","t3","r1","r2","r3"]
            ref = abs(md[k])
            if ref > 1e-10
                err = abs(jd[k] - md[k]) / ref * 100
                max_err[k] = max(max_err[k], err)
                if err > 1.0
                    println("  Grid $gid $k: Julia=$(round(jd[k],sigdigits=6)) MYSTRAN=$(round(md[k],sigdigits=6)) err=$(round(err,digits=2))%")
                end
            end
        end
    end
    println("  Max errors: ", join(["$k=$(round(v,digits=2))%" for (k,v) in sort(collect(max_err))], "  "))
end

function compare_spc(ju, mys, label)
    println("\n=== $label: SPC FORCES ===")
    lines = split(mys, "\n")
    mys_spc = Dict{Int,Dict{String,Float64}}()
    in_spc = false
    for line in lines
        if occursin("S P C   F O R C E S", line)
            in_spc = true; continue
        end
        if in_spc && occursin("SYS", line); continue; end
        if in_spc && (occursin("-----", line) || occursin("SPC FORCE", line) || occursin("AUTOSPC", line))
            in_spc = false; continue
        end
        if in_spc
            parts = split(strip(line))
            if length(parts) >= 7
                gid = tryparse(Int, parts[1])
                if gid !== nothing && tryparse(Int, parts[2]) !== nothing
                    mys_spc[gid] = Dict(
                        "t1" => parse(Float64, parts[3]),
                        "t2" => parse(Float64, parts[4]),
                        "t3" => parse(Float64, parts[5]),
                        "r1" => parse(Float64, parts[6]),
                        "r2" => parse(Float64, parts[7]),
                        "r3" => parse(Float64, parts[8])
                    )
                end
            end
        end
    end

    ju_spc = Dict{Int,Any}()
    for d in ju["spc_forces"]; ju_spc[d["grid_id"]] = d; end

    for (gid, ms) in sort(collect(mys_spc), by=first)
        if !haskey(ju_spc, gid)
            println("  Grid $gid: MISSING in Julia")
            continue
        end
        js = ju_spc[gid]
        for k in ["t1","t2","t3"]
            ref = abs(ms[k])
            if ref > 1e-6
                err = abs(js[k] - ms[k]) / ref * 100
                if err > 0.1
                    println("  Grid $gid $k: Julia=$(round(js[k],sigdigits=6)) MYSTRAN=$(round(ms[k],sigdigits=6)) err=$(round(err,digits=2))%")
                end
            end
        end
    end
    # Total SPC
    ju_total = [sum(s[k] for s in values(ju_spc)) for k in ["t1","t2","t3"]]
    mys_total = [sum(s[k] for s in values(mys_spc)) for k in ["t1","t2","t3"]]
    println("  SPC Total: Julia=[$(round.(ju_total,sigdigits=6))]  MYSTRAN=[$(round.(mys_total,sigdigits=6))]")
end

function compare_quad_forces(ju, mys, label)
    println("\n=== $label: QUAD4 FORCES ===")
    lines = split(mys, "\n")
    mys_qf = Dict{Int,Dict{String,Float64}}()
    in_qf = false
    for line in lines
        if occursin("E L E M E N T   E N G I N E E R I N G   F O R C E S", line) && occursin("Q U A D 4", lines[min(end, findfirst(==(line), lines)+1)])
            in_qf = true; continue
        end
        if occursin("Q U A D 4", line) && in_qf; continue; end
        if in_qf && (occursin("Nxx", line) || occursin("Element Location", line) || occursin("Shear", line)); continue; end
        if in_qf && occursin("-----", line); in_qf = false; continue; end
        if in_qf && occursin("MAX*", line); in_qf = false; continue; end
        if in_qf
            parts = split(strip(line))
            if length(parts) >= 9
                eid = tryparse(Int, parts[1])
                if eid !== nothing
                    mys_qf[eid] = Dict(
                        "fx" => parse(Float64, parts[2]),
                        "fy" => parse(Float64, parts[3]),
                        "fxy" => parse(Float64, parts[4]),
                        "mx" => parse(Float64, parts[5]),
                        "my" => parse(Float64, parts[6]),
                        "mxy" => parse(Float64, parts[7]),
                        "qx" => parse(Float64, parts[8]),
                        "qy" => parse(Float64, parts[9])
                    )
                end
            end
        end
    end

    ju_qf = Dict{Int,Any}()
    for f in ju["forces"]["quad4"]; ju_qf[f["eid"]] = f; end

    for (eid, mf) in sort(collect(mys_qf), by=first)
        if !haskey(ju_qf, eid)
            println("  Elem $eid: MISSING in Julia")
            continue
        end
        jf = ju_qf[eid]
        for k in ["fx","fy","fxy"]
            ref = abs(mf[k])
            if ref > 1e-6
                err = abs(jf[k] - mf[k]) / ref * 100
                println("  Elem $eid $k: Julia=$(round(jf[k],sigdigits=6)) MYSTRAN=$(round(mf[k],sigdigits=6)) err=$(round(err,digits=2))%")
            else
                if abs(jf[k]) > 1e-6
                    println("  Elem $eid $k: Julia=$(round(jf[k],sigdigits=6)) MYSTRAN=~0  [MISMATCH]")
                end
            end
        end
    end
end

function compare_quad_stress(ju, mys, label)
    println("\n=== $label: QUAD4 STRESSES ===")
    lines = split(mys, "\n")
    mys_qs = Dict{Int,Dict{String,Float64}}()
    in_qs = false
    for (i, line) in enumerate(lines)
        if occursin("E L E M E N T   S T R E S S E S", line) && i+1 <= length(lines) && occursin("Q U A D 4", lines[i+1])
            in_qs = true; continue
        end
        if in_qs && (occursin("Q U A D 4", line) || occursin("Elem", line) || occursin("Distance", line) || occursin("max through", line))
            continue
        end
        if in_qs && occursin("-----", line); in_qs = false; continue; end
        if in_qs && occursin("MAX*", line); in_qs = false; continue; end
        if in_qs
            parts = split(strip(line))
            if length(parts) >= 10
                eid = tryparse(Int, parts[1])
                if eid !== nothing && parts[2] == "CENTER"
                    # First line (bottom face)
                    mys_qs[eid] = Dict(
                        "normal_x" => parse(Float64, parts[4]),
                        "normal_y" => parse(Float64, parts[5]),
                        "shear_xy" => parse(Float64, parts[6]),
                        "von_mises" => parse(Float64, parts[10])
                    )
                end
            end
        end
    end

    ju_qs = Dict{Int,Any}()
    for s in ju["stresses"]["quad4"]; ju_qs[s["eid"]] = s; end

    for (eid, ms) in sort(collect(mys_qs), by=first)
        if !haskey(ju_qs, eid)
            println("  Elem $eid: MISSING in Julia")
            continue
        end
        js = ju_qs[eid]
        # Julia has z1/z2 layers
        jl = js["z1"]
        for k in ["normal_x","normal_y","shear_xy","von_mises"]
            ref = abs(ms[k])
            if ref > 1e-6
                err = abs(jl[k] - ms[k]) / ref * 100
                println("  Elem $eid $k: Julia=$(round(jl[k],sigdigits=6)) MYSTRAN=$(round(ms[k],sigdigits=6)) err=$(round(err,digits=2))%")
            end
        end
    end
end

function compare_crod(ju, mys, label)
    println("\n=== $label: CROD FORCES & STRESSES ===")
    lines = split(mys, "\n")

    # Parse CROD forces
    mys_rf = Dict{Int,Float64}()
    in_rf = false
    for line in lines
        if occursin("F O R   E L E M E N T   T Y P E   R O D", line) && occursin("F O R C E S", lines[max(1,findfirst(==(line), lines)-1)])
            in_rf = true; continue
        end
        if in_rf && occursin("Element", line); continue; end
        if in_rf && occursin("-----", line); in_rf = false; continue; end
        if in_rf && occursin("MAX*", line); in_rf = false; continue; end
        if in_rf
            parts = split(strip(line))
            i = 1
            while i + 1 <= length(parts)
                eid = tryparse(Int, parts[i])
                if eid !== nothing
                    mys_rf[eid] = parse(Float64, parts[i+1])
                    i += 3  # skip eid, axial, torque
                else
                    i += 1
                end
            end
        end
    end

    # Parse CROD stresses
    mys_rs = Dict{Int,Float64}()
    in_rs = false
    for line in lines
        if occursin("F O R   E L E M E N T   T Y P E   R O D", line) && occursin("S T R E S S E S", lines[max(1,findfirst(==(line), lines)-1)])
            in_rs = true; continue
        end
        if in_rs && occursin("Element", line); continue; end
        if in_rs && occursin("-----", line); in_rs = false; continue; end
        if in_rs && occursin("MAX*", line); in_rs = false; continue; end
        if in_rs
            parts = split(strip(line))
            i = 1
            while i + 1 <= length(parts)
                eid = tryparse(Int, parts[i])
                if eid !== nothing
                    mys_rs[eid] = parse(Float64, parts[i+1])
                    i += 4  # skip eid, axial_stress, margin, torsional_stress
                else
                    i += 1
                end
            end
        end
    end

    ju_rf = Dict{Int,Any}()
    for f in ju["forces"]["crod"]; ju_rf[f["eid"]] = f; end

    ju_rs = Dict{Int,Any}()
    for s in ju["stresses"]["crod"]; ju_rs[s["eid"]] = s; end

    for (eid, mf) in sort(collect(mys_rf), by=first)
        jf = get(ju_rf, eid, nothing)
        if jf === nothing
            println("  Elem $eid: MISSING in Julia forces")
            continue
        end
        ref = abs(mf)
        err_f = ref > 1e-6 ? abs(jf["axial"] - mf) / ref * 100 : 0.0
        ms = get(mys_rs, eid, 0.0)
        js = get(ju_rs, eid, nothing)
        err_s = 0.0
        if js !== nothing && abs(ms) > 1e-6
            err_s = abs(js["axial"] - ms) / abs(ms) * 100
        end
        println("  Elem $eid: Force Julia=$(round(jf["axial"],sigdigits=6)) MYSTRAN=$(round(mf,sigdigits=6)) err=$(round(err_f,digits=2))%  |  Stress Julia=$(js !== nothing ? round(js["axial"],sigdigits=6) : "N/A") MYSTRAN=$(round(ms,sigdigits=6)) err=$(round(err_s,digits=2))%")
    end
end

# ============================================================
# Run all comparisons
# ============================================================
base_ju = joinpath(@__DIR__, "..", "..", "output") * "/"
base_mys = joinpath(@__DIR__, "..", "..", "..", "MYSTRAN", "MSC") * "/"

tests = [
    ("quad4_membrane_verify", "QUAD4 MEMBRANE"),
    ("crod_axial_verify", "CROD AXIAL"),
    ("rbe2_verify", "RBE2 RIGID LINK"),
    ("rbe3_verify", "RBE3 LOAD DIST"),
]

for (name, label) in tests
    println("\n" * "="^70)
    println(" TEST: $label ($name)")
    println("="^70)

    ju = JSON.parsefile(base_ju * name * ".JU.JSON")
    mys = read(base_mys * name * ".F06", String)

    compare_disp(ju, mys, label)
    compare_spc(ju, mys, label)

    if haskey(ju["forces"], "quad4") && length(ju["forces"]["quad4"]) > 0
        compare_quad_forces(ju, mys, label)
        compare_quad_stress(ju, mys, label)
    end
    if haskey(ju["forces"], "crod") && length(ju["forces"]["crod"]) > 0
        compare_crod(ju, mys, label)
    end
end
