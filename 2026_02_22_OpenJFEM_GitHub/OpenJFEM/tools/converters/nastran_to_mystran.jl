"""
    nastran_to_mystran.jl

Converts a NASTRAN-formatted .bdf file to a MYSTRAN-formatted .bdf file,
runs MYSTRAN, and collects results into a dedicated output folder.

Usage:
    julia nastran_to_mystran.jl <input.bdf> [mystran_exe_path]

Output:
    <input>.MYSTRAN.bdf          - converted MYSTRAN input file (next to original)
    <input_name>_MYST_RESULTS/   - folder with all MYSTRAN output files

The MYSTRAN executable path can be supplied as the second argument.
If omitted, the script looks for MYSTRAN_EXE environment variable,
then falls back to the default project location.
"""

using Printf

# ─────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────

const DEFAULT_MYSTRAN_EXE = joinpath(@__DIR__, "..", "..", "..",
    "MYSTRAN", "mystran-17.0.0-windows-x86_64.exe")

# MYSTRAN output extensions to collect into the results folder
const MYSTRAN_OUTPUT_EXTS = [
    ".F06", ".ERR", ".OP2", ".NEU", ".BUG", ".ANS",
    ".L1F", ".L1I", ".L1K", ".L1N", ".L1O", ".L1P",
    ".L1Q", ".L1S", ".L1T", ".L1U", ".L1V", ".L1W", ".L1X",
]

# SOL number mapping: NASTRAN → MYSTRAN
const SOL_MAP = Dict(
    "101" => "1",   # Linear statics
    "103" => "3",   # Normal modes
    "105" => "5",   # Buckling
    "111" => "11",  # Modal frequency response
    "112" => "12",  # Modal transient response
)

# ─────────────────────────────────────────────────
# Case Control output request conversion
# ─────────────────────────────────────────────────

"""
Maps NASTRAN case control output keywords to MYSTRAN equivalents.
Each entry: nastran_keyword => (mystran_keyword, mystran_options)
"""
const OUTPUT_REQUEST_MAP = Dict(
    "DISPLACEMENT" => ("DISP",     "(PRINT,PLOT,PUNCH)"),
    "SPCFORCES"    => ("SPCFORCE", "(PRINT,PLOT,PUNCH)"),
    "FORCE"        => ("ELFORCE",  "(PRINT,PLOT,PUNCH)"),
    "STRESS"       => ("STRESS",   "(PRINT,PLOT,PUNCH)"),
    "STRAIN"       => ("STRAIN",   "(PRINT,PLOT,PUNCH)"),
    "ESE"          => ("ESE",      "(PRINT,PLOT,PUNCH)"),
    "MPCFORCES"    => ("MPCFORCE", "(PRINT,PLOT,PUNCH)"),
    "OLOAD"        => ("OLOAD",    "(PRINT,PLOT,PUNCH)"),
    "GPFORCE"      => ("GPFORCE",  "(PRINT,PLOT)"),
)

# Additional MYSTRAN output requests not commonly in NASTRAN files
const MYSTRAN_EXTRA_OUTPUTS = [
    "OLOAD(PRINT,PLOT,PUNCH) = ALL",
    "GPFORCE(PRINT,PLOT) = ALL",
    "MPCFORCE(PRINT,PLOT,PUNCH) = ALL",
]

# MYSTRAN PARAM/DEBUG cards to inject into Bulk Data (if missing)
const MYSTRAN_PARAMS = [
    "PARAM   SOLLIB   IntMKL",
    "PARAM   GRDPNT   0",
    "PARAM   POST     -1",
    "PARAM   AUTOSPC  YES",
    "PARAM   BAILOUT  -1",
]

const MYSTRAN_DEBUG = [
    "DEBUG   192     2                                                       GPFO summary",
    "DEBUG   200     1                                                       ANS",
]

# ─────────────────────────────────────────────────
# Parsing helpers
# ─────────────────────────────────────────────────

"""Return the trimmed uppercase version of a line for keyword matching."""
function norm(line::AbstractString)
    return uppercase(strip(line))
end

"""Check if a line is a comment or blank."""
function is_comment_or_blank(line::AbstractString)
    s = lstrip(line)
    return isempty(s) || startswith(s, "\$")
end

"""
Extract the keyword from a case control line like:
  "  DISPLACEMENT(SORT1,REAL) = ALL"  →  "DISPLACEMENT"
  "  FORCE = ALL"                     →  "FORCE"
  "  SUBTITLE = TIP LOAD"            →  "SUBTITLE"
"""
function extract_cc_keyword(line::AbstractString)
    s = strip(line)
    if startswith(s, "\$") || isempty(s)
        return nothing
    end
    # Match: KEYWORD, KEYWORD(, KEYWORD =, KEYWORD=
    m = match(r"^([A-Za-z]+)", s)
    return m === nothing ? nothing : uppercase(m.captures[1])
end

# ─────────────────────────────────────────────────
# SPC1 THRU reformatter
# ─────────────────────────────────────────────────

"""
    reformat_spc1_thru(line) -> Vector{String}

Splits SPC1 cards that mix individual grid IDs with THRU ranges into
separate cards. MYSTRAN cannot parse mixed-format SPC1 cards like:
    SPC1  2  123  G1  G2  G3  THRU  G_end
This is split into:
    SPC1  2  123  G1  G2
    SPC1  2  123  G3  THRU  G_end
Simple THRU-only cards (SPC1 SID C G_start THRU G_end) pass through unchanged.
"""
function reformat_spc1_thru(line::AbstractString)
    f = bdf_fields(line)
    if isempty(f) || uppercase(strip(f[1])) != "SPC1"
        return [line]
    end

    # Find THRU keyword position (search from field 4 onward)
    thru_pos = 0
    for k in 4:length(f)
        if uppercase(strip(f[k])) == "THRU"
            thru_pos = k
            break
        end
    end

    # No THRU — pass through unchanged
    if thru_pos == 0
        return [line]
    end

    # Simple THRU format: field 5 = THRU (e.g. SPC1 SID C G_start THRU G_end)
    if thru_pos == 5
        return [line]
    end

    # Mixed format: individual grids in fields 4..(thru_pos-2), then
    # field (thru_pos-1) is the THRU range start, field (thru_pos+1) is range end
    sid = strip(f[2])
    dof = strip(f[3])

    # Individual grids: fields 4 to (thru_pos - 2)
    individual_grids = String[]
    for k in 4:(thru_pos - 2)
        g = strip(f[k])
        !isempty(g) && push!(individual_grids, g)
    end

    # THRU range: field (thru_pos-1) THRU field (thru_pos+1)
    thru_start = strip(f[thru_pos - 1])
    thru_end   = thru_pos + 1 <= length(f) ? strip(f[thru_pos + 1]) : ""

    result = String[]

    # Card 1: individual grids (if any)
    if !isempty(individual_grids)
        card = rpad("SPC1", 8) * lpad(sid, 8) * lpad(dof, 8)
        for g in individual_grids
            card *= lpad(g, 8)
        end
        push!(result, card)
    end

    # Card 2: THRU range
    if !isempty(thru_start) && !isempty(thru_end)
        card = rpad("SPC1", 8) * lpad(sid, 8) * lpad(dof, 8) *
               lpad(thru_start, 8) * rpad("    THRU", 8) * lpad(thru_end, 8)
        push!(result, card)
    end

    return result
end

"""
    filter_spc1_mpc_grids(line, mpc_grids) -> Vector{String}

Removes MPC/RBE-dependent grid IDs from an SPC1 card.
MYSTRAN cannot have the same DOF in both SPC and MPC displacement sets.
Returns empty vector if all grids are filtered out.
Handles both individual-grid and THRU-range formats.
"""
function filter_spc1_mpc_grids(line::AbstractString, mpc_grids::Set{Int})
    isempty(mpc_grids) && return [line]

    f = bdf_fields(line)
    if isempty(f) || uppercase(strip(f[1])) != "SPC1"
        return [line]
    end

    sid = strip(f[2])
    dof = strip(f[3])

    # Check for THRU format (field 5 = THRU)
    if length(f) >= 6 && uppercase(strip(get(f, 5, ""))) == "THRU"
        start_id = parse_bdf_int(f[4])
        end_id   = parse_bdf_int(f[6])
        # For THRU ranges, expand only if a few MPC grids are in the range
        # and re-emit as individual IDs minus the MPC ones
        conflict_grids = [g for g in mpc_grids if start_id <= g <= end_id]
        if isempty(conflict_grids)
            return [line]  # No conflicts, pass through
        end
        # Expand the range and filter out MPC grids
        remaining = [g for g in start_id:end_id if !(g in mpc_grids)]
        if isempty(remaining)
            return String[]  # All grids filtered out
        end
        # Re-emit as individual grid SPC1 cards (max 5 grids per card)
        result = String[]
        for chunk in Iterators.partition(remaining, 5)
            card = rpad("SPC1", 8) * lpad(sid, 8) * lpad(dof, 8)
            for g in chunk
                card *= lpad(string(g), 8)
            end
            push!(result, card)
        end
        return result
    end

    # Individual grid format: filter out MPC-dependent grids
    remaining_grids = String[]
    for k in 4:length(f)
        gs = strip(f[k])
        isempty(gs) && continue
        gid = tryparse(Int, gs)
        if gid === nothing || !(gid in mpc_grids)
            push!(remaining_grids, gs)
        end
    end

    if isempty(remaining_grids)
        return String[]  # All grids filtered out
    end

    # Re-emit with remaining grids
    card = rpad("SPC1", 8) * lpad(sid, 8) * lpad(dof, 8)
    for g in remaining_grids
        card *= lpad(g, 8)
    end
    return [card]
end

# ─────────────────────────────────────────────────
# Core conversion
# ─────────────────────────────────────────────────

"""
    convert_nastran_to_mystran(input_path::String) -> String

Reads a NASTRAN .bdf file at `input_path`, converts it to MYSTRAN format,
and returns the converted content as a string.
"""
function convert_nastran_to_mystran(input_path::String)
    lines = readlines(input_path; keep=false)

    # ── Phase 1: identify section boundaries ──
    cend_idx = 0
    begin_bulk_idx = 0
    enddata_idx = 0
    has_id_card = false

    for (i, line) in enumerate(lines)
        n = norm(line)
        if n == "CEND"
            cend_idx = i
        elseif startswith(n, "BEGIN") && occursin("BULK", n)
            begin_bulk_idx = i
        elseif n == "ENDDATA"
            enddata_idx = i
        elseif startswith(n, "ID ")
            has_id_card = true
        end
    end

    if cend_idx == 0 || begin_bulk_idx == 0
        error("Cannot find CEND or BEGIN BULK in $input_path — is this a valid NASTRAN .bdf?")
    end
    if enddata_idx == 0
        enddata_idx = length(lines) + 1
    end

    exec_lines = lines[1:cend_idx]           # includes CEND
    cc_lines   = lines[cend_idx+1:begin_bulk_idx-1]
    bulk_lines = lines[begin_bulk_idx+1 : (enddata_idx > 0 ? enddata_idx-1 : length(lines))]

    # ── Phase 2: convert Executive Control ──
    out_exec = String[]
    model_name = basename(input_path)
    model_name = replace(model_name, r"\.[^.]+$" => "")  # strip extension

    for line in exec_lines
        n = norm(line)
        if startswith(n, "ID ")
            # Keep existing ID card
            push!(out_exec, line)
            has_id_card = true
        elseif startswith(n, "SOL")
            # Insert ID card before SOL if not present
            if !has_id_card
                push!(out_exec, "ID $(uppercase(model_name)),STATIC")
                has_id_card = true
            end
            # Convert SOL number
            m = match(r"^(\s*SOL\s+)(\d+)(.*)", line, )
            if m !== nothing
                sol_num = strip(m.captures[2])
                mystran_sol = get(SOL_MAP, sol_num, sol_num)
                push!(out_exec, "$(m.captures[1])$(mystran_sol)$(m.captures[3])")
            else
                push!(out_exec, line)
            end
        else
            push!(out_exec, line)
        end
    end

    # ── Phase 3: convert Case Control ──
    out_cc = String[]
    seen_keywords = Set{String}()   # track which output requests are present
    has_echo = false
    has_spc_global = false
    global_spc_value = ""

    # First pass: detect existing features
    for line in cc_lines
        kw = extract_cc_keyword(line)
        if kw !== nothing
            push!(seen_keywords, kw)
            if kw == "ECHO"
                has_echo = true
            end
        end
    end

    # Add ECHO = UNSORT if not present (before other output requests)
    if !has_echo
        push!(out_cc, "ECHO = UNSORT")
    end

    for line in cc_lines
        n = norm(line)

        # Skip if comment/blank — pass through
        if is_comment_or_blank(line)
            push!(out_cc, line)
            continue
        end

        kw = extract_cc_keyword(line)

        # Convert SUBTITLE → LABEL
        if kw == "SUBTITLE"
            new_line = replace(line, r"SUBTITLE"i => "LABEL")
            push!(out_cc, new_line)
            continue
        end

        # Convert output request keywords
        if kw !== nothing && haskey(OUTPUT_REQUEST_MAP, kw)
            mystran_kw, mystran_opts = OUTPUT_REQUEST_MAP[kw]
            # Extract RHS (= ALL, = 1, etc.)
            m = match(r"=\s*(.+)$", line)
            rhs = m !== nothing ? strip(m.captures[1]) : "ALL"
            # Preserve indentation
            indent = match(r"^(\s*)", line).captures[1]
            push!(out_cc, "$(indent)$(mystran_kw)$(mystran_opts) = $(rhs)")
            continue
        end

        # Pass through everything else (TITLE, SPC, LOAD, SUBCASE, SET, etc.)
        push!(out_cc, line)
    end

    # Add missing MYSTRAN output requests at the global level (before first SUBCASE)
    # Find insertion point: before first SUBCASE or at end of case control
    insert_idx = length(out_cc) + 1
    for (i, line) in enumerate(out_cc)
        n = norm(line)
        if startswith(n, "SUBCASE")
            insert_idx = i
            break
        end
    end

    extras_to_add = String[]
    for extra_line in MYSTRAN_EXTRA_OUTPUTS
        extra_kw = extract_cc_keyword(extra_line)
        if extra_kw !== nothing && !(extra_kw in seen_keywords)
            push!(extras_to_add, extra_line)
        end
    end

    if !isempty(extras_to_add)
        for (j, extra) in enumerate(extras_to_add)
            insert!(out_cc, insert_idx + j - 1, extra)
        end
    end

    # ── Phase 4: convert Bulk Data ──

    # Pre-scan: resolve SPCADD cards (MYSTRAN does not support SPCADD)
    # SPCADD target_sid  sub_sid1  sub_sid2 ...
    # → renumber SPC/SPC1 cards from sub_sid to target_sid and drop SPCADD
    spcadd_map = Dict{Int, Int}()  # sub_sid => target_sid
    for line in bulk_lines
        n = norm(line)
        if startswith(n, "SPCADD")
            f = bdf_fields(line)
            if length(f) >= 3
                target = parse_bdf_int(f[2])
                for k in 3:length(f)
                    sub = parse_bdf_int(f[k])
                    sub > 0 && target > 0 && (spcadd_map[sub] = target)
                end
            end
        end
    end
    if !isempty(spcadd_map)
        println("  Resolving SPCADD: $(join(["$sub→$tgt" for (sub,tgt) in spcadd_map], ", "))")
    end

    # Pre-scan: collect MPC/RBE-dependent grids to filter from SPC1 cards
    # MYSTRAN cannot have the same DOF in both SPC and MPC/RBE sets
    mpc_dependent_grids = Set{Int}()
    for line in bulk_lines
        n = norm(line)
        f = nothing
        if startswith(n, "RBE3")
            # RBE3: dependent (reference) grid is in field 4
            f = bdf_fields(line)
            if length(f) >= 4
                gid = parse_bdf_int(f[4])
                gid > 0 && push!(mpc_dependent_grids, gid)
            end
        elseif startswith(n, "RBE2")
            # RBE2: dependent grids are in fields 5+ (field 3=independent, field 4=components)
            f = bdf_fields(line)
            for k in 5:length(f)
                gid = parse_bdf_int(f[k])
                gid > 0 && push!(mpc_dependent_grids, gid)
            end
        elseif startswith(n, "MPC") && !startswith(n, "MPCADD")
            # MPC: dependent grid is in field 3 (f[1]=MPC, f[2]=SID, f[3]=G1_dep)
            f = bdf_fields(line)
            if length(f) >= 3
                gid = parse_bdf_int(f[3])
                gid > 0 && push!(mpc_dependent_grids, gid)
            end
        end
    end
    if !isempty(mpc_dependent_grids)
        println("  Found $(length(mpc_dependent_grids)) MPC/RBE-dependent grids to filter from SPC1")
    end

    # Pre-scan: detect K6ROT and collect data for drilling spring injection
    k6rot_value = 0.0
    k6rot_grids    = Dict{Int, NTuple{3,Float64}}()   # nid => (x,y,z)
    k6rot_mat1     = Dict{Int, Tuple{Float64,Float64}}()  # mid => (E, ν)
    k6rot_mat8_G   = Dict{Int, Float64}()                 # mid => G12 (in-plane shear for MAT8)
    k6rot_pshell   = Dict{Int, Tuple{Int,Float64}}()   # pid => (mid, thickness)
    k6rot_pcomp    = Dict{Int, Tuple{Int,Float64}}()   # pid => (mid_ply1, total_thickness)
    k6rot_quads    = Vector{Tuple{Int,Int,Int,Int,Int}}()  # (pid, g1, g2, g3, g4)
    k6rot_trias    = Vector{Tuple{Int,Int,Int,Int}}()      # (pid, g1, g2, g3)
    k6rot_max_eid  = 0
    k6rot_max_pid  = 0

    for line in bulk_lines
        n = norm(line)
        f = bdf_fields(line)
        isempty(f) && continue
        card = uppercase(f[1])

        if startswith(n, "PARAM") && occursin("K6ROT", n)
            # Extract K6ROT value
            if length(f) >= 3
                k6rot_value = parse_nastran_float(f[3])
            end
        end

        # Collect data needed for K6ROT spring computation
        if card == "GRID" && length(f) >= 6
            nid = parse_bdf_int(f[2])
            x = parse_nastran_float(length(f) >= 4 ? f[4] : "0")
            y = parse_nastran_float(length(f) >= 5 ? f[5] : "0")
            z = parse_nastran_float(length(f) >= 6 ? f[6] : "0")
            nid > 0 && (k6rot_grids[nid] = (x, y, z))
        elseif card == "MAT1" && length(f) >= 5
            mid = parse_bdf_int(f[2])
            E = parse_nastran_float(f[3])
            nu = parse_nastran_float(f[5])
            mid > 0 && E > 0 && (k6rot_mat1[mid] = (E, nu))
        elseif card == "MAT8" && length(f) >= 6
            mid = parse_bdf_int(f[2])
            E1 = parse_nastran_float(f[3])
            G12 = length(f) >= 6 ? parse_nastran_float(f[6]) : E1 / 6.0
            mid > 0 && (k6rot_mat8_G[mid] = G12 > 0 ? G12 : E1 / 6.0)
        elseif card == "PSHELL" && length(f) >= 4
            pid = parse_bdf_int(f[2]); mid = parse_bdf_int(f[3])
            t = parse_nastran_float(f[4])
            pid > 0 && (k6rot_pshell[pid] = (mid, t))
            pid > k6rot_max_pid && (k6rot_max_pid = pid)
        elseif card == "PCOMP" && length(f) >= 2
            pid = parse_bdf_int(f[2])
            pid > k6rot_max_pid && (k6rot_max_pid = pid)
            # Use Z0 field (field 3) to estimate total thickness: for SYM, total_t ≈ 2*|Z0|
            z0 = length(f) >= 3 ? abs(parse_nastran_float(f[3])) : 0.0
            t_est = z0 > 0 ? 2.0 * z0 : 0.0
            is_sym = any(uppercase(strip(fi)) == "SYM" for fi in f)
            # Try to get first ply MID from continuation (fields 10+)
            mid_ply = length(f) >= 10 ? parse_bdf_int(f[10]) : 0
            # Also get ply thickness from continuation for sum
            if t_est <= 0 && length(f) >= 11
                ply_t = parse_nastran_float(f[11])
                n_plies_here = div(length(f) - 9, 4)  # 4 fields per ply
                t_est = ply_t * n_plies_here * (is_sym ? 2 : 1)
            end
            pid > 0 && (k6rot_pcomp[pid] = (mid_ply, t_est))
        elseif card == "CQUAD4" && length(f) >= 7
            eid = parse_bdf_int(f[2]); pid = parse_bdf_int(f[3])
            g1 = parse_bdf_int(f[4]); g2 = parse_bdf_int(f[5])
            g3 = parse_bdf_int(f[6]); g4 = parse_bdf_int(f[7])
            push!(k6rot_quads, (pid, g1, g2, g3, g4))
            eid > k6rot_max_eid && (k6rot_max_eid = eid)
        elseif card == "CTRIA3" && length(f) >= 6
            eid = parse_bdf_int(f[2]); pid = parse_bdf_int(f[3])
            g1 = parse_bdf_int(f[4]); g2 = parse_bdf_int(f[5]); g3 = parse_bdf_int(f[6])
            push!(k6rot_trias, (pid, g1, g2, g3))
            eid > k6rot_max_eid && (k6rot_max_eid = eid)
        end
        # Track max element/property IDs for any card
        if card in ("CBAR","CROD","CELAS1","CONROD","CELAS2") && length(f) >= 2
            eid = parse_bdf_int(f[2])
            eid > k6rot_max_eid && (k6rot_max_eid = eid)
        end
        if card in ("PBAR","PBARL","PROD","PELAS") && length(f) >= 2
            pid = parse_bdf_int(f[2])
            pid > k6rot_max_pid && (k6rot_max_pid = pid)
        end
    end

    # Generate K6ROT drilling springs if K6ROT parameter is present
    k6rot_bulk_lines = String[]
    if k6rot_value > 0 && (!isempty(k6rot_quads) || !isempty(k6rot_trias))
        # Compute per-node spring stiffness using Hughes-Brezzi formula:
        #   alpha_drill = (k6rot / 1e5) * G * h   (per unit area)
        #   K_node = alpha_drill * A_elem / n_nodes_per_elem  (accumulated)

        # Helper: get G (shear modulus), h (thickness) for a shell property
        function _get_shell_Gh(pid)
            # Helper to get G from a material ID (MAT1 or MAT8)
            function _get_G(mid)
                if haskey(k6rot_mat1, mid)
                    E, nu = k6rot_mat1[mid]
                    return E / (2.0 * (1.0 + nu))
                elseif haskey(k6rot_mat8_G, mid)
                    return k6rot_mat8_G[mid]
                end
                return 0.0
            end

            if haskey(k6rot_pshell, pid)
                mid, h = k6rot_pshell[pid]
                G = _get_G(mid)
                G > 0 && return G, h
            end
            # PCOMP: use total thickness from Z0 estimate and ply material
            if haskey(k6rot_pcomp, pid)
                mid, t_total = k6rot_pcomp[pid]
                G = mid > 0 ? _get_G(mid) : 0.0
                if G <= 0
                    # Fallback: use any available material
                    if !isempty(k6rot_mat1)
                        E, nu = first(values(k6rot_mat1))
                        G = E / (2.0 * (1.0 + nu))
                    elseif !isempty(k6rot_mat8_G)
                        G = first(values(k6rot_mat8_G))
                    end
                end
                G > 0 && return G, max(t_total, 0.1)
            end
            return 0.0, 0.0
        end

        # Helper: triangle area from 3D coords
        function _tri_area(p1, p2, p3)
            e1 = (p2[1]-p1[1], p2[2]-p1[2], p2[3]-p1[3])
            e2 = (p3[1]-p1[1], p3[2]-p1[2], p3[3]-p1[3])
            cx = e1[2]*e2[3] - e1[3]*e2[2]
            cy = e1[3]*e2[1] - e1[1]*e2[3]
            cz = e1[1]*e2[2] - e1[2]*e2[1]
            return 0.5 * sqrt(cx^2 + cy^2 + cz^2)
        end

        # Accumulate per-node spring stiffness
        node_K = Dict{Int, Float64}()  # nid => accumulated drilling stiffness

        for (pid, g1, g2, g3, g4) in k6rot_quads
            G, h = _get_shell_Gh(pid)
            (G <= 0 || h <= 0) && continue
            alpha = (k6rot_value / 1e5) * G * h
            p1 = get(k6rot_grids, g1, nothing); p2 = get(k6rot_grids, g2, nothing)
            p3 = get(k6rot_grids, g3, nothing); p4 = get(k6rot_grids, g4, nothing)
            (p1 === nothing || p2 === nothing || p3 === nothing || p4 === nothing) && continue
            A = _tri_area(p1, p2, p3) + _tri_area(p1, p3, p4)
            K_per_node = alpha * A / 4.0
            for g in (g1, g2, g3, g4)
                node_K[g] = get(node_K, g, 0.0) + K_per_node
            end
        end

        for (pid, g1, g2, g3) in k6rot_trias
            G, h = _get_shell_Gh(pid)
            (G <= 0 || h <= 0) && continue
            alpha = (k6rot_value / 1e5) * G * h
            p1 = get(k6rot_grids, g1, nothing); p2 = get(k6rot_grids, g2, nothing)
            p3 = get(k6rot_grids, g3, nothing)
            (p1 === nothing || p2 === nothing || p3 === nothing) && continue
            A = _tri_area(p1, p2, p3)
            K_per_node = alpha * A / 3.0
            for g in (g1, g2, g3)
                node_K[g] = get(node_K, g, 0.0) + K_per_node
            end
        end

        # Filter out MPC/RBE-dependent grids (can't add springs there)
        for g in mpc_dependent_grids
            delete!(node_K, g)
        end

        if !isempty(node_K)
            # Generate CELAS2 cards (3 springs per node: DOFs 4, 5, 6)
            # connected to a fixed anchor GRID node (MYSTRAN doesn't support G2=0 grounded springs)
            eid_start = k6rot_max_eid + 1
            eid = eid_start
            spring_count = 0

            # Create anchor grid at origin, fully constrained
            max_nid = maximum(keys(k6rot_grids))
            anchor_nid = max_nid + 999999

            # Find SPC set ID from case control
            anchor_spc_sid = 0
            for ccl in cc_lines
                m_spc = match(r"(?i)^\s*SPC\s*=\s*(\d+)", ccl)
                if m_spc !== nothing
                    anchor_spc_sid = parse(Int, m_spc.captures[1])
                end
            end
            # After SPCADD resolution, the SPC set might have been renumbered
            if haskey(spcadd_map, anchor_spc_sid)
                anchor_spc_sid = spcadd_map[anchor_spc_sid]
            end
            if anchor_spc_sid <= 0; anchor_spc_sid = 1; end

            push!(k6rot_bulk_lines, "\$")
            push!(k6rot_bulk_lines, "\$ --- K6ROT drilling springs (added by converter, K6ROT=$(k6rot_value)) ---")
            push!(k6rot_bulk_lines, rpad("GRID", 8) * lpad(string(anchor_nid), 8) *
                  "              0.      0.      0.")
            push!(k6rot_bulk_lines, rpad("SPC1", 8) * lpad(string(anchor_spc_sid), 8) *
                  lpad("123456", 8) * lpad(string(anchor_nid), 8))

            # Format float into ≤8 chars for NASTRAN fixed-field (e.g., "1.234+3" shorthand)
            function _fmt8(x)
                x == 0.0 && return "      0."
                s = @sprintf("%.3E", x)   # e.g., "1.234E+03"
                # Convert to Nastran shorthand: "1.234+3"
                m = match(r"^([+-]?\d+\.\d+)E([+-])(\d+)$", s)
                if m !== nothing
                    mant = m.captures[1]; sign = m.captures[2]; exp_s = string(parse(Int, m.captures[3]))
                    ns = mant * sign * exp_s
                    return length(ns) <= 8 ? lpad(ns, 8) : lpad(s[1:min(8,length(s))], 8)
                end
                return lpad(s[1:min(8,length(s))], 8)
            end

            for (nid, K_val) in sort(collect(node_K), by=x->x[1])
                K_val <= 0 && continue
                # Split stiffness across 3 rotational DOFs
                K_per_dof = K_val / 3.0
                K_str = _fmt8(K_per_dof)
                # 3 CELAS2 cards (fixed 8-char format): EID, K, G1, C1, G2=anchor, C2=same_dof
                for dof in (4, 5, 6)
                    push!(k6rot_bulk_lines, rpad("CELAS2", 8) * lpad(string(eid), 8) *
                          K_str * lpad(string(nid), 8) *
                          lpad(string(dof), 8) * lpad(string(anchor_nid), 8) *
                          lpad(string(dof), 8))
                    eid += 1
                end
                spring_count += 1
            end

            println("  K6ROT=$(k6rot_value): Injecting $(spring_count * 3) CELAS2 drilling springs at $(spring_count) shell nodes (anchor NID=$(anchor_nid))")
        end
    elseif k6rot_value > 0
        println("  K6ROT=$(k6rot_value) detected but no shell elements found — skipping drilling springs")
    end

    out_bulk = String[]
    has_param_sollib  = false
    has_param_grdpnt  = false
    has_param_post    = false
    has_param_autospc = false
    has_param_bailout = false
    has_debug_192     = false
    has_debug_200     = false

    for line in bulk_lines
        n = norm(line)
        if startswith(n, "PARAM") && occursin("SOLLIB", n)
            has_param_sollib = true
        elseif startswith(n, "PARAM") && occursin("GRDPNT", n)
            has_param_grdpnt = true
        elseif startswith(n, "PARAM") && occursin("POST", n)
            has_param_post = true
        elseif startswith(n, "PARAM") && occursin("AUTOSPC", n)
            has_param_autospc = true
            # Force AUTOSPC YES for MYSTRAN (overrides NO)
            push!(out_bulk, "PARAM   AUTOSPC  YES")
            continue
        elseif startswith(n, "PARAM") && occursin("BAILOUT", n)
            has_param_bailout = true
            # Force BAILOUT -1 to continue past non-fatal errors
            push!(out_bulk, "PARAM   BAILOUT  -1")
            continue
        elseif startswith(n, "DEBUG") && occursin("192", n)
            has_debug_192 = true
        elseif startswith(n, "DEBUG") && occursin("200", n)
            has_debug_200 = true
        end

        # Drop SPCADD cards (resolved by renumbering SPC1 cards)
        if startswith(n, "SPCADD")
            continue
        end

        # Drop PARAM K6ROT (replaced by CELAS drilling springs)
        if startswith(n, "PARAM") && occursin("K6ROT", n) && k6rot_value > 0
            continue
        end

        # Handle SPC/SPC1 cards: renumber for SPCADD, reformat THRU, filter MPC grids
        if startswith(n, "SPC1") || (startswith(n, "SPC") && !startswith(n, "SPC1") &&
                                     !startswith(n, "SPCADD") && !startswith(n, "SPCFORCE"))
            # Renumber SPC set ID if referenced by SPCADD
            if !isempty(spcadd_map)
                f = bdf_fields(line)
                if length(f) >= 2
                    sid = parse_bdf_int(f[2])
                    if haskey(spcadd_map, sid)
                        new_sid = spcadd_map[sid]
                        # Replace the SID field in the line
                        line = rpad(f[1], 8) * lpad(string(new_sid), 8) * line[17:end]
                    end
                end
            end

            if startswith(uppercase(strip(line)), "SPC1")
                # Reformat mixed THRU cards
                spc_lines = occursin("THRU", uppercase(line)) ? reformat_spc1_thru(line) : [line]
                # Filter out MPC-dependent grids from each SPC1 card
                for spc_line in spc_lines
                    filtered = filter_spc1_mpc_grids(spc_line, mpc_dependent_grids)
                    append!(out_bulk, filtered)
                end
            else
                push!(out_bulk, line)
            end
        else
            push!(out_bulk, line)
        end
    end

    # Add missing PARAM/DEBUG cards before ENDDATA
    params_to_add = String[]
    if !has_param_sollib;  push!(params_to_add, MYSTRAN_PARAMS[1]); end
    if !has_param_grdpnt;  push!(params_to_add, MYSTRAN_PARAMS[2]); end
    if !has_param_post;    push!(params_to_add, MYSTRAN_PARAMS[3]); end
    if !has_param_autospc; push!(params_to_add, MYSTRAN_PARAMS[4]); end
    if !has_param_bailout; push!(params_to_add, MYSTRAN_PARAMS[5]); end
    if !has_debug_192;     push!(params_to_add, MYSTRAN_DEBUG[1]); end
    if !has_debug_200;     push!(params_to_add, MYSTRAN_DEBUG[2]); end

    if !isempty(params_to_add)
        push!(out_bulk, "\$")
        push!(out_bulk, "\$ --- MYSTRAN parameters (added by converter) ---")
        append!(out_bulk, params_to_add)
    end

    # Append K6ROT drilling springs (if generated)
    if !isempty(k6rot_bulk_lines)
        append!(out_bulk, k6rot_bulk_lines)
    end

    # ── Phase 5: assemble output ──
    output_lines = String[]
    append!(output_lines, out_exec)
    append!(output_lines, out_cc)
    push!(output_lines, "BEGIN BULK")
    append!(output_lines, out_bulk)
    push!(output_lines, "ENDDATA")

    return join(output_lines, "\n") * "\n"
end

# ─────────────────────────────────────────────────
# MYSTRAN runner
# ─────────────────────────────────────────────────

"""
    run_mystran(mystran_bdf::String, mystran_exe::String) -> Bool

Runs the MYSTRAN executable on the given .bdf file.
MYSTRAN reads the input filename from stdin.
Returns true if the run completed (exit code 0).
"""
function run_mystran(mystran_bdf::String, mystran_exe::String)
    if !isfile(mystran_exe)
        error("MYSTRAN executable not found: $mystran_exe")
    end
    if !isfile(mystran_bdf)
        error("MYSTRAN input file not found: $mystran_bdf")
    end

    bdf_dir  = dirname(abspath(mystran_bdf))
    bdf_name = basename(mystran_bdf)

    println("  Running MYSTRAN on: $bdf_name")
    println("  Working directory:  $bdf_dir")

    # MYSTRAN reads input filename from stdin and writes outputs to the cwd
    proc = run(pipeline(`$mystran_exe`, stdin=IOBuffer(bdf_name), stdout=devnull, stderr=devnull),
               wait=true)

    return proc.exitcode == 0
end

"""
    collect_results(mystran_bdf::String, results_dir::String)

Moves all MYSTRAN output files (matching the base name of the .bdf)
from the working directory into `results_dir`.
"""
function collect_results(mystran_bdf::String, results_dir::String)
    bdf_dir  = dirname(abspath(mystran_bdf))
    bdf_base = replace(basename(mystran_bdf), r"\.[^.]+$" => "")

    mkpath(results_dir)

    moved = 0
    for f in readdir(bdf_dir)
        fpath = joinpath(bdf_dir, f)
        isfile(fpath) || continue

        fname_base = replace(f, r"\.[^.]+$" => "")
        fname_ext  = uppercase(something(match(r"\.[^.]+$", f), (match="",)).match)

        if fname_base == bdf_base && fname_ext in MYSTRAN_OUTPUT_EXTS
            dest = joinpath(results_dir, f)
            mv(fpath, dest; force=true)
            moved += 1
        end
    end
    println("  Collected $moved output files into: $(basename(results_dir))/")
end

# ─────────────────────────────────────────────────
# F06 Validation
# ─────────────────────────────────────────────────

"""
    validate_f06(f06_path) -> (success::Bool, errors::Vector{String})

Checks a MYSTRAN .F06 file for successful completion.
Returns (true, []) if MYSTRAN END is found, or (false, error_messages) otherwise.
"""
function validate_f06(f06_path::String)
    if !isfile(f06_path)
        return (false, ["F06 file not found: $f06_path"])
    end

    raw = read(f06_path)
    cleaned = UInt8[b < 0x20 && b != UInt8('\n') && b != UInt8('\r') ? UInt8(' ') :
                    b > 0x7E ? UInt8(' ') : b for b in raw]
    content = String(cleaned)
    lines = split(content, '\n')

    found_mystran_end = false
    errors = String[]

    for line in lines
        u = uppercase(strip(String(line)))
        if occursin("MYSTRAN END", u)
            found_mystran_end = true
        end
        if startswith(u, "*ERROR") || startswith(u, "* ERROR")
            push!(errors, strip(String(line)))
        end
        if occursin("PROCESSING TERMINATED", u)
            push!(errors, strip(String(line)))
        end
    end

    return (found_mystran_end, errors)
end

# ─────────────────────────────────────────────────
# BDF Geometry Parser
# ─────────────────────────────────────────────────

"""Parse a NASTRAN float string, handling '1.234-5' shorthand for '1.234E-5'."""
function parse_nastran_float(s::AbstractString)
    s = strip(s)
    isempty(s) && return 0.0
    val = tryparse(Float64, s)
    val !== nothing && return val
    # Handle "1.234-5" -> "1.234E-5"
    for i in length(s):-1:2
        c = s[i]; prev = s[i-1]
        if (c == '+' || c == '-') && prev != 'E' && prev != 'e'
            return parse(Float64, s[1:i-1] * "E" * s[i:end])
        end
    end
    return 0.0
end

"""Parse an integer from a BDF field, returning 0 if blank/invalid."""
function parse_bdf_int(s::AbstractString)
    s = strip(s)
    isempty(s) && return 0
    val = tryparse(Int, s)
    return val !== nothing ? val : 0
end

"""Split a BDF line into fixed-width 8-char fields, or free-format comma fields."""
function bdf_fields(line::AbstractString)
    if occursin(',', line)
        return [strip(f) for f in split(line, ',')]
    else
        # Fixed 8-char field format
        fields = String[]
        for i in 1:8:length(line)
            j = min(i+7, length(line))
            push!(fields, strip(line[i:j]))
        end
        return fields
    end
end

"""
    parse_bdf_geometry(bdf_path) -> NamedTuple

Parses GRID, element, property, constraint, and load cards from a NASTRAN .bdf file.
Also parses case-control section for subcase SPC/LOAD IDs.
Returns sorted element/node lists plus v3 constraint/load data for JFEM binary output.
"""
function parse_bdf_geometry(bdf_path::String)
    lines = readlines(bdf_path; keep=false)

    # Raw storage — structural elements
    grid_map = Dict{Int, NTuple{3,Float64}}()   # nid => (x,y,z)
    quad_raw = Dict{Int, Tuple{Int,NTuple{4,Int}}}()  # eid => (pid, (g1,g2,g3,g4))
    tria_raw = Dict{Int, Tuple{Int,NTuple{3,Int}}}()  # eid => (pid, (g1,g2,g3))
    bar_raw  = Dict{Int, Tuple{Int,Int,Int}}()         # eid => (pid, ga, gb)
    rod_raw  = Dict{Int, Tuple{Int,Int,Int}}()         # eid => (pid, ga, gb)

    # Properties
    pshell_t = Dict{Int, Float64}()  # pid => thickness
    pbar_a   = Dict{Int, Float64}()  # pid => area (PBAR, PBARL)
    prod_a   = Dict{Int, Float64}()  # pid => area (PROD)
    pcomp_t  = Dict{Int, Float64}()  # pid => total thickness

    # v3: Constraint & load raw storage
    celas_raw = Dict{Int, NTuple{5,Int}}()      # eid => (pid, g1, c1, g2, c2)
    pelas_k   = Dict{Int, Float64}()             # pid => stiffness K
    rbe2_raw  = Vector{NamedTuple{(:eid,:gn,:cm,:slaves),Tuple{Int,Int,Int,Vector{Int}}}}()
    rbe3_raw  = Vector{NamedTuple{(:eid,:refgrid,:refc,:deps),Tuple{Int,Int,Int,Vector{Int}}}}()
    spc1_raw  = Vector{NamedTuple{(:sid,:c,:nids),Tuple{Int,String,Vector{Int}}}}()
    spcadd_map = Dict{Int, Vector{Int}}()        # target_sid => [sub_sids]
    force_raw  = Vector{NamedTuple{(:sid,:gid,:cid,:mag,:dir),Tuple{Int,Int,Int,Float64,Vector{Float64}}}}()
    moment_raw = Vector{NamedTuple{(:sid,:gid,:cid,:mag,:dir),Tuple{Int,Int,Int,Float64,Vector{Float64}}}}()
    load_combos = Vector{NamedTuple{(:sid,:S,:comps),Tuple{Int,Float64,Vector{Tuple{Float64,Int}}}}}()
    cord2r_raw = Dict{Int, NTuple{9,Float64}}()  # cid => (A1,A2,A3, B1,B2,B3, C1,C2,C3)

    # ── Parse case control section for subcase SPC/LOAD IDs ──
    # subcase_info: maps subcase_id => (spc_id, load_id)
    subcase_info = Dict{Int, NamedTuple{(:spc_id,:load_id),Tuple{Int,Int}}}()
    cend_idx = 0; begin_bulk_idx = 0
    for (i, line) in enumerate(lines)
        n = uppercase(strip(line))
        if n == "CEND"; cend_idx = i; end
        if startswith(n, "BEGIN") && occursin("BULK", n); begin_bulk_idx = i; end
    end

    if cend_idx > 0 && begin_bulk_idx > cend_idx
        cc_lines = lines[cend_idx+1 : begin_bulk_idx-1]
        current_sc = 0
        global_spc = 0; global_load = 0
        local_spc = Dict{Int,Int}(); local_load = Dict{Int,Int}()
        for ccl in cc_lines
            n = uppercase(strip(ccl))
            (isempty(n) || startswith(n, "\$")) && continue
            m_sc = match(r"^SUBCASE\s+(\d+)", n)
            if m_sc !== nothing
                current_sc = parse(Int, m_sc.captures[1])
                continue
            end
            m_spc = match(r"^SPC\s*=\s*(\d+)", n)
            if m_spc !== nothing
                val = parse(Int, m_spc.captures[1])
                if current_sc == 0; global_spc = val; else; local_spc[current_sc] = val; end
                continue
            end
            m_ld = match(r"^LOAD\s*=\s*(\d+)", n)
            if m_ld !== nothing
                val = parse(Int, m_ld.captures[1])
                if current_sc == 0; global_load = val; else; local_load[current_sc] = val; end
                continue
            end
        end
        # Build subcase_info: merge global with local overrides
        all_scs = unique(vcat(collect(keys(local_spc)), collect(keys(local_load))))
        if isempty(all_scs); all_scs = [1]; end
        for sc in all_scs
            spc_val  = get(local_spc, sc, global_spc)
            load_val = get(local_load, sc, global_load)
            subcase_info[sc] = (spc_id=spc_val, load_id=load_val)
        end
        # If no subcases were explicitly defined but global SPC/LOAD exist
        if isempty(subcase_info) && (global_spc > 0 || global_load > 0)
            subcase_info[1] = (spc_id=global_spc, load_id=global_load)
        end
    end

    # ── Parse bulk data ──
    in_bulk = false
    for line in lines
        n = uppercase(strip(line))
        if startswith(n, "BEGIN") && occursin("BULK", n)
            in_bulk = true; continue
        end
        n == "ENDDATA" && break
        !in_bulk && continue
        s = lstrip(line)
        (isempty(s) || startswith(s, "\$")) && continue

        f = bdf_fields(line)
        isempty(f) && continue
        card = uppercase(f[1])

        if card == "GRID" && length(f) >= 6
            nid = parse_bdf_int(f[2])
            x = parse_nastran_float(length(f) >= 4 ? f[4] : "0")
            y = parse_nastran_float(length(f) >= 5 ? f[5] : "0")
            z = parse_nastran_float(length(f) >= 6 ? f[6] : "0")
            nid > 0 && (grid_map[nid] = (x, y, z))

        elseif card == "CQUAD4" && length(f) >= 7
            eid = parse_bdf_int(f[2]); pid = parse_bdf_int(f[3])
            g1 = parse_bdf_int(f[4]); g2 = parse_bdf_int(f[5])
            g3 = parse_bdf_int(f[6]); g4 = parse_bdf_int(f[7])
            eid > 0 && (quad_raw[eid] = (pid, (g1, g2, g3, g4)))

        elseif card == "CTRIA3" && length(f) >= 6
            eid = parse_bdf_int(f[2]); pid = parse_bdf_int(f[3])
            g1 = parse_bdf_int(f[4]); g2 = parse_bdf_int(f[5]); g3 = parse_bdf_int(f[6])
            eid > 0 && (tria_raw[eid] = (pid, (g1, g2, g3)))

        elseif card == "CBAR" && length(f) >= 5
            eid = parse_bdf_int(f[2]); pid = parse_bdf_int(f[3])
            ga = parse_bdf_int(f[4]); gb = parse_bdf_int(f[5])
            eid > 0 && (bar_raw[eid] = (pid, ga, gb))

        elseif card == "CROD" && length(f) >= 5
            eid = parse_bdf_int(f[2]); pid = parse_bdf_int(f[3])
            ga = parse_bdf_int(f[4]); gb = parse_bdf_int(f[5])
            eid > 0 && (rod_raw[eid] = (pid, ga, gb))

        elseif card == "PSHELL" && length(f) >= 4
            pid = parse_bdf_int(f[2])
            t = parse_nastran_float(f[4])
            pid > 0 && (pshell_t[pid] = t)

        elseif card == "PROD" && length(f) >= 4
            pid = parse_bdf_int(f[2])
            a = parse_nastran_float(f[4])
            pid > 0 && (prod_a[pid] = a)

        elseif card == "PBAR" && length(f) >= 4
            pid = parse_bdf_int(f[2])
            a = parse_nastran_float(f[4])
            pid > 0 && (pbar_a[pid] = a)

        elseif card == "PBARL"
            pid = parse_bdf_int(f[2])
            pid > 0 && !haskey(pbar_a, pid) && (pbar_a[pid] = 0.0)

        # ── v3: Constraint & load cards ──
        elseif card == "CELAS1" && length(f) >= 7
            eid = parse_bdf_int(f[2]); pid = parse_bdf_int(f[3])
            g1 = parse_bdf_int(f[4]); c1 = parse_bdf_int(f[5])
            g2 = parse_bdf_int(f[6]); c2 = parse_bdf_int(f[7])
            eid > 0 && (celas_raw[eid] = (pid, g1, c1, g2, c2))

        elseif card == "PELAS" && length(f) >= 3
            pid = parse_bdf_int(f[2])
            k = parse_nastran_float(f[3])
            pid > 0 && (pelas_k[pid] = k)

        elseif card == "RBE2" && length(f) >= 5
            eid = parse_bdf_int(f[2]); gn = parse_bdf_int(f[3]); cm = parse_bdf_int(f[4])
            slaves = Int[]
            for k in 5:length(f)
                gid = parse_bdf_int(f[k])
                gid > 0 && push!(slaves, gid)
            end
            eid > 0 && !isempty(slaves) && push!(rbe2_raw, (eid=eid, gn=gn, cm=cm, slaves=slaves))

        elseif card == "RBE3" && length(f) >= 5
            eid = parse_bdf_int(f[2])
            refgrid = parse_bdf_int(f[4]); refc = parse_bdf_int(f[5])
            # Dependent grids: skip WT/C fields, collect grid IDs from field 9 onward
            deps = Int[]
            k = 6
            while k <= length(f)
                tok = strip(f[k])
                if isempty(tok); k += 1; continue; end
                gid = tryparse(Int, tok)
                if gid !== nothing && gid > 0
                    # Could be WT(float), C(dof int like 123456), or grid ID
                    # Heuristic: grid IDs are large (> 6), WT/C are small or have decimal
                    if occursin('.', tok)
                        # This is a WT (weight) float — skip it
                    elseif gid <= 6 || (gid >= 10 && gid <= 999999 && all(c -> c in "0123456", tok))
                        # Looks like a DOF component (e.g., 123, 123456) — skip
                    else
                        push!(deps, gid)
                    end
                end
                k += 1
            end
            eid > 0 && refgrid > 0 && push!(rbe3_raw, (eid=eid, refgrid=refgrid, refc=refc, deps=deps))

        elseif card == "SPC1" && length(f) >= 4
            sid = parse_bdf_int(f[2]); c_str = strip(f[3])
            nids = Int[]
            # Check for THRU
            thru_pos = 0
            for k in 4:length(f)
                if uppercase(strip(f[k])) == "THRU"; thru_pos = k; break; end
            end
            if thru_pos > 0 && thru_pos + 1 <= length(f)
                g_start = parse_bdf_int(f[thru_pos - 1 >= 4 ? thru_pos - 1 : 4])
                g_end   = parse_bdf_int(f[thru_pos + 1])
                # Collect individual grids before THRU range start
                for k in 4:(thru_pos - 2)
                    gid = parse_bdf_int(f[k])
                    gid > 0 && push!(nids, gid)
                end
                # Expand THRU range
                if g_start > 0 && g_end >= g_start
                    append!(nids, g_start:g_end)
                end
            else
                for k in 4:length(f)
                    gid = parse_bdf_int(f[k])
                    gid > 0 && push!(nids, gid)
                end
            end
            sid > 0 && !isempty(nids) && push!(spc1_raw, (sid=sid, c=c_str, nids=nids))

        elseif card == "SPCADD" && length(f) >= 3
            target = parse_bdf_int(f[2])
            subs = Int[]
            for k in 3:length(f)
                sub = parse_bdf_int(f[k])
                sub > 0 && push!(subs, sub)
            end
            target > 0 && !isempty(subs) && (spcadd_map[target] = subs)

        elseif card == "FORCE" && length(f) >= 8
            sid = parse_bdf_int(f[2]); gid = parse_bdf_int(f[3]); cid = parse_bdf_int(f[4])
            mag = parse_nastran_float(f[5])
            n1 = parse_nastran_float(f[6]); n2 = parse_nastran_float(f[7]); n3 = parse_nastran_float(f[8])
            sid > 0 && gid > 0 && push!(force_raw, (sid=sid, gid=gid, cid=cid, mag=mag, dir=[n1, n2, n3]))

        elseif card == "MOMENT" && length(f) >= 8
            sid = parse_bdf_int(f[2]); gid = parse_bdf_int(f[3]); cid = parse_bdf_int(f[4])
            mag = parse_nastran_float(f[5])
            n1 = parse_nastran_float(f[6]); n2 = parse_nastran_float(f[7]); n3 = parse_nastran_float(f[8])
            sid > 0 && gid > 0 && push!(moment_raw, (sid=sid, gid=gid, cid=cid, mag=mag, dir=[n1, n2, n3]))

        elseif card == "LOAD" && length(f) >= 5
            sid = parse_bdf_int(f[2]); overall_s = parse_nastran_float(f[3])
            comps = Tuple{Float64,Int}[]
            k = 4
            while k + 1 <= length(f)
                si = parse_nastran_float(f[k]); li = parse_bdf_int(f[k+1])
                li > 0 && push!(comps, (si, li))
                k += 2
            end
            sid > 0 && !isempty(comps) && push!(load_combos, (sid=sid, S=overall_s, comps=comps))

        elseif card == "CORD2R" && length(f) >= 10
            cid = parse_bdf_int(f[2])
            a1 = parse_nastran_float(f[4]); a2 = parse_nastran_float(f[5]); a3 = parse_nastran_float(f[6])
            b1 = parse_nastran_float(f[7]); b2 = parse_nastran_float(f[8]); b3 = parse_nastran_float(f[9])
            # C point is on the continuation (field 10+)
            c1 = length(f) >= 10 ? parse_nastran_float(f[10]) : 0.0
            c2 = length(f) >= 11 ? parse_nastran_float(f[11]) : 0.0
            c3 = length(f) >= 12 ? parse_nastran_float(f[12]) : 0.0
            cid > 0 && (cord2r_raw[cid] = (a1, a2, a3, b1, b2, b3, c1, c2, c3))
        end
    end

    # Build sorted output arrays matching JFEM format
    # Nodes: sorted by nid
    node_ids = sort(collect(keys(grid_map)))
    nodes = [(nid, grid_map[nid]...) for nid in node_ids]

    # Helper to get thickness/area from property
    get_thickness(pid) = get(pshell_t, pid, get(pcomp_t, pid, 0.0))
    get_bar_area(pid)  = get(pbar_a, pid, 0.0)
    get_rod_area(pid)  = get(prod_a, pid, 0.0)

    # Elements: sorted by eid — (eid, pid, nodes..., thickness_or_area)
    quad_eids = sort(collect(keys(quad_raw)))
    quads = [(eid, quad_raw[eid][1], quad_raw[eid][2]..., Float32(get_thickness(quad_raw[eid][1])))
             for eid in quad_eids]

    tria_eids = sort(collect(keys(tria_raw)))
    trias = [(eid, tria_raw[eid][1], tria_raw[eid][2]..., Float32(get_thickness(tria_raw[eid][1])))
             for eid in tria_eids]

    bar_eids = sort(collect(keys(bar_raw)))
    bars = [(eid, bar_raw[eid][1], bar_raw[eid][2], bar_raw[eid][3], Float32(get_bar_area(bar_raw[eid][1])))
            for eid in bar_eids]

    rod_eids_sorted = sort(collect(keys(rod_raw)))
    rods = [(eid, rod_raw[eid][1], rod_raw[eid][2], rod_raw[eid][3], Float32(get_rod_area(rod_raw[eid][1])))
            for eid in rod_eids_sorted]

    return (nodes=nodes, quads=quads, trias=trias, bars=bars, rods=rods,
            node_ids=node_ids, quad_eids=quad_eids, tria_eids=tria_eids,
            bar_eids=bar_eids, rod_eids=rod_eids_sorted,
            # v3 data
            celas_raw=celas_raw, pelas_k=pelas_k,
            rbe2_raw=rbe2_raw, rbe3_raw=rbe3_raw,
            spc1_raw=spc1_raw, spcadd_map=spcadd_map,
            force_raw=force_raw, moment_raw=moment_raw,
            load_combos=load_combos, cord2r_raw=cord2r_raw,
            subcase_info=subcase_info)
end

# ─────────────────────────────────────────────────
# v3 Helper Functions — Constraints, Loads, CORD2R
# ─────────────────────────────────────────────────

"""Build 3×3 rotation matrix from CORD2R A (origin), B (z-axis point), C (xz-plane point)."""
function build_cord2r_matrix(A, B, C)
    z = B .- A; z = z ./ max(sqrt(sum(z .^ 2)), 1e-30)
    xp = C .- A; xp = xp ./ max(sqrt(sum(xp .^ 2)), 1e-30)
    y = [z[2]*xp[3]-z[3]*xp[2], z[3]*xp[1]-z[1]*xp[3], z[1]*xp[2]-z[2]*xp[1]]
    yn = max(sqrt(sum(y .^ 2)), 1e-30); y = y ./ yn
    x = [y[2]*z[3]-y[3]*z[2], y[3]*z[1]-y[1]*z[3], y[1]*z[2]-y[2]*z[1]]
    return hcat(x, y, z)  # 3×3 rotation matrix: columns are x,y,z unit vectors
end

"""Transform vector from local CORD2R to global coordinates."""
function coord_transform(cord2r_raw, cid, vec)
    cid == 0 && return vec
    !haskey(cord2r_raw, cid) && return vec
    c = cord2r_raw[cid]
    A = [c[1], c[2], c[3]]; B = [c[4], c[5], c[6]]; C = [c[7], c[8], c[9]]
    R = build_cord2r_matrix(A, B, C)
    return R * vec
end

"""Build CELAS, RBE2, RBE3 tables from parsed BDF geometry data."""
function build_constraint_tables(geom)
    # CELAS1: (eid, g1, c1, g2, c2, K)
    jfem_celas = Tuple{Int,Int,Int,Int,Int,Float32}[]
    for (eid, (pid, g1, c1, g2, c2)) in geom.celas_raw
        K = Float32(get(geom.pelas_k, pid, 0.0))
        push!(jfem_celas, (eid, g1, c1, g2, c2, K))
    end
    sort!(jfem_celas, by=x->x[1])

    # RBE2: (eid, gn, cm, slaves[])
    jfem_rbe2s = [(eid=r.eid, gn=r.gn, cm=r.cm, slaves=r.slaves) for r in geom.rbe2_raw]
    sort!(jfem_rbe2s, by=x->x.eid)

    # RBE3: (eid, refgrid, refc, deps[])
    jfem_rbe3s = [(eid=r.eid, refgrid=r.refgrid, refc=r.refc, deps=r.deps) for r in geom.rbe3_raw]
    sort!(jfem_rbe3s, by=x->x.eid)

    return jfem_celas, jfem_rbe2s, jfem_rbe3s
end

"""Collect SPC constraints for a subcase, resolving SPCADD."""
function collect_spc_for_subcase(geom, spc_id)
    spc_id <= 0 && return Tuple{Int32, UInt32}[]

    # Resolve SPCADD
    active_sids = Set{Int}()
    if haskey(geom.spcadd_map, spc_id)
        union!(active_sids, geom.spcadd_map[spc_id])
    else
        push!(active_sids, spc_id)
    end

    # Collect SPC1 cards matching active SIDs
    spc_nodes = Dict{Int, Set{Int}}()  # nid => set of DOFs
    for spc in geom.spc1_raw
        spc.sid in active_sids || continue
        for nid in spc.nids
            if !haskey(spc_nodes, nid); spc_nodes[nid] = Set{Int}(); end
            for ch in spc.c
                if isdigit(ch); push!(spc_nodes[nid], parse(Int, string(ch))); end
            end
        end
    end

    # Convert to (nid, dof_mask) tuples
    result = Tuple{Int32, UInt32}[]
    for (nid, dofs) in sort(collect(spc_nodes), by=x->x[1])
        mask = 0
        for d in sort(collect(dofs))
            mask = mask * 10 + d
        end
        push!(result, (Int32(nid), UInt32(mask)))
    end
    return result
end

"""Recursively collect FORCE/MOMENT cards, handling LOAD combos."""
function _collect_loads_recursive!(geom, sid, scale, forces_acc, moments_acc)
    # Collect FORCE cards
    for frc in geom.force_raw
        frc.sid == sid || continue
        dir_len = sqrt(sum(frc.dir .^ 2))
        dir_norm = dir_len > 1e-30 ? frc.dir ./ dir_len : frc.dir
        global_dir = coord_transform(geom.cord2r_raw, frc.cid, dir_norm)
        fvec = global_dir .* (frc.mag * scale)
        if !haskey(forces_acc, frc.gid); forces_acc[frc.gid] = zeros(3); end
        forces_acc[frc.gid] .+= fvec
    end

    # Collect MOMENT cards
    for mom in geom.moment_raw
        mom.sid == sid || continue
        dir_len = sqrt(sum(mom.dir .^ 2))
        dir_norm = dir_len > 1e-30 ? mom.dir ./ dir_len : mom.dir
        global_dir = coord_transform(geom.cord2r_raw, mom.cid, dir_norm)
        mvec = global_dir .* (mom.mag * scale)
        if !haskey(moments_acc, mom.gid); moments_acc[mom.gid] = zeros(3); end
        moments_acc[mom.gid] .+= mvec
    end

    # Recurse through LOAD combos
    for combo in geom.load_combos
        combo.sid == sid || continue
        for (si, li) in combo.comps
            _collect_loads_recursive!(geom, li, scale * combo.S * si, forces_acc, moments_acc)
        end
    end
end

"""Collect applied forces and moments for a subcase's LOAD set."""
function collect_forces_for_subcase(geom, load_id)
    safe_f32(x) = (v = Float32(x); isnan(v) || isinf(v) ? Float32(0) : v)
    force_data = Tuple{Int32, Float32, Float32, Float32}[]
    moment_data = Tuple{Int32, Float32, Float32, Float32}[]
    load_id <= 0 && return force_data, moment_data

    forces_acc = Dict{Int, Vector{Float64}}()
    moments_acc = Dict{Int, Vector{Float64}}()
    _collect_loads_recursive!(geom, load_id, 1.0, forces_acc, moments_acc)

    for (nid, fvec) in sort(collect(forces_acc), by=x->x[1])
        if sqrt(sum(fvec .^ 2)) > 1e-30
            push!(force_data, (Int32(nid), safe_f32(fvec[1]), safe_f32(fvec[2]), safe_f32(fvec[3])))
        end
    end
    for (nid, mvec) in sort(collect(moments_acc), by=x->x[1])
        if sqrt(sum(mvec .^ 2)) > 1e-30
            push!(moment_data, (Int32(nid), safe_f32(mvec[1]), safe_f32(mvec[2]), safe_f32(mvec[3])))
        end
    end

    return force_data, moment_data
end

# ─────────────────────────────────────────────────
# MYSTRAN F06 Parser
# ─────────────────────────────────────────────────

"""Check if a line (with spaces removed) contains a keyword."""
function f06_header_match(line::AbstractString, keyword::String)
    return occursin(keyword, replace(line, " " => ""))
end

"""
    parse_mystran_f06(f06_path, geom) -> Vector of per-subcase result NamedTuples

Parses a MYSTRAN .F06 file and returns displacement, force, and stress data
organized per subcase, matching the JFEM binary layout.
"""
function parse_mystran_f06(f06_path::String, geom)
    raw = read(f06_path)
    # Replace non-ASCII and null bytes with spaces to avoid uppercase() errors
    cleaned = UInt8[b < 0x20 && b != UInt8('\n') && b != UInt8('\r') ? UInt8(' ') :
                    b > 0x7E ? UInt8(' ') : b for b in raw]
    content = String(cleaned)
    lines = split(content, '\n')

    node_ids  = geom.node_ids
    quad_eids = geom.quad_eids
    tria_eids = geom.tria_eids
    bar_eids  = geom.bar_eids
    rod_eids  = geom.rod_eids

    nNodes = length(node_ids)
    nQuads = length(quad_eids)
    nTrias = length(tria_eids)
    nBars  = length(bar_eids)
    nRods  = length(rod_eids)

    # Index maps for fast lookup
    node_idx = Dict(nid => i for (i, nid) in enumerate(node_ids))
    quad_idx = Dict(eid => i for (i, eid) in enumerate(quad_eids))
    tria_idx = Dict(eid => i for (i, eid) in enumerate(tria_eids))
    bar_idx  = Dict(eid => i for (i, eid) in enumerate(bar_eids))
    rod_idx  = Dict(eid => i for (i, eid) in enumerate(rod_eids))

    safe_f32(x) = (v = Float32(x); isnan(v) || isinf(v) ? Float32(0) : v)

    # Per-subcase accumulators
    subcases = Dict{Int, NamedTuple{(:sid,:disp,:shell,:bar,:rod),
                Tuple{Int,Vector{Float32},Vector{Float32},Vector{Float32},Vector{Float32}}}}()

    function get_sc(sid)
        if !haskey(subcases, sid)
            subcases[sid] = (sid=sid,
                disp  = zeros(Float32, nNodes * 6),
                shell = zeros(Float32, (nQuads + nTrias) * 7),
                bar   = zeros(Float32, nBars * 7),
                rod   = zeros(Float32, nRods * 2))
        end
        return subcases[sid]
    end

    current_subcase = 1
    section = :none  # :disp, :bar_forces, :quad_forces, :tria_forces, :rod_forces,
                     # :bar_stresses, :quad_stresses, :tria_stresses, :rod_stresses
    skip_header_lines = 0

    i = 1
    while i <= length(lines)
        line = lines[i]
        clean = replace(line, " " => "")
        uclean = uppercase(clean)

        # Detect subcase
        m_sc = match(r"OUTPUT\s+FOR\s+SUBCASE\s+(\d+)", uppercase(line))
        if m_sc !== nothing
            current_subcase = parse(Int, m_sc.captures[1])
            i += 1; continue
        end

        # Detect section headers (two-line headers for forces/stresses)
        if f06_header_match(uppercase(line), "DISPLACEMENTS")
            section = :disp; skip_header_lines = 3; i += 1; continue
        end

        if f06_header_match(uppercase(line), "ELEMENTENGINEERINGFORCES")
            # Next line has element type
            if i + 1 <= length(lines)
                next_clean = uppercase(replace(lines[i+1], " " => ""))
                if occursin("TYPEBAR", next_clean)
                    section = :bar_forces; skip_header_lines = 2; i += 2; continue
                elseif occursin("TYPEQUAD4", next_clean) || occursin("TYPEQUAD", next_clean)
                    section = :quad_forces; skip_header_lines = 3; i += 2; continue
                elseif occursin("TYPETRIA3", next_clean) || occursin("TYPETRIA", next_clean)
                    section = :tria_forces; skip_header_lines = 3; i += 2; continue
                elseif occursin("TYPEROD", next_clean)
                    section = :rod_forces; skip_header_lines = 2; i += 2; continue
                end
            end
            i += 1; continue
        end

        if f06_header_match(uppercase(line), "ELEMENTSTRESSESINLOCAL")
            if i + 1 <= length(lines)
                next_clean = uppercase(replace(lines[i+1], " " => ""))
                if occursin("TYPEBAR", next_clean)
                    section = :bar_stresses; skip_header_lines = 3; i += 2; continue
                elseif occursin("TYPEQUAD4", next_clean) || occursin("TYPEQUAD", next_clean)
                    section = :quad_stresses; skip_header_lines = 4; i += 2; continue
                elseif occursin("TYPETRIA3", next_clean) || occursin("TYPETRIA", next_clean)
                    section = :tria_stresses; skip_header_lines = 4; i += 2; continue
                elseif occursin("TYPEROD", next_clean)
                    section = :rod_stresses; skip_header_lines = 2; i += 2; continue
                end
            end
            i += 1; continue
        end

        # Skip header/column label lines
        if skip_header_lines > 0
            skip_header_lines -= 1; i += 1; continue
        end

        # End of data section detection
        stripped = strip(line)
        if section != :none
            if isempty(stripped) || startswith(stripped, "---") ||
               startswith(uppercase(stripped), "MAX*") || startswith(uppercase(stripped), "MIN*") ||
               startswith(uppercase(stripped), "ABS*") || startswith(uppercase(stripped), "*FOR") ||
               startswith(uppercase(stripped), "OUTPUT FOR")
                if startswith(uppercase(stripped), "OUTPUT FOR") ||
                   (isempty(stripped) && i+1 <= length(lines) && occursin("OUTPUT FOR", uppercase(lines[i+1])))
                    section = :none
                end
                i += 1; continue
            end
        end

        sc = get_sc(current_subcase)

        # ── Parse displacement data rows ──
        if section == :disp
            fields = split(stripped)
            if length(fields) >= 8
                nid = tryparse(Int, fields[1])
                if nid !== nothing && haskey(node_idx, nid)
                    idx = node_idx[nid]
                    base = (idx - 1) * 6
                    # fields: nid, coord_sys, t1, t2, t3, r1, r2, r3
                    for k in 1:6
                        sc.disp[base + k] = safe_f32(parse_nastran_float(fields[2 + k]))
                    end
                end
            end

        # ── Parse QUAD4/TRIA3 force data rows ──
        elseif section == :quad_forces || section == :tria_forces
            fields = split(stripped)
            # Format: eid, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, Qx, Qy (9 fields)
            if length(fields) >= 9
                eid = tryparse(Int, fields[1])
                if eid !== nothing
                    idx_map = section == :quad_forces ? quad_idx : tria_idx
                    offset = section == :tria_forces ? nQuads : 0
                    if haskey(idx_map, eid)
                        idx = idx_map[eid] + offset
                        base = (idx - 1) * 7
                        # fx=Nxx, fy=Nyy, fxy=Nxy, mx=Mxx, my=Myy, mxy=Mxy
                        sc.shell[base+1] = safe_f32(parse_nastran_float(fields[2]))  # Nxx -> fx
                        sc.shell[base+2] = safe_f32(parse_nastran_float(fields[3]))  # Nyy -> fy
                        sc.shell[base+3] = safe_f32(parse_nastran_float(fields[4]))  # Nxy -> fxy
                        sc.shell[base+4] = safe_f32(parse_nastran_float(fields[5]))  # Mxx -> mx
                        sc.shell[base+5] = safe_f32(parse_nastran_float(fields[6]))  # Myy -> my
                        sc.shell[base+6] = safe_f32(parse_nastran_float(fields[7]))  # Mxy -> mxy
                        # vonmises filled from stresses section
                    end
                end
            end

        # ── Parse QUAD4/TRIA3 stress data rows (2 lines per element) ──
        elseif section == :quad_stresses || section == :tria_stresses
            fields = split(stripped)
            # First line: eid, "CENTER"/"Anywhere", fibre_dist, sx, sy, sxy, angle, major, minor, vm [, txz, tyz]
            if length(fields) >= 10
                eid = tryparse(Int, fields[1])
                if eid !== nothing
                    idx_map = section == :quad_stresses ? quad_idx : tria_idx
                    offset = section == :tria_stresses ? nQuads : 0
                    if haskey(idx_map, eid)
                        idx = idx_map[eid] + offset
                        base = (idx - 1) * 7
                        vm_z1 = abs(parse_nastran_float(fields[10]))
                        # Read second line (z2)
                        if i + 1 <= length(lines)
                            fields2 = split(strip(lines[i+1]))
                            if length(fields2) >= 8
                                vm_z2 = abs(parse_nastran_float(fields2[8]))
                                sc.shell[base+7] = safe_f32(max(vm_z1, vm_z2))
                                i += 1  # skip z2 line
                            else
                                sc.shell[base+7] = safe_f32(vm_z1)
                            end
                        end
                    end
                end
            end

        # ── Parse BAR force data rows ──
        elseif section == :bar_forces
            fields = split(stripped)
            # Format: eid, BendA1, BendA2, BendB1, BendB2, Shear1, Shear2, AxialForce, Torque (9 fields)
            if length(fields) >= 9
                eid = tryparse(Int, fields[1])
                if eid !== nothing && haskey(bar_idx, eid)
                    idx = bar_idx[eid]
                    base = (idx - 1) * 7
                    sc.bar[base+1] = safe_f32(parse_nastran_float(fields[8]))  # axial
                    sc.bar[base+2] = safe_f32(parse_nastran_float(fields[6]))  # shear_1
                    sc.bar[base+3] = safe_f32(parse_nastran_float(fields[7]))  # shear_2
                    sc.bar[base+4] = safe_f32(parse_nastran_float(fields[9]))  # torque
                    sc.bar[base+5] = safe_f32(parse_nastran_float(fields[2]))  # moment_a1 (BendA Plane1)
                    sc.bar[base+6] = safe_f32(parse_nastran_float(fields[3]))  # moment_a2 (BendA Plane2)
                end
            end

        # ── Parse BAR stress data rows (2 lines per element) ──
        elseif section == :bar_stresses
            fields = split(stripped)
            # First line (end_a): eid, SA1, SA2, SA3, SA4, Axial, SA-Max, SA-Min [, MS-T]
            if length(fields) >= 7
                eid = tryparse(Int, fields[1])
                if eid !== nothing && haskey(bar_idx, eid)
                    idx = bar_idx[eid]
                    base = (idx - 1) * 7
                    axial_stress = parse_nastran_float(fields[6])
                    sa = [abs(parse_nastran_float(fields[k])) for k in 2:5]
                    vm = abs(axial_stress)
                    for s in sa; vm = max(vm, s); end
                    # Read second line (end_b)
                    if i + 1 <= length(lines)
                        fields2 = split(strip(lines[i+1]))
                        if length(fields2) >= 4 && tryparse(Int, fields2[1]) === nothing
                            sb = [abs(parse_nastran_float(fields2[k])) for k in 1:4]
                            for s in sb; vm = max(vm, s); end
                            i += 1
                        end
                    end
                    sc.bar[base+7] = safe_f32(vm)
                end
            end

        # ── Parse ROD force data rows (up to 3 elements per line) ──
        elseif section == :rod_forces
            fields = split(stripped)
            # Format: [eid, axial, torque] repeated up to 3x per line
            j = 1
            while j + 2 <= length(fields)
                eid = tryparse(Int, fields[j])
                if eid !== nothing && haskey(rod_idx, eid)
                    idx = rod_idx[eid]
                    base = (idx - 1) * 2
                    sc.rod[base+1] = safe_f32(parse_nastran_float(fields[j+1]))  # axial
                    sc.rod[base+2] = safe_f32(parse_nastran_float(fields[j+2]))  # torque
                end
                j += 3
            end

        # ── Parse ROD stress data rows (up to 2 elements per line) ──
        elseif section == :rod_stresses
            # We don't need rod stresses for JFEM (only axial/torque from forces)
            # Skip - vonmises not stored for rods in JFEM format
        end

        i += 1
    end

    # Return sorted by subcase ID
    sc_ids = sort(collect(keys(subcases)))
    if isempty(sc_ids)
        sc_ids = [1]
        get_sc(1)
    end
    return [subcases[sid] for sid in sc_ids]
end

# ─────────────────────────────────────────────────
# JFEM Binary Writer
# ─────────────────────────────────────────────────

"""
    write_jfem_binary(jfem_path, geom, subcases_data; jfem_celas=[], jfem_rbe2s=[], jfem_rbe3s=[])

Writes a .jfem binary file compatible with postv11.html viewer.
Format version 3 with constraint tables and per-subcase SPC/force/moment data.
"""
function write_jfem_binary(jfem_path::String, geom, subcases_data;
                           jfem_celas=[], jfem_rbe2s=[], jfem_rbe3s=[])
    nNodes = length(geom.nodes)
    nQuads = length(geom.quads)
    nTrias = length(geom.trias)
    nBars  = length(geom.bars)
    nRods  = length(geom.rods)
    nSC    = length(subcases_data)
    nCelas = length(jfem_celas)
    nRBE2  = length(jfem_rbe2s)
    nRBE3  = length(jfem_rbe3s)

    open(jfem_path, "w") do io
        # Magic: 'JFEM'
        write(io, UInt8('J')); write(io, UInt8('F')); write(io, UInt8('E')); write(io, UInt8('M'))
        # Header (v3: extended with constraint counts)
        write(io, UInt32(3))       # version 3
        write(io, UInt32(nNodes))
        write(io, UInt32(nQuads))
        write(io, UInt32(nTrias))
        write(io, UInt32(nBars))
        write(io, UInt32(nRods))
        write(io, UInt32(nSC))
        write(io, UInt32(nCelas))  # v3
        write(io, UInt32(nRBE2))   # v3
        write(io, UInt32(nRBE3))   # v3

        # Node table: nid(i32), x(f32), y(f32), z(f32)
        for (nid, x, y, z) in geom.nodes
            write(io, Int32(nid))
            write(io, Float32(x)); write(io, Float32(y)); write(io, Float32(z))
        end

        # CQUAD4 table: eid(i32), pid(i32), g1-g4(i32), thickness(f32)
        for (eid, pid, g1, g2, g3, g4, t) in geom.quads
            write(io, Int32(eid)); write(io, Int32(pid))
            write(io, Int32(g1)); write(io, Int32(g2)); write(io, Int32(g3)); write(io, Int32(g4))
            write(io, t)
        end

        # CTRIA3 table: eid(i32), pid(i32), g1-g3(i32), thickness(f32)
        for (eid, pid, g1, g2, g3, t) in geom.trias
            write(io, Int32(eid)); write(io, Int32(pid))
            write(io, Int32(g1)); write(io, Int32(g2)); write(io, Int32(g3))
            write(io, t)
        end

        # CBAR table: eid(i32), pid(i32), ga(i32), gb(i32), area(f32)
        for (eid, pid, ga, gb, a) in geom.bars
            write(io, Int32(eid)); write(io, Int32(pid))
            write(io, Int32(ga)); write(io, Int32(gb))
            write(io, a)
        end

        # CROD table: eid(i32), pid(i32), ga(i32), gb(i32), area(f32)
        for (eid, pid, ga, gb, a) in geom.rods
            write(io, Int32(eid)); write(io, Int32(pid))
            write(io, Int32(ga)); write(io, Int32(gb))
            write(io, a)
        end

        # v3: CELAS table: eid(i32), g1(i32), c1(i32), g2(i32), c2(i32), stiffness(f32), pad(f32)
        for (eid, g1, c1, g2, c2, K_stiff) in jfem_celas
            write(io, Int32(eid)); write(io, Int32(g1)); write(io, Int32(c1))
            write(io, Int32(g2)); write(io, Int32(c2)); write(io, K_stiff); write(io, Float32(0))
        end

        # v3: RBE2 table (variable-length): eid(i32), gn(i32), cm(i32), nSlaves(u32), [slave_nid(i32) × nSlaves]
        for rbe in jfem_rbe2s
            write(io, Int32(rbe.eid)); write(io, Int32(rbe.gn)); write(io, Int32(rbe.cm))
            write(io, UInt32(length(rbe.slaves)))
            for s in rbe.slaves; write(io, Int32(s)); end
        end

        # v3: RBE3 table (variable-length): eid(i32), refgrid(i32), refc(i32), nDep(u32), [dep_nid(i32) × nDep]
        for rbe in jfem_rbe3s
            write(io, Int32(rbe.eid)); write(io, Int32(rbe.refgrid)); write(io, Int32(rbe.refc))
            write(io, UInt32(length(rbe.deps)))
            for d in rbe.deps; write(io, Int32(d)); end
        end

        # Per-subcase data
        for sc in subcases_data
            write(io, UInt32(sc.sid))
            write(io, sc.disp)    # nNodes * 6 Float32
            write(io, sc.shell)   # (nQuads + nTrias) * 7 Float32
            write(io, sc.bar)     # nBars * 7 Float32
            write(io, sc.rod)     # nRods * 2 Float32

            # v3: SPC data
            spc = hasproperty(sc, :spc) ? sc.spc : Tuple{Int32, UInt32}[]
            write(io, UInt32(length(spc)))
            for (nid, mask) in spc
                write(io, nid); write(io, mask)
            end

            # v3: Applied forces
            forces = hasproperty(sc, :forces) ? sc.forces : Tuple{Int32, Float32, Float32, Float32}[]
            write(io, UInt32(length(forces)))
            for (nid, fx, fy, fz) in forces
                write(io, nid); write(io, fx); write(io, fy); write(io, fz)
            end

            # v3: Applied moments
            moments = hasproperty(sc, :moments) ? sc.moments : Tuple{Int32, Float32, Float32, Float32}[]
            write(io, UInt32(length(moments)))
            for (nid, mx, my, mz) in moments
                write(io, nid); write(io, mx); write(io, my); write(io, mz)
            end
        end
    end

    println("  JFEM v3: $nNodes nodes, $(nQuads)Q+$(nTrias)T shells, $nBars bars, $nRods rods, $(nCelas) springs, $(nRBE2) RBE2, $(nRBE3) RBE3, $nSC subcases")
end

# ─────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────

function main()
    if length(ARGS) < 1
        println("""
        Usage: julia nastran_to_mystran.jl <input.bdf> [mystran_exe_path]

        Converts a NASTRAN .bdf to MYSTRAN format, runs MYSTRAN, and
        collects results into <model_name>_MYST_RESULTS/

        Arguments:
          input.bdf         Path to the NASTRAN .bdf input file
          mystran_exe_path  (optional) Path to the MYSTRAN executable
        """)
        return
    end

    input_bdf  = abspath(ARGS[1])
    mystran_exe = length(ARGS) >= 2 ? abspath(ARGS[2]) : ""

    # Resolve MYSTRAN executable
    if isempty(mystran_exe)
        mystran_exe = get(ENV, "MYSTRAN_EXE", "")
    end
    if isempty(mystran_exe) || !isfile(mystran_exe)
        mystran_exe = abspath(DEFAULT_MYSTRAN_EXE)
    end

    if !isfile(input_bdf)
        error("Input file not found: $input_bdf")
    end
    if !isfile(mystran_exe)
        error("MYSTRAN executable not found: $mystran_exe\n" *
              "Set MYSTRAN_EXE environment variable or pass as second argument.")
    end

    input_dir  = dirname(input_bdf)
    input_name = basename(input_bdf)
    model_name = replace(input_name, r"\.[^.]+$" => "")

    # Output paths
    mystran_bdf = joinpath(input_dir, model_name * ".MYSTRAN.bdf")
    results_dir = joinpath(input_dir, model_name * "_MYST_RESULTS")

    println("╔══════════════════════════════════════════════════╗")
    println("║         NASTRAN → MYSTRAN Converter + Runner     ║")
    println("╚══════════════════════════════════════════════════╝")
    println()
    println("  Input:    $input_name")
    println("  Output:   $(basename(mystran_bdf))")
    println("  Results:  $(basename(results_dir))/")
    println("  MYSTRAN:  $(basename(mystran_exe))")
    println()

    # Step 1: Convert
    println("─── Step 1: Converting NASTRAN → MYSTRAN ───")
    content = convert_nastran_to_mystran(input_bdf)
    write(mystran_bdf, content)
    println("  Written: $(basename(mystran_bdf))")
    println()

    # Step 2: Run MYSTRAN (from the directory containing the .bdf)
    println("─── Step 2: Running MYSTRAN ───")
    original_dir = pwd()
    cd(input_dir)
    try
        success = run_mystran(mystran_bdf, mystran_exe)
        if success
            println("  MYSTRAN process exited OK.")
        else
            println("  WARNING: MYSTRAN exited with non-zero code.")
        end
    catch e
        println("  WARNING: MYSTRAN execution error: $e")
        println("  Attempting to collect any partial results...")
    finally
        cd(original_dir)
    end
    println()

    # Step 3: Collect results
    println("─── Step 3: Collecting results ───")
    collect_results(mystran_bdf, results_dir)
    println()

    # Step 4: Validate F06 and generate JFEM binary
    println("─── Step 4: Validating results & generating JFEM binary ───")
    mystran_bdf_base = replace(basename(mystran_bdf), r"\.[^.]+$" => "")
    f06_path = joinpath(results_dir, mystran_bdf_base * ".F06")
    jfem_path = joinpath(results_dir, model_name * ".MYSTRAN.jfem")

    if !isfile(f06_path)
        println("  ERROR: F06 file not found: $f06_path")
        println("  MYSTRAN run failed — no output generated.")
        println()
        println("FAILED.")
        return
    end

    # Validate F06 for successful completion
    f06_ok, f06_errors = validate_f06(f06_path)

    if !isempty(f06_errors)
        println("  MYSTRAN errors found:")
        for err in f06_errors
            println("    ✗ $err")
        end
    end

    if !f06_ok
        println()
        println("  ERROR: MYSTRAN did not complete successfully (no MYSTRAN END found).")
        println("  The analysis failed — .jfem binary NOT generated.")
        println("  Check the F06 file for details: $(basename(f06_path))")
        println()
        println("FAILED.")
        return
    end

    println("  ✓ MYSTRAN completed successfully (MYSTRAN END found)")

    # Parse BDF geometry (includes v3 constraint/load data)
    println("  Parsing BDF geometry: $(basename(input_bdf))")
    geom = parse_bdf_geometry(input_bdf)

    # Build v3 constraint tables
    jfem_celas, jfem_rbe2s, jfem_rbe3s = build_constraint_tables(geom)

    # Parse MYSTRAN F06 results
    println("  Parsing MYSTRAN F06:  $(basename(f06_path))")
    subcases_data = parse_mystran_f06(f06_path, geom)

    # Augment each subcase with SPC/force/moment data from BDF
    augmented = []
    for sc in subcases_data
        sid = sc.sid
        sc_info = get(geom.subcase_info, sid, (spc_id=0, load_id=0))
        spc_data = collect_spc_for_subcase(geom, sc_info.spc_id)
        force_data, moment_data = collect_forces_for_subcase(geom, sc_info.load_id)
        push!(augmented, (sid=sc.sid, disp=sc.disp, shell=sc.shell, bar=sc.bar, rod=sc.rod,
                          spc=spc_data, forces=force_data, moments=moment_data))
    end

    # Write v3 JFEM binary
    println("  Writing JFEM v3 binary: $(basename(jfem_path))")
    write_jfem_binary(jfem_path, geom, augmented;
                      jfem_celas=jfem_celas, jfem_rbe2s=jfem_rbe2s, jfem_rbe3s=jfem_rbe3s)
    println()

    println("Done.")
end

main()
