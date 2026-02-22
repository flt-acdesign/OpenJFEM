# main.jl — JFEM solver entry point and orchestrator


using LinearAlgebra
using SparseArrays
using Printf
using WriteVTK
using Statistics
using JSON

println(">>> Loading Modules...")

include("FEMKernels.jl"); using .FEM
include("parsing/NastranParser.jl"); using .NastranParser
include("solver/Solver.jl"); using .Solver
include("ModelBuilder.jl")
include("Export.jl")

function main(filename::String; output_dir::Union{String,Nothing}=nothing)
    if !isfile(filename)
        println("ERROR: File not found at: $filename")
        return
    end

    script_dir = @__DIR__
    if isnothing(output_dir)
        output_dir = joinpath(script_dir, "..", "output")
    end

    if !isdir(output_dir)
        mkdir(output_dir)
    end

    # --- Parse BDF ---
    println(">>> Reading BDF file: $filename")
    lines = readlines(filename)

    println(">>> Checking format...")
    lines = NastranParser.convert_mystran_to_nastran(lines)

    println(">>> Parsing Bulk Data...")
    cc, bulk = NastranParser.read_bulk_and_case(lines)
    cards = NastranParser.process_cards(bulk)

    # --- Card inventory ---
    export_card_inventory(cards, output_dir, filename)

    # --- Build model ---
    model = build_model(cards, cc)
    transform_geometry!(model)

    # --- Assemble global stiffness ---
    K, id_map, X, ndof, node_R, max_elem_stiff, rbe3_map, snorm_normals = Solver.assemble_stiffness(model)

    # --- Prepare results accumulators ---
    global_results = Dict(
        "displacements" => [],
        "spc_forces" => [],
        "forces" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => [], "conrod" => [], "celas1" => []),
        "stresses" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => [], "conrod" => [], "celas1" => []),
        "strains" => Dict("cbar" => [], "quad4" => [], "tria3" => [], "crod" => [], "conrod" => [], "celas1" => [])
    )

    sorted_sids = sort(collect(keys(cc["SUBCASES"])))

    # --- Build JFEM binary element tables ---
    jfem_node_ids, jfem_quads, jfem_trias, jfem_bars, jfem_rods = build_jfem_element_tables(model, id_map)
    jfem_celas, jfem_rbe2s, jfem_rbe3s = build_jfem_constraint_tables(model, id_map)
    jfem_subcases_data = []

    # --- Solve each subcase ---
    for sid in sorted_sids
        sub = cc["SUBCASES"][sid]
        println("\n>>> Solving Subcase $sid...")
        load_id = get(sub, "LOAD", nothing)
        spc_id = get(sub, "SPC", nothing)
        u, stresses, sub_res = Solver.solve_case(K, ndof, model, id_map, X, load_id, spc_id, node_R; max_elem_stiff=max_elem_stiff, rbe3_map=rbe3_map, snorm_normals=snorm_normals)

        # Accumulate results
        append!(global_results["displacements"], sub_res["displacements"])
        append!(global_results["spc_forces"], sub_res["spc_forces"])
        for k in keys(sub_res["forces"]); append!(global_results["forces"][k], sub_res["forces"][k]); end
        for k in keys(sub_res["stresses"]); append!(global_results["stresses"][k], sub_res["stresses"][k]); end
        for k in keys(sub_res["strains"]); append!(global_results["strains"][k], sub_res["strains"][k]); end

        # Collect JFEM subcase data (v3: includes SPC, forces, moments)
        sc_data = collect_jfem_subcase_data(u, sub_res, id_map, jfem_node_ids, jfem_quads, jfem_trias, jfem_bars, jfem_rods; model=model, spc_id=spc_id, load_id=load_id)
        push!(jfem_subcases_data, (sid=sid, disp=sc_data.disp, shell=sc_data.shell, bar=sc_data.bar, rod=sc_data.rod, spc=sc_data.spc, forces=sc_data.forces, moments=sc_data.moments))

        # Export VTK for this subcase
        export_vtk_subcase(filename, output_dir, sid, model, id_map, X, u, stresses)
    end

    # --- Export aggregated results ---
    export_json(filename, output_dir, global_results)
    export_jfem_binary(filename, output_dir, id_map, X, jfem_node_ids, jfem_quads, jfem_trias, jfem_bars, jfem_rods, jfem_subcases_data; jfem_celas=jfem_celas, jfem_rbe2s=jfem_rbe2s, jfem_rbe3s=jfem_rbe3s)
end

if !isempty(ARGS)
    target_file = ARGS[1]
    if !isabspath(target_file)
        target_file = joinpath(@__DIR__, "..", target_file)
    end
    target_file = normpath(target_file)
    out_dir = length(ARGS) >= 2 ? normpath(ARGS[2]) : nothing
    main(target_file; output_dir=out_dir)
else
    target_file = joinpath(@__DIR__, "..", "models", "OpenJFEM.bdf")
    target_file = normpath(target_file)
    main(target_file)
end
