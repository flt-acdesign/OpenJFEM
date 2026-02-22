using HDF5
using WriteVTK
using LinearAlgebra

"""
    convert_h5_full()
    
Reads Nastran HDF5 and extracts Displacements, SPC Forces, Stresses, Strains, and Bar Forces.
"""
function convert_h5_full()
    
    src_dir = @__DIR__
    h5_filename = normpath(joinpath(src_dir, "..", "..", "models", "ala2.h5"))
    output_filename = replace(h5_filename, ".h5" => "_full_results")

    if !isfile(h5_filename)
        println("ERROR: File not found: $h5_filename")
        return
    end

    println(">>> Opening HDF5 file: $h5_filename")
    
    h5open(h5_filename, "r") do file
        
        
        
        println("--- Reading Nodes ---")
        if !haskey(file, "NASTRAN/INPUT/NODE/GRID")
            println("ERROR: No GRID data found.")
            return
        end
        grids = read(file["NASTRAN/INPUT/NODE/GRID"])
        
        
        node_ids = [g.ID for g in grids]
        coords_vec = [g.X for g in grids] 
        coords = reduce(hcat, coords_vec) 
        
        
        node_map = Dict(id => i for (i, id) in enumerate(node_ids))
        num_nodes = length(node_ids)

        
        
        
        println("--- Reading Connectivity ---")
        cells = MeshCell[]
        
        
        
        elem_map = Dict{Int, Int}() 
        current_cell_idx = 0

        function process_elements(path, vtk_type)
            if haskey(file, path)
                data = read(file[path])
                added_count = 0
                for row in data
                    try
                        nids_nas = row.G 
                        if all(id -> haskey(node_map, id), nids_nas)
                            nids_vtk = [node_map[id] for id in nids_nas]
                            push!(cells, MeshCell(vtk_type, nids_vtk))
                            
                            current_cell_idx += 1
                            elem_map[row.EID] = current_cell_idx
                            added_count += 1
                        end
                    catch; end
                end
                println("    Loaded $added_count elements from $path")
            end
        end

        
        process_elements("NASTRAN/INPUT/ELEMENT/CQUAD4", VTKCellTypes.VTK_QUAD)
        process_elements("NASTRAN/INPUT/ELEMENT/CTRIA3", VTKCellTypes.VTK_TRIANGLE)
        process_elements("NASTRAN/INPUT/ELEMENT/CBAR",   VTKCellTypes.VTK_LINE)
        

        num_cells = length(cells)

        
        
        
        println("--- Reading Nodal Results ---")
        
        
        disp_vec = zeros(3, num_nodes)
        spc_vec  = zeros(3, num_nodes)

        
        function read_nodal_vector(path, target_array)
            if haskey(file, path)
                data = read(file[path])
                if !isempty(data)
                    
                    target_domain = data[1].DOMAIN_ID
                    count = 0
                    for row in data
                        if row.DOMAIN_ID == target_domain && haskey(node_map, row.ID)
                            idx = node_map[row.ID]
                            target_array[1, idx] = row.X
                            target_array[2, idx] = row.Y
                            target_array[3, idx] = row.Z
                            count += 1
                        end
                    end
                    println("    Loaded $(count) entries from $path (Subcase $target_domain)")
                end
            end
        end

        read_nodal_vector("NASTRAN/RESULT/NODAL/DISPLACEMENT", disp_vec)
        read_nodal_vector("NASTRAN/RESULT/NODAL/SPC_FORCE", spc_vec)

        
        
        
        println("--- Reading Element Results ---")

        
        
        stress_vm = zeros(num_cells)
        strain_vm = zeros(num_cells)
        bar_axial_force = zeros(num_cells)

        
        
        function read_element_scalar(base_path, result_array, field_candidates)
            
            
            
            for sub_table in ["CQUAD4", "CTRIA3", "CBAR", "CQUADR", "CTRIAR"]
                full_path = "$base_path/$sub_table"
                if haskey(file, full_path)
                    data = read(file[full_path])
                    if isempty(data); continue; end
                    
                    
                    valid_field = nothing
                    row_sample = data[1]
                    for f in field_candidates
                        if hasproperty(row_sample, f)
                            valid_field = f
                            break
                        end
                    end

                    if !isnothing(valid_field)
                        
                        target_domain = data[1].DOMAIN_ID
                        for row in data
                            
                            if row.DOMAIN_ID == target_domain && haskey(elem_map, row.EID)
                                idx = elem_map[row.EID]
                                result_array[idx] = getproperty(row, valid_field)
                            end
                        end
                    end
                end
            end
        end

        
        # Nastran H5 usually uses "VON_MISES" or "VONM"
        read_element_scalar("NASTRAN/RESULT/ELEMENT/STRESS", stress_vm, [:VON_MISES, :VONM, :VM])

        
        read_element_scalar("NASTRAN/RESULT/ELEMENT/STRAIN", strain_vm, [:VON_MISES, :VONM, :VM, :MAX_SHEAR])

        
        # For CBAR, the axial force field is typically "AXIAL"
        
        read_element_scalar("NASTRAN/RESULT/ELEMENT/FORCE", bar_axial_force, [:AXIAL, :FORCE])

        
        
        
        println(">>> Writing VTK: $output_filename.vtu")
        vtk = vtk_grid(output_filename, coords, cells)
        
        
        vtk["Node_ID", VTKPointData()] = node_ids
        vtk["Displacement", VTKPointData()] = disp_vec
        vtk["SPC_Reaction_Forces", VTKPointData()] = spc_vec
        
        
        # Note: Elements that don't have a value (e.g., bars generally don't have Von Mises in the same way shells do)
        
        vtk["Stress_VonMises", VTKCellData()] = stress_vm
        vtk["Strain_VonMises", VTKCellData()] = strain_vm
        vtk["Bar_Axial_Force", VTKCellData()] = bar_axial_force
        
        vtk_save(vtk)
        println(">>> Done! Loaded Displacements, SPCs, Stress, Strain, and Forces.")
    end
end


    convert_h5_full()
