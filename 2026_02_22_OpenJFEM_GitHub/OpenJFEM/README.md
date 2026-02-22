# OpenJFEM

A linear static finite element solver written in Julia with an integrated HTML5 3D post-processor.

OpenJFEM reads industry-standard Bulk Data Format (`.bdf`) input files, solves linear static problems (SOL 101), and exports results in multiple formats including its own binary JFEM format for interactive visualization.

## Features

- **Linear static analysis** (SOL 101) with direct and iterative solvers (AMG-preconditioned GMRES)
- **Multi-threaded assembly** using Julia threads for parallel element stiffness computation
- **BDF input format** compatible with standard `.bdf` files
- **MYSTRAN interop** — built-in converter to run models with the open-source [MYSTRAN](https://www.mystran.com/) solver for cross-validation
- **JFEM binary format** (v3) — compact binary results for fast post-processing
- **LUDWIG post-processor** — HTML5/Babylon.js 3D viewer with dual-model comparison, transparency, color contouring, and interactive controls
- **Comprehensive test suite** — 60+ verification models with automated regression testing

## Supported Elements

| Element | Description |
|---------|-------------|
| CQUAD4  | 4-node quadrilateral shell (membrane + bending + MITC4 transverse shear) |
| CTRIA3  | 3-node triangular shell (membrane + bending) |
| CBAR    | 2-node beam element (axial, bending, shear, torsion) with PBARL sections |
| CROD    | 2-node rod element (axial + torsion) |
| CONROD  | Rod with inline properties |
| CELAS1/2/3 | Scalar spring elements |
| CONM2   | Concentrated mass element |
| RBE2    | Rigid body element (master-slave constraints) |
| RBE3    | Interpolation element (weighted averaging) |

## Supported Properties and Materials

- **PSHELL** / **PCOMP** — isotropic and composite laminate shell properties
- **PBAR** / **PBARL** — beam cross-section properties (BOX, I, T, L, ROD, TUBE, BAR, CROSS, H, HAT, Z, CHAN, CHAN1, CHAN2)
- **PROD** — rod/tube cross-section properties
- **MAT1** — isotropic material
- **MAT2** — anisotropic 2D material
- **MAT8** — orthotropic material for shells
- **CORD2R** / **CORD2C** — rectangular and cylindrical coordinate systems

## Supported Loads and Constraints

- **FORCE** / **MOMENT** — concentrated nodal loads
- **PLOAD4** / **PLOAD2** — pressure loads on shell elements
- **PLOAD1** — distributed loads on bar elements
- **GRAV** — gravity body loads
- **LOAD** — load combination cards
- **SPC** / **SPC1** / **SPCADD** — single-point constraints
- **MPC** — general multipoint constraints
- **AUTOSPC** — automatic constraint of singular DOFs

## Installation

Requires **Julia 1.9** or later.

```bash
cd OpenJFEM
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Quick Start

### 1. Run the Solver

```bash
# Basic run
julia --project=. src/main.jl models/fwing/fwing.bdf

# Multi-threaded (recommended for large models)
julia -t auto --project=. src/main.jl models/OpenJFEM/OpenJFEM.bdf
```

The solver produces three output files next to the input:

| Extension   | Description |
|-------------|-------------|
| `.jfem`     | Binary results — load into LUDWIG post-processor |
| `.vtu`      | VTK unstructured grid — open in [ParaView](https://www.paraview.org/) |
| `.JU.JSON`  | JSON with displacements, forces, stresses, strains |
| `.CARDS.JSON` | Card inventory summary |

### 2. Visualize Results

Open `POST/postv11.html` in any modern browser (Chrome, Edge, Firefox). No server required — just double-click the file.

1. Click the **Struct** tab
2. Load a `.jfem` file into **Model A**
3. Adjust deformation scale, color field, and display options
4. Optionally load a second `.jfem` into **Model B** for side-by-side comparison
5. Use the **Vis** tab to control vertical offset, transparency, and model visibility

See `POST/POST_GUIDE.html` for detailed instructions.

### 3. Cross-validate with MYSTRAN

```bash
# Convert BDF to MYSTRAN format, run MYSTRAN, and collect results
julia --project=. tools/converters/nastran_to_mystran.jl models/fwing/fwing.bdf

# Compare OpenJFEM vs MYSTRAN results
julia --project=. tools/testing/verify_all.jl
```

See `tools/GUIDE.html` for full tool documentation.

## Project Structure

```
OpenJFEM/
├── src/                             # Core solver
│   ├── main.jl                      # Entry point and orchestrator
│   ├── FEMKernels.jl                # Element stiffness matrices (QUAD4, TRIA3, CBAR, CROD)
│   ├── ModelBuilder.jl              # Model construction and coordinate transforms
│   ├── Export.jl                    # VTK, JSON, JFEM binary export
│   ├── parsing/                     # BDF parser
│   │   ├── NastranParser.jl         # Parser module hub
│   │   ├── card_processor.jl        # Card parsing and field extraction
│   │   ├── extract_geometry.jl      # GRID, coordinate systems
│   │   ├── extract_elements.jl      # Element definitions
│   │   ├── extract_properties.jl    # Property cards (PSHELL, PBAR, etc.)
│   │   ├── extract_materials.jl     # Material cards (MAT1, MAT2, MAT8)
│   │   ├── extract_loads.jl         # Load cards (FORCE, PLOAD4, GRAV, etc.)
│   │   ├── extract_constraints.jl   # Constraint cards (SPC, RBE2, RBE3, MPC)
│   │   ├── mystran_converter.jl     # MYSTRAN format writer
│   │   └── utilities.jl             # Parsing helpers
│   └── solver/                      # Solver subsystem
│       ├── Solver.jl                # Module hub
│       ├── assembly.jl              # Global stiffness assembly (threaded)
│       ├── boundary_conditions.jl   # BC processing
│       ├── constraints.jl           # MPC/RBE constraint handling
│       ├── loads.jl                 # Load vector assembly
│       ├── solve_case.jl            # Per-subcase solver
│       ├── stress_recovery.jl       # Element stress/strain recovery
│       ├── snorm.jl                 # Shell normal computation
│       └── helpers.jl               # Solver utilities
├── tools/                           # Utility scripts
│   ├── benchmark.jl                 # Performance benchmarking
│   ├── GUIDE.html                   # Tool documentation
│   ├── comparison/                  # Result comparison tools
│   │   ├── compare_results.jl       # Full statistical comparison
│   │   ├── compare_OpenJFEM.jl        # Displacement comparison
│   │   ├── quick_compare.jl         # Quick summary comparison
│   │   └── detailed_compare.jl      # Per-node detailed comparison
│   ├── converters/                  # Format converters
│   │   ├── f06_2_json.jl            # F06 text output to JSON
│   │   ├── HDF5_2_VTK.jl           # HDF5 to VTK converter
│   │   └── nastran_to_mystran.jl    # BDF to MYSTRAN format converter
│   └── testing/                     # Test runners
│       ├── run_mystran_tests.jl     # Automated batch testing
│       ├── verify_all.jl            # Comprehensive regression suite
│       └── verify_tests.jl          # Single-test verification
├── models/                          # 60+ verification test cases
│   ├── ala3/                        # Aircraft wing FEM (574 shell nodes)
│   ├── fwing/                       # Fighter wing model
│   └── ...                          # Cantilevers, trusses, springs, RBE, composites, etc.
├── POST/                            # Post-processor
│   ├── postv11.html                 # LUDWIG — HTML5/Babylon.js 3D viewer
│   └── POST_GUIDE.html              # Viewer documentation
├── Project.toml                     # Julia project manifest
└── .gitignore
```

## Validation

The solver has been validated against MYSTRAN on models ranging from single-element tests to full assemblies with several types of card.

- **Displacements**: correlation > 0.994 across all DOFs
- **QUAD4 stresses/strains**: correlation > 0.995
- **CROD/CBAR forces**: correlation > 0.998
- **SPC reaction forces**: validated per-DOF

## JFEM Binary Format (v3)

The `.jfem` format is a compact little-endian binary format containing:

- **Header**: magic number, version, node/element counts
- **Mesh**: node coordinates, QUAD4/TRIA3/CBAR/CROD connectivity
- **Constraints**: SPC nodes, CELAS springs, RBE2/RBE3 elements
- **Subcases**: per-subcase displacements, shell forces/strains, bar forces, rod forces, applied forces/moments

The format is designed for fast loading in the LUDWIG post-processor and is fully documented in `src/Export.jl`.

## License  

GNU AFFERO GENERAL PUBLIC LICENSE (AGPL-3.0).


