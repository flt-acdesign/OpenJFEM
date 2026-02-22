# Solver.jl — Module hub for the JFEM linear static solver
# Includes: helpers, SNORM normals, load resolution, constraint assembly,
#           stiffness assembly, boundary conditions, stress recovery, and
#           the top-level solve_case orchestrator.

module Solver

using LinearAlgebra
using SparseArrays
using Statistics
using IterativeSolvers
using AlgebraicMultigrid
using Dates
using StaticArrays
using ..FEM

include("helpers.jl")
include("snorm.jl")
include("loads.jl")
include("constraints.jl")
include("assembly.jl")
include("boundary_conditions.jl")
include("stress_recovery.jl")
include("solve_case.jl")

end # module Solver
