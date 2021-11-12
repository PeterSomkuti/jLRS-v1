module jLRS

# standard library
using Printf

# third party libraries
using Distributions
using HDF5
using Proj4
using Glob
using AstroLib
using DataFrames
using Dates
using LinearAlgebra
using Random
using SatelliteToolbox
using SQLite
using SparseArrays
using StaticArrays


jLRS_seed = 123456
# Seed the generator at the beginning. So the *same* code with the *same* inputs
# will always produce the same numbers. If even one parameter changes, we will
# get different results!
Random.seed!(jLRS_seed)

# Type definitions required by the rest of the module
include("types.jl")

# Common math-y functions
include("misc.jl")

# Solar data required for many calculations
include("solar.jl")

# Surface samplings (BRDF, land cover, etc.)
include("surface.jl")

# SIF-related functions
include("SIF.jl")

# Functions for OCO-type
include("OCO.jl")

# Functions for Geostationary-type
include("Geostationary.jl")

# Functions for SwathNadir-type
include("SwathNadirOrbiter.jl")

# Main calculations
include("light_response.jl")

# Plotting tools
include("plotting.jl")

@info "You have loaded jLRS - thanks!"

end # module
