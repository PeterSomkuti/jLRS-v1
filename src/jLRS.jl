module jLRS

# standard library
using Printf

# third party libraries
using HDF5
using Proj4
using Glob
using AstroLib
using DataFrames
using SQLite

# Type definitions required by the rest of the module
include("types.jl")

# Common math-y functions (may use jLRS types)
include("misc.jl")

# Functions for OCO-type
include("OCO.jl")

# Main calculations
include("light_response.jl")

# Plotting tools
include("plotting.jl")

@info "You have loaded jLRS - thanks!"

end # module
