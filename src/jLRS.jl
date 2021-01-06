module jLRS

using Printf
using HDF5
using Proj4
using Glob


# Type definitions required by the rest of the module
include("types.jl")
# Common math-y functions (may use jLRS types)
include("misc.jl")


# Functions for OCO-type
include("OCO.jl")

print("You have loaded jLRS")


end # module
