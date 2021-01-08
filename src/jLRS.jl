module jLRS

# standard library
using Printf

# third party libraries
using HDF5
using Proj4
using Glob
using AstroLib


# Type definitions required by the rest of the module
include("types.jl")
# Common math-y functions (may use jLRS types)
include("misc.jl")


# Functions for OCO-type
include("OCO.jl")

# Main calculations
include("light_response.jl")

print("You have loaded jLRS - thanks!")


end # module
