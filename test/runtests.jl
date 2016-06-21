if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end
using Compat
using Dopri

include("fortran.jl")
include("lowlevel.jl")
include("highlevel.jl")
