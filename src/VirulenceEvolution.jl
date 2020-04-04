module VirulenceEvolution

export Dynamics, record!, gillespie_single, gillespie_meta

include("sampling.jl")

include("utilities.jl")

include("singlepopulation.jl")

include("metapopulation.jl")

end # module
