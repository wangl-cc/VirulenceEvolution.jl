module VirulenceEvolution

export gillespie_single, gillespie_meta

include("sampling.jl")

include("utility.jl")

include("singlepopulation.jl")

include("metapopulation.jl")

end # module
