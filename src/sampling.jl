# sampling functions for simulation
# Copy and modified from pkg StatsBase
# Link: https://github.com/JuliaStats/StatsBase.jl/blob/master/src/sampling.jl
# [MIT License](https://github.com/JuliaStats/StatsBase.jl/blob/master/LICENSE.md)

using Random

"""
    sample([rng], wv::AbstractArray)

Select a single random integer in `wa` with probabilities
proportional to the weights given in `wa`.


Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).

"""
function sample(rng::AbstractRNG, wv::AbstractVector)
    t = rand(rng) * sum(wv)
    n = length(wv)
    i = 1
    cw = wv[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += wv[i]
    end
    i
end

function sample(rng::AbstractRNG, wa::AbstractArray)
    idx = sample(rng, vec(wa))
    ind2sub(axes(wa), idx)
end

sample(wa::AbstractArray) = sample(Random.GLOBAL_RNG, wa)

const ind2sub = Base._ind2sub
