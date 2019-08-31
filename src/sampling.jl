# sampling functions for simulation
# Copy and modified from pkg StatsBase
# Link: https://github.com/JuliaStats/StatsBase.jl/blob/master/src/sampling.jl#L413
# [MIT License](https://github.com/JuliaStats/StatsBase.jl/blob/master/LICENSE.md)

using Random

function sample(rng::AbstractRNG, wv::AbstractVector{R})where R <: Real
    t = rand(rng) * sum(wv)
    n = length(wv)
    i = 1
    cw = wv[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += wv[i]
    end
    return i
end

function sample(rng::AbstractRNG, wa::AbstractArray{R})where R <: Real
    h = size(wa, 1)
    i = sample(rng, vec(wa))
    d = fld(i-1, h)
    i-(d*h), d+1
end

sample(wv) = sample(Random.GLOBAL_RNG, wv)
