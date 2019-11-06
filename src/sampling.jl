# sampling functions for simulation
# Copy and modified from pkg StatsBase
# Link: https://github.com/JuliaStats/StatsBase.jl/blob/master/src/sampling.jl#L413
# [MIT License](https://github.com/JuliaStats/StatsBase.jl/blob/master/LICENSE.md)

using Random

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
    linear2cart(idx, collect(size(wa)))
end

sample(wa) = sample(Random.GLOBAL_RNG, wa)

function linear2cart(idx::Integer, sizes::AbstractVector{<:Integer})
    a = accumulate(*, sizes[1:end-1])
    d = similar(a)
    m = similar(a)
    d[end], m[end] = fldmod1(idx, a[end])
    @inbounds for i = length(a) - 1:-1:1
        (d[i], m[i]) = fldmod1(m[i+1], a[i])
    end
    m[1], d...
end
