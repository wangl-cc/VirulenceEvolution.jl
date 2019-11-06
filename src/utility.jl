"""
    choosereaction(reactionrates::AbstractArray{R}...)where R<:Real

Calculate the reaction time `τ` and randomly choose an reaction from reactions with weights `reactions`.
"""
function choosereaction(reactions::AbstractArray{<:Real}...)
    sum_rs = [sum(i) for i in reactions]
    sum_all = sum(sum_rs)
    τ = -log(rand()) / sum_all
    r_idx = sample(sum_rs) # reaction index
    idx = sample(reactions[r_idx]) # index in reaction
    return τ, r_idx, idx
end

"""
    Dynamics{R<:Real, I}

A struct records changes over time.
"""
struct Dynamics{R <: Real,T}
    history::Vector{Tuple{R,T}}
end

Dynamics(t::Real, n) = Dynamics([(t, n)])

function record!(d::Dynamics, t::Real, n)
    push!(d.history, (t, copy(n)))
    d
end

"""
    zeroatdiag!(m::AbstractMatrix{T}) where T

Make diagonal elements of given matrix `m` zero.

# Example
```jldoctest
julia> m = rand(2, 2)
2×2 Array{Float64,2}:
 0.301516  0.101581
 0.216389  0.672865

julia> zeroatdiag!(m)
2×2 Array{Float64,2}:
 0.0       0.101581
 0.216389  0.0
```
"""
function zeroatdiag!(m::AbstractMatrix{T})where T
    zero_t = zero(T)
    for i in 1:size(m, 1)
        m[i, i] = zero_t
    end
    m
end

"""
    apply(f, vars::Symbol...)

Apply the `f` to  all `vars` in place.

# Example
```jldoctest
julia> @macroexpand @apply copy a b c
quote
    a = copy(a)
    b = copy(b)
    c = copy(c)
end

julia> a = 1
1

julia> @apply x->x+1 a
2

julia> a
2
```
"""
macro apply(f, vars::Symbol...)
    ex = Expr(:block)
    for var in vars
        push!(ex.args, :($var = $f($var)))
    end
    esc(ex)
end

@inline function slicecol(m::AbstractMatrix, i::Integer, colsize::Integer)
    m[(i-1)*colsize+1:i*colsize]
end

@inline slicecol(m::AbstractMatrix, i::Integer) = slicecol(m, i, size(m, 1))

Base.getindex(a::Array, t::Tuple) = Core.arrayref(true, a, t...)
Base.setindex!(a::Array{T}, x::T, t::Tuple) where {T} = Core.arrayset(true, a, x, t...)
Base.setindex!(a::Array{T}, x, t::Tuple) where {T} = Core.arrayset(true, a, convert(T, x)::T, t...)
