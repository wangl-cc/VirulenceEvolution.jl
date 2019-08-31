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
struct Dynamics{R<:Real, T}
    history::Vector{Tuple{R, T}}
end

Dynamics(t::Real, n) = Dynamics([(t, n)])

function record!(d::Dynamics, t::Real, n)
    push!(d.history, (t, n))
end
