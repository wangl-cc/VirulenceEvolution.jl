using LinearAlgebra

function gillespie(mutfunc, # virulence mutation function
                   T::Real; # max time
                   # variable parameters
                   S::AbstractVector{IntType}, # susceptible host population sizes
                   I::AbstractMatrix{IntType}, # infected host population sizes
                   R::AbstractMatrix{IntType}, # infected host population sizes
                   v::AbstractVector{RealType}, # virulences
                   # fixed parameters
                   β::RealType, # infect coefficient
                   b::RealType, # birth rate
                   d::RealType, # death rate
                   r::RealType, # recovery rate
                   μ::RealType, # immunity loss rate
                   mt::RealType, # matation rate
                   mg::RealType, # immigration rate
                   maxepoch::Integer = 10_000_000) where {IntType <: Integer,RealType <: Real}
    length(S) == size(I, 2) == size(R, 2) || error("S, I, R must have the same metapopulation size!")
    length(v) == size(I, 1) == size(R, 1) || error("S, I, R must have the same number of viruses species!")

    t =  zero(RealType) # initial time
    epoch = 0x0000000000000000 # init epoch

    @apply copy S I R v

    dynamics_S = Dynamics.(t, S)
    dynamics_I = Dynamics.(t, sum(I, dims=1))
    dynamics_R = Dynamics.(t, sum(R, dims=1))

    maxepoch = UInt(maxepoch)
    while t <= T && epoch <= maxepoch && length(I) > 0
        epoch += 0x0000000000000001

        # calculate reactions

        population_num = length(S)
        viruses_num = length(v)

        reproduction_S = b * S
        reproduction_I = b * I
        reproduction_R = b * R

        death_S = d * S
        death_I = d * I
        death_R = d * R

        immigration_S = mg * S
        immigration_I = mg * I
        immigration_R = mg * R

        infect_S = Array{RealType}(undef, viruses_num, population_num)
        infect_R = Array{RealType}(undef, viruses_num, viruses_num, population_num)
        for i in 1:population_num
            Ii = slicecol(I, i, viruses_num)
            Ri = slicecol(R, i, viruses_num)
            infect_S[:, i] = β * S[i] * Ii
            infect_R[:, :, i] = zeroatdiag!(β * Ii * transpose(Ri))
        end

        recovery = r * I
        immuneloss = μ * R

        kill = v .* I

        pathogenmutation = mt * I

        τ, r_idx, idx = choosereaction(reproduction_S, reproduction_I, reproduction_R,
                                       death_S, death_I, death_R,
                                       immigration_S, immigration_I, immigration_R,
                                       infect_S, infect_R, recovery, immuneloss,
                                       kill, pathogenmutation)

        t += τ # time goes

        if r_idx == 1 # S reproduction
            S[idx] += 1
        elseif r_idx == 2 # I reproduction
            S[idx[2]] += 1
        elseif r_idx == 3 # R reproduction
            S[idx[2]] += 1
        elseif r_idx == 4 # S death
            S[idx] -= 1
        elseif r_idx == 5 # I death
            i, j = idx
            I[i, j] -= 1
            if minimum(I[i, :] .<= 0) && minimum(R[i, :] .<= 0)
                I = I[1:end .!= i, :]
                R = R[1:end .!= i, :]
                deleteat!(v, i)
            end
        elseif r_idx == 6 # R death
            i, j = idx
            R[i, j] -= 1
            if minimum(I[i, :] .<= 0) && minimum(R[i, :] .<= 0)
                I = I[1:end .!= i, :]
                R = R[1:end .!= i, :]
                deleteat!(v, i)
            end
        elseif r_idx == 7 # S immigration
            S[idx] -= 1
            S[rand(1:population_num)] += 1
        elseif r_idx == 8 # I immigration
            I[idx] -= 1
            I[idx[1], rand(1:population_num)] += 1
        elseif r_idx == 9 # R immigration
            R[idx] -= 1
            R[idx[1], rand(1:population_num)] += 1
        elseif r_idx == 10 # I_i infect S
            S[idx[2]] -= 1
            I[idx] += 1
        elseif r_idx == 11 # I_i infect R_j
            i, j, k = idx
            I[i, k] += 1

            R[j, k] -= 1
            if minimum(I[j, :] .<= 0) && minimum(R[j, :] .<= 0)
                I = I[1:end .!= j, :]
                R = R[1:end .!= j, :]
                deleteat!(v, j)
            end
        elseif r_idx == 12 # recovery
            I[idx] -= 1
            R[idx] += 1
        elseif r_idx == 13 # immunity loss
            i, j = idx
            S[j] += 1
            R[i, j] -= 1
            if minimum(I[i, :] .<= 0) && minimum(R[i, :] .<= 0)
                I = I[1:end .!= i, :]
                R = R[1:end .!= i, :]
                deleteat!(v, i)
            end
        elseif r_idx == 14 # pathogen kill host
            i, j = idx
            I[i, j] -= 1
            if minimum(I[i, :] .<= 0) && minimum(R[i, :] .<= 0)
                I = I[1:end .!= i, :]
                R = R[1:end .!= i, :]
                deleteat!(v, i)
            end
        elseif r_idx == 15 # pathogen mutation
            i, j = idx
            I = vcat(I, repeat([zero(IntType)], outer = (1, population_num)))
            R = vcat(R, repeat([zero(IntType)], outer = (1, population_num)))
            I[i, j] -= 1
            I[end, j] += 1
            new_v = mutfunc(v[i])
            push!(v, new_v)

            if minimum(I[i, :] .<= 0) && minimum(R[i, :] .<= 0)
                I = I[1:end .!= i, :]
                R = R[1:end .!= i, :]
                deleteat!(v, i)
            end
        end
        record!.(dynamics_S, t, S)
        record!.(dynamics_I, t, sum(I, dims=1))
        record!.(dynamics_R, t, sum(R, dims=1))
    end
    println("Simulation ends at epoch ", Int(epoch), "!")
    dynamics_S, dynamics_I, dynamics_R
end
