function gillespie(S::Integer, # susceptible host population size
                   I::AbstractVector, # infected host population sizes
                   R::AbstractVector, # reoverd host population sizes
                   T::Real; # max time
                   # parameters
                   β::Real, # infect coefficient
                   b::Real, # birth rate
                   d::Real, # death rate
                   v::AbstractArray, # virulence
                   r::Real, #recovery rate
                   μ::Real, # immunity loss rate
                   m::Real, # matation rate
                   mutfunc, # virulence mutation function
                   maxepoch::Integer=10000000)
    length(R) == length(I) == length(v) || error("R, I and v must have same length!")
    t = 0.0 # initial time
    epoch = 0x0000000000000000 # init epoch

    dynamics_S = Dynamics(t, S)

    I = copy(I)
    dynamics_I_current = [Dynamics(t, i) for i in I]
    dynamics_I_all = copy(dynamics_I_current)

    R = copy(R)
    dynamics_R_current = [Dynamics(t, i) for i in R]
    dynamics_R_all = copy(dynamics_R_current)

    v = copy(v)
    dynamics_v = Dynamics(t, (copy(I), copy(v)))
    v_all = copy(v)

    maxepoch = UInt(maxepoch)
    while t <= T && epoch <= maxepoch
        epoch += 0x0000000000000001
        length(I) <= 0 && break # check whether or not pathogen extinct

        # calculate reactions

        reproduction_S = b * S
        reproduction_I = b * I
        reproduction_R = b * R

        death_S = d * S
        death_I = d * I
        death_R = d * R

        infect_S = β * S * I # I_i infect S
        infect_R = β * I * transpose(R) # I_i infect R_j
        for i in 1:size(infect_R, 1)
            infect_R[i, i] = 0
        end
        recovery = r * I
        immuneloss = μ * R

        kill = v .* I

        pathogenmutation = m * I

        τ, r_idx, idx = choosereaction([reproduction_S, death_S], reproduction_I, reproduction_R,
                                       death_I, death_R, infect_S,
                                       infect_R, recovery, immuneloss,
                                       kill, pathogenmutation)

        t += τ # time goes

        if r_idx == 1 # S reproduction and death
            if idx::Int == 1
                S += 1
            elseif idx == 2
                S -= 1
            else
                error("Unexpected idx at reaction $r_idx.")
            end
            record!(dynamics_S, t, S)
        elseif r_idx == 2 || r_idx == 3 # I and R reproduction
            S += 1
            record!(dynamics_S, t, S)
        elseif r_idx == 4 # I death
            n = I[idx] -=1
            record!(dynamics_I_current[idx], t, n)
            record!(dynamics_v, t, (copy(I), copy(v)))
            if n <= 0 && R[idx] <= 0
                deleteat!(I, idx)
                deleteat!(R, idx)
                deleteat!(dynamics_I_current, idx)
                deleteat!(dynamics_R_current, idx)
                deleteat!(v, idx)
            end
        elseif r_idx == 5 # R death
            n = R[idx] -= 1
            record!(dynamics_R_current[idx], t, n)
            if n <= 0 && I[idx] <= 0
                deleteat!(I, idx)
                deleteat!(R, idx)
                deleteat!(dynamics_I_current, idx)
                deleteat!(dynamics_R_current, idx)
                deleteat!(v, idx)
            end
        elseif r_idx == 6 # I_i infect S
            S -= 1
            record!(dynamics_S, t, S)

            n = I[idx] += 1
            record!(dynamics_I_current[idx], t, n)
            record!(dynamics_v, t, (copy(I), copy(v)))
        elseif r_idx == 7 # I_i infect R_j
            i, j = idx
            n = I[i] += 1
            record!(dynamics_I_current[i], t, n)
            record!(dynamics_v, t, (copy(I), copy(v)))

            n = R[j] -= 1
            record!(dynamics_R_current[j], t, n)
            if n <= 0 && I[j] <= 0
                deleteat!(I, j)
                deleteat!(R, j)
                deleteat!(dynamics_I_current, j)
                deleteat!(dynamics_R_current, j)
                deleteat!(v, j)
            end
        elseif r_idx == 8 # recovery
            n = I[idx] -=1
            record!(dynamics_I_current[idx], t, n)
            record!(dynamics_v, t, (copy(I), copy(v)))

            n = R[idx] += 1
            record!(dynamics_R_current[idx], t, n)
        elseif r_idx == 9 # immunity loss
            S += 1
            record!(dynamics_S, t, S)

            n = R[idx] -= 1
            record!(dynamics_R_current[idx], t, n)
            if n <= 0 && I[idx] <= 0
                deleteat!(I, idx)
                deleteat!(R, idx)
                deleteat!(dynamics_I_current, idx)
                deleteat!(dynamics_R_current, idx)
                deleteat!(v, idx)
            end
        elseif r_idx == 10 # pathogen kill host
            n = I[idx] -= 1
            record!(dynamics_I_current[idx], t, n)
            record!(dynamics_v, t, (copy(I), copy(v)))
            if n <= 0 && R[idx] <= 0
                deleteat!(I, idx)
                deleteat!(R, idx)
                deleteat!(dynamics_I_current, idx)
                deleteat!(dynamics_R_current, idx)
                deleteat!(v, idx)
            end
        elseif r_idx == 11 # pathogen mutation
            n = I[idx] -= 1
            record!(dynamics_I_current[idx],t , n)

            new_dynamics_I = Dynamics(t, 1)
            new_dynamics_R = Dynamics(t, 0)
            push!(I, 1)
            push!(R, 0)
            push!(dynamics_I_all, new_dynamics_I)
            push!(dynamics_R_all, new_dynamics_R)
            push!(dynamics_I_current, new_dynamics_I)
            push!(dynamics_R_current, new_dynamics_R)
            new_v = mutfunc(v[idx])
            push!(v, new_v)
            push!(v_all, new_v)

            record!(dynamics_v, t, (copy(I), copy(v)))

            if n <= 0 && R[idx] <= 0
                deleteat!(I, idx)
                deleteat!(R, idx)
                deleteat!(dynamics_I_current, idx)
                deleteat!(dynamics_R_current, idx)
                deleteat!(v, idx)
            end
        end
    end
    println("Simulation ends at epoch ", Int(epoch), "!")
    return dynamics_S, dynamics_I_all, dynamics_R_all, v_all, dynamics_v
end
