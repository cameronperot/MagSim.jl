"""
    simulate_Ising_parallel(
        algorithm ::Function,
        Ls        ::Array{Int, 1},
        βs        ::Union{Array{Float64}, AbstractRange{Float64}};
        n_sweeps  ::Int     =10^5,
        cutoff    ::Float64 =0.2,
        start_type::Symbol  =:cold,
        seed      ::Int     =8,
        )

Runs multiple Ising model simulations in parallel.

Arguments
* `algorithm`  : The algorithm with which to run the simulation
* `Ls`         : The side lengths for which to run the simulations with
* `βs`         : The inverse temperatures for which to run the simulations with
* `n_sweeps`   : The number of update sweeps for each simulation
* `cutoff`     : How much of the initial data to throw away (i.e. to ignore the observations
    when the system was thermalizing
* `start_type` : Whether or not to start the simulation with a cold or hot start, available
    options are `[:cold, :hot]`
* `seed`       : Random number generator seed value
Returns
* A dictionary with all of the simulations, where the keys are the specific `(L, β)` values
"""
function simulate_Ising_parallel(
    algorithm::Function,
    Ls::Array{Int,1},
    βs::Union{Array{Float64},AbstractRange{Float64}};
    n_sweeps::Int = 10^5,
    cutoff::Float64 = 0.2,
    start_type::Symbol = :cold,
    seed::Int = 8,
)

    out = Dict{Tuple,Ising}()
    Threads.@threads for β in βs
        Threads.@threads for L in Ls
            out[(L, β)] = algorithm(Ising(
                L,
                β,
                n_sweeps = n_sweeps,
                cutoff = cutoff,
                start_type = start_type,
                seed = seed,
            ))
        end
    end

    return out
end


"""
    simulate_Potts_parallel(
        algorithm ::Function,
        Ls        ::Array{Int, 1},
        βs        ::Union{Array{Float64}, AbstractRange{Float64}};
        n_sweeps  ::Int     =10^5,
        cutoff    ::Float64 =0.2,
        start_type::Symbol  =:cold,
        seed      ::Int     =8,
        )

Runs multiple Potts model simulations in parallel.

Arguments
* `q`          : The number of possible states
* `algorithm`  : The algorithm with which to run the simulation
* `Ls`         : The side lengths for which to run the simulations with
* `βs`         : The inverse temperatures for which to run the simulations with
* `n_sweeps`   : The number of update sweeps for each simulation
* `cutoff`     : How much of the initial data to throw away (i.e. to ignore the observations
    when the system was thermalizing
* `start_type` : Whether or not to start the simulation with a cold or hot start, available
    options are `[:cold, :hot]`
* `seed`       : Random number generator seed value
Returns
* A dictionary with all of the simulations, where the keys are the specific `(L, β)` values
"""
function simulate_Potts_parallel(
    q::Int,
    algorithm::Function,
    Ls::Array{Int,1},
    βs::Union{Array{Float64},AbstractRange{Float64}};
    n_sweeps::Int = 10^5,
    cutoff::Float64 = 0.2,
    start_type::Symbol = :cold,
    seed::Int = 8,
)

    out = Dict{Tuple,Potts}()
    Threads.@threads for β in βs
        Threads.@threads for L in Ls
            out[(L, β)] = algorithm(Potts(
                q,
                L,
                β,
                n_sweeps = n_sweeps,
                cutoff = cutoff,
                start_type = start_type,
                seed = seed,
            ))
        end
    end

    return out
end
