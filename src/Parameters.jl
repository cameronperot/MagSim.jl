"""
    Parameters

Type containing the parameters of a magnetism simulation.

Attributes
* `L`         : The side length of the lattice, e.g. for a square lattice there will be
    `L^2` spins
* `β`         : `1 / kT`, i.e. the inverse temperature
* `n_sweeps`  : The number of update sweeps to run over the lattice
* `cutoff`    : The first percentage of observable measurements to discard, useful to get
    rid of values during thermalization period
* `start_type`: `:cold` starts the lattice with all spins aligned, `:hot` all spins random
* `seed`      : Seed value for the random number generator
"""
struct Parameters
    L::Int
    β::Float64
    n_sweeps::Int
    cutoff::Float64
    start_type::Symbol
    seed::Int
    q::Int
end


function show(io::IO, params::Parameters)
    println("""
            L          = $(params.L)
            β          = $(params.β)
            n_sweeps   = $(params.n_sweeps)
            start_type = $(params.start_type)
            seed       = $(params.seed)
            q          = $(params.q)
            """)
end
