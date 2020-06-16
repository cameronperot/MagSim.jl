"""
    initialize_σ_Ising(L::Int, start_type::Symbol, rng::AbstractRNG)

Initialize the Ising model lattice `σ` with values ∈ {-1, 1}.

Arguments
* `L`         : Side length of the square lattice
* `start_type`: Either `:hot` or `:cold` to initialize the lattice with a hot or cold start
* `rng`       : Random number generator
Returns
* `σ`         : 2D array representing a square lattice
"""
function initialize_σ_Ising(L::Int, start_type::Symbol, rng::AbstractRNG)
    if start_type == :cold
        σ = ones(Int8, (L, L))
    elseif start_type == :hot
        σ = rand(rng, Int8[-1, 1], (L, L))
    else
        error("invalid start_type provided, valid choices are [:cold, :hot]")
    end

    return σ
end


"""
    initialize_σ_Potts(q::Int, L::Int, start_type::Symbol, rng::AbstractRNG)


Initialize the Potts model lattice `σ` with values ∈ {1, ..., q}.

Arguments
* `q`: Number of possible spin values
* `L`         : Side length of the square lattice
* `start_type`: Either `:hot` or `:cold` to initialize the lattice with a hot or cold start
* `rng`       : Random number generator
Returns
* `σ`         : 2D array representing a square lattice
"""
function initialize_σ_Potts(q::Int, L::Int, start_type::Symbol, rng::AbstractRNG)
    if start_type == :cold
        σ = ones(Int8, (L, L))
    elseif start_type == :hot
        σ = rand(rng, UnitRange{Int8}(1, q), (L, L))
    else
        error("invalid start_type provided, valid choices are [:cold, :hot]")
    end

    return σ
end


"""
    initialize_σ_XY(L::Int, start_type::Symbol, rng::AbstractRNG)


Initialize the XY model lattice `σ` with unit vectors represented by `NTuple{2, Float64}`.

Arguments
* `L`         : Side length of the square lattice
* `start_type`: Either `:hot` or `:cold` to initialize the lattice with a hot or cold start
* `rng`       : Random number generator
Returns
* `σ`         : 2D array representing a square lattice
    represented by a two-tuple
"""
function initialize_σ_XY(L::Int, start_type::Symbol, rng::AbstractRNG)
    if start_type == :cold
        σ = [(0.0, 1.0) for i = 1:L, j = 1:L]
    elseif start_type == :hot
        σ = [random_XYVector(rng) for i = 1:L, j = 1:L]
    else
        error("invalid start_type provided, valid choices are [:cold, :hot]")
    end

    return σ
end


"""
    random_XYVector(rng::AbstractRNG)

Generates a random unit vector represented by `NTuple{2, Float64}`.

Arguments
* `rng` : Random number generator
Returns
* A two-tuple representing a unit-vector with random angular orientation
"""
function random_XYVector(rng::AbstractRNG)
    θ = rand(rng) * 2π
    return (cos(θ), sin(θ))
end
