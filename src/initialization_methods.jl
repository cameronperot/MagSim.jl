"""
	initialize_σ_Ising(L::Int, start_type::Symbol, rng::MersenneTwister)

Initialize the Ising model lattice `σ` with values ∈ {-1, 1}.
"""
function initialize_σ_Ising(L::Int, start_type::Symbol, rng::MersenneTwister)
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
	initialize_σ_Potts(q::Int, L::Int, start_type::Symbol, rng::MersenneTwister)


Initialize the Potts model lattice `σ` with values ∈ {1, ..., q}.
"""
function initialize_σ_Potts(q::Int, L::Int, start_type::Symbol, rng::MersenneTwister)
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
	initialize_σ_XY(L::Int, start_type::Symbol, rng::MersenneTwister)


Initialize the XY model lattice `σ` with unit vectors represented by `NTuple{2, Float64}`.
"""
function initialize_σ_XY(L::Int, start_type::Symbol, rng::MersenneTwister)
	if start_type == :cold
		σ = [(0.0, 1.0) for i in 1:L, j in 1:L]
	elseif start_type == :hot
		σ = [random_XYVector(rng) for i in 1:L, j in 1:L]
	else
		error("invalid start_type provided, valid choices are [:cold, :hot]")
	end

	return σ
end


"""
	random_XYVector(rng::MersenneTwister)

Generates a random unit vector represented by `NTuple{2, Float64}`.
"""
function random_XYVector(rng::MersenneTwister)
	θ = rand(rng) * 2π
	return (cos(θ), sin(θ))
end
