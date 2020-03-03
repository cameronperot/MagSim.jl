"""
	compute_E_Ising(σ::Array{Int8, 2})

Computes the Ising model internal energy for the input spin configuration `σ`.
"""
function compute_E_Ising(σ::Array{Int8, 2})
	E::Int = 0
	L = size(σ, 1)
	for j in 1:L, i in 1:L
		E -= σ[i, j] * (σ[plus(i, L), j] + σ[i, plus(j, L)])
	end

	return E
end


"""
	compute_E_Potts(σ::Array{Int8, 2})

Computes the Potts model internal energy for the input spin configuration `σ`.
"""
function compute_E_Potts(σ::Array{Int8, 2})
	E::Int = 0
	L = size(σ, 1)
	for j in 1:L, i in 1:L
		E -= Int(σ[i, j] == σ[plus(i, L), j]) + Int(σ[i, j] == σ[i, plus(j, L)])
	end

	return E
end


"""
	compute_E_M_XY(σ::Array{NTuple{2, Float64}, 2})

Computes the XY model internal energy for the input spin configuration `σ`.
"""
function compute_E_M_XY(σ::Array{NTuple{2, Float64}, 2})
	E = 0.
	Mx = 0.
	My = 0.
	L = size(σ, 1)
	for j in 1:L, i in 1:L
		E  -= dot(σ[i, j], σ[plus(i, L), j]) + dot(σ[i, j], σ[i, plus(j, L)])
		Mx += σ[i, j][1]
		My += σ[i, j][2]
	end

	return (E, Mx, My)
end


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
