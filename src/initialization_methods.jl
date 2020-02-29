"""
	compute_E_Ising(σ::Array{Int, 2})

Computes the Ising model internal energy for the input spin configuration `σ`.
"""
function compute_E_Ising(σ::Array{Int8, 2})
	E₀::Int = 0
	L = size(σ, 1)
	for j in 1:L
		for i in 1:L
			E₀ -= σ[i, j] * (σ[plus(i, L), j] + σ[i, plus(j, L)])
		end
	end

	return E₀
end


"""
	compute_E_Potts(σ::Array{Int, 2})

Computes the Potts model internal energy for the input spin configuration `σ`.
"""
function compute_E_Potts(σ::Array{Int8, 2})
	E₀::Int = 0
	L = size(σ, 1)
	for j in 1:L
		for i in 1:L
			E₀ -= Int(σ[i, j] == σ[plus(i, L), j]) + Int(σ[i, j] == σ[i, plus(j, L)])
		end
	end

	return E₀
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
	initialize_σ_Ising(q::Int, L::Int, start_type::Symbol, rng::MersenneTwister)


Initialize the Potts model lattice `σ` with values ∈ {1, ..., q}
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
