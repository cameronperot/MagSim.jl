"""
	compute_E_Ising_2D(σ::Array{Int, 2}, L::Int)

"""
function compute_E_Ising_2D(σ::Array{Int8, 2}, L::Int)
	E₀::Int = 0
	for j in 1:L
		for i in 1:L
			E₀ -= σ[i, j] * (σ[plus(i, L), j] + σ[i, plus(j, L)])
		end
	end

	return E₀
end


"""
	initialize_σ_Ising_2D(L::Int, start_type::Symbol, rng::MersenneTwister)

"""
function initialize_σ_Ising_2D(L::Int, start_type::Symbol, rng::MersenneTwister)
	if start_type == :cold
		σ = ones(Int8, (L, L))
	elseif start_type == :hot
		σ = rand(rng, Int8[-1, 1], (L, L))
	else
		error("invalid start_type provided, valid choices are [:cold, :hot]")
	end

	return σ
end
