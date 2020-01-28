"""
	AbstractMagnetismModel

Abstract type

Subtypes
* Ising
* Potts
* XY
"""
abstract type AbstractMagnetismModel end


"""
	Parameters

"""
mutable struct Parameters
	L         ::Int
	β         ::Float64
	n_sweeps  ::Int
	cutoff    ::Float64
	start_type::String
	seed      ::Int
end


"""
	Observables()

"""
mutable struct Observables
	E_current        ::Int
	E                ::Array
	M_current        ::Int
	M                ::Array
	avg_cluster_size ::Float64
	statistics       ::Dict

	function Observables()
		new(0, [], 0, [], 0., Dict())
	end
end


"""
	initialize_σ_2D(L::Int, start_type::String, rng::MersenneTwister)

"""
function initialize_σ_2D(L::Int, start_type::String, rng::MersenneTwister)
	if start_type == "cold"
		σ = ones(Int, (L, L))
	elseif start_type == "hot"
		σ = rand(rng, Int[-1, 1], (L, L))
	end

	return σ
end


"""
	compute_E_2D(σ::Array{Int, 2}, L::Int)

"""
function compute_E_2D(σ::Array{Int, 2}, L::Int)
	E₀ = 0
	for j in 1:L
		for i in 1:L
			E₀ -= σ[i, j] * (σ[plus(i, L), j] + σ[i, plus(j, L)])
		end
	end

	return E₀
end


"""
	Ising_2D(
		L         ::Int,
		β         ::Float64;
		J         ::Float64=1.,
		n_sweeps  ::Int    =10^4,
		start_type::String ="cold",
		seed      ::Int    =8,
	)

"""
mutable struct Ising_2D <: AbstractMagnetismModel
	σ          ::Array{Int, 2}
	params     ::Parameters
	observables::Observables
	rng        ::MersenneTwister

	function Ising_2D(
		L         ::Int,
		β         ::Float64;
		n_sweeps  ::Int    =10^4,
		cutoff    ::Float64=0.2,
		start_type::String ="cold",
		seed      ::Int    =8,
	)

		params = Parameters(L, β, n_sweeps, cutoff, start_type, seed)
		rng = MersenneTwister(seed)
		σ = initialize_σ_2D(L, start_type, rng)
		observables = Observables()
		observables.E_current = compute_E_2D(σ, L)
		observables.M_current = sum(σ)

		new(σ, params, observables, rng)
	end
end


mutable struct Potts_2D <: AbstractMagnetismModel
end


mutable struct XY_2D <: AbstractMagnetismModel
end
