"""
	AbstractMagnetismModel

Abstract type for various magnetism models.

Subtypes
* `Ising_2D`
* `Potts_2D`
* `XY_2D`
"""
abstract type AbstractMagnetismModel end


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
	L         ::Int
	β         ::Float64
	n_sweeps  ::Int
	cutoff    ::Float64
	start_type::Symbol
	seed      ::Int
end


"""
	Observables()

Type containing the observables data for a magnetism simulation.

Attributes
* `E_current`       : The current energy of the system
* `E`               : Array tracking the energy of the system over the updates
* `M_current`       : The current magnetization of the system
* `M`               : Array tracking the magnetization of the system over the updates
* `avg_cluster_size`: Average cluster size (applies only to cluster algorithms)
* `statistics`      : Dictionary containing the per spin values of average energy, average
	magnetization, heat capacity, magnetic susceptibility, Binder parameter, and the
	integrated autocorrelation time (computed using the binning method)
"""
mutable struct Observables
	E_current        ::Int
	E                ::Array{Int, 1}
	M_current        ::Int
	M                ::Array{Int, 1}
	avg_cluster_size ::Float64
	statistics       ::Dict

	function Observables()
		new(0, [], 0, [], 0, Dict())
	end
end


"""
	Ising_2D(
		L         ::Int,
		β         ::Float64;
		n_sweeps  ::Int    =10^4,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

Type representing a 2D Ising model.

Arguments
* `L`: The side length of the lattice, i.e. for a square lattice there will be `L^2` spins
* `β`: `1 / kT`, i.e. the inverse temperature

Keyword Arguments
* `n_sweeps`  : The number of update sweeps to run over the lattice
* `cutoff`    : The percentage of measurements to throw away, i.e. throw away values during
	thermalization
* `start_type`: `:cold` starts the lattice with all spins aligned, `:hot` all spins random
* `seed`      : Seed value for the random number generator

Returns
* A new instance of `Ising_2D`

Attributes
* `σ`          : Lattice containing the spins, values are ±1
* `params`     : Parameters type, contains relevant parameters for the simulation
* `observables`: Observables type, contains relevant observables data for the simulation
* `rng`        : MersenneTwister type random number generator

# Examples
```julia-repl
julia> Ising_2D(32, log(1 + √2) / 2)
Ising_2D
L          = 32
β          = 0.3
n_sweeps   = 100000
start_type = cold
seed       = 8
```
"""
struct Ising_2D <: AbstractMagnetismModel
	σ          ::Array{Int8, 2}
	params     ::Parameters
	observables::Observables
	rng        ::MersenneTwister

	function Ising_2D(params::Parameters)
		rng                   = MersenneTwister(params.seed)
		σ                     = initialize_σ_Ising_2D(params.L, params.start_type, rng)
		observables           = Observables()
		observables.E_current = compute_E_Ising_2D(σ, params.L)
		observables.M_current = sum(σ)

		new(σ, params, observables, rng)
	end

	function Ising_2D(
		L         ::Int,
		β         ::Float64;
		n_sweeps  ::Int    =10^4,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

		params = Parameters(L, β, n_sweeps, cutoff, start_type, seed)

		Ising_2D(params)
	end
end


function show(io::IO, model::AbstractMagnetismModel)
	println(typeof(model))
	show(model.params)
end


function show(io::IO, params::Parameters)
	println("""
			L          = $(params.L)
			β          = $(params.β)
			n_sweeps   = $(params.n_sweeps)
			start_type = $(params.start_type)
			seed       = $(params.seed)
			"""
			)
end


"""
	Potts_2D

"""
struct Potts_2D <: AbstractMagnetismModel
end


"""
	XY_2D

"""
struct XY_2D <: AbstractMagnetismModel
end
