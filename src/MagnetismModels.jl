"""
	AbstractMagnetismModel

Abstract type for various magnetism models.

Subtypes
* `Ising`
* `Potts`
* `XY`
"""
abstract type AbstractMagnetismModel end


function show(io::IO, model::AbstractMagnetismModel)
	println(typeof(model))
	show(model.params)
end


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
	q         ::Int
end


function show(io::IO, params::Parameters)
	println(
		"""
		L          = $(params.L)
		β          = $(params.β)
		n_sweeps   = $(params.n_sweeps)
		start_type = $(params.start_type)
		seed       = $(params.seed)
		q          = $(params.q)
		"""
		)
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
	Ising(
		L         ::Int,
		β         ::Float64;
		n_sweeps  ::Int    =10^5,
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
* A new instance of `Ising`

Attributes
* `σ`          : Lattice containing the spins, values are ±1
* `params`     : Parameters type, contains relevant parameters for the simulation
* `observables`: Observables type, contains relevant observables data for the simulation
* `rng`        : MersenneTwister type random number generator

# Examples
```julia-repl
julia> Ising(32, log(1 + √2) / 2)
Ising
L          = 32
β          = 0.3
n_sweeps   = 100000
start_type = cold
seed       = 8
q          = 0
```
"""
struct Ising <: AbstractMagnetismModel
	σ          ::Array{Int8, 2}
	params     ::Parameters
	observables::Observables
	rng        ::MersenneTwister

	function Ising(params::Parameters)
		rng                   = MersenneTwister(params.seed)
		σ                     = initialize_σ_Ising(params.L, params.start_type, rng)
		observables           = Observables()
		observables.E_current = compute_E_Ising(σ)
		observables.M_current = sum(σ)

		new(σ, params, observables, rng)
	end

	function Ising(
		L         ::Int,
		β         ::Real;
		n_sweeps  ::Int    =10^5,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

		params = Parameters(L, β, n_sweeps, cutoff, start_type, seed, 0)

		Ising(params)
	end
end


"""
	Potts(
		q         ::Int,
		L         ::Int,
		β         ::Float64;
		n_sweeps  ::Int    =10^5,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

Type representing a 2D Potts q-state model.

Arguments
* `q`: The number of different states each spin can take on
* `L`: The side length of the lattice, i.e. for a square lattice there will be `L^2` spins
* `β`: `1 / kT`, i.e. the inverse temperature

Keyword Arguments
* `n_sweeps`  : The number of update sweeps to run over the lattice
* `cutoff`    : The percentage of measurements to throw away, i.e. throw away values during
	thermalization
* `start_type`: `:cold` starts the lattice with all spins aligned, `:hot` all spins random
* `seed`      : Seed value for the random number generator

Returns
* A new instance of `Potts`

Attributes
* `σ`          : Lattice containing the spins, values are ±1
* `params`     : Parameters type, contains relevant parameters for the simulation
* `observables`: Observables type, contains relevant observables data for the simulation
* `rng`        : MersenneTwister type random number generator

# Examples
```julia-repl
julia> Potts(8, 32, log(1 + √8))
Potts
L          = 32
β          = 0.3
n_sweeps   = 100000
start_type = cold
seed       = 8
q          = 8
```
"""
struct Potts <: AbstractMagnetismModel
	σ          ::Array{Int8, 2}
	params     ::Parameters
	observables::Observables
	rng        ::MersenneTwister

	function Potts(params::Parameters)
		rng                   = MersenneTwister(params.seed)
		σ                     = initialize_σ_Potts(params.q, params.L, params.start_type, rng)
		observables           = Observables()
		observables.E_current = compute_E_Potts(σ)

		new(σ, params, observables, rng)
	end

	function Potts(
		q         ::Int,
		L         ::Int,
		β         ::Real;
		n_sweeps  ::Int    =10^5,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

		params = Parameters(L, β, n_sweeps, cutoff, start_type, seed, q)

		Potts(params)
	end
end


"""
	XY

"""
struct XY <: AbstractMagnetismModel
end
