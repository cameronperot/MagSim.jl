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
	Ising(
		L         ::Int,
		β         ::Real;
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

	function Ising(
		L         ::Int,
		β         ::Real;
		n_sweeps  ::Int    =10^5,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

		params                = Parameters(L, β, n_sweeps, cutoff, start_type, seed, 0)
		rng                   = MersenneTwister(params.seed)
		σ                     = initialize_σ_Ising(params.L, params.start_type, rng)
		observables           = Observables()
		observables.E_current = compute_E_Ising(σ)
		observables.M_current = sum(σ)

		return new(σ, params, observables, rng)
	end
end


"""
	Potts(
		q         ::Int,
		L         ::Int,
		β         ::Real;
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
* `σ`          : Lattice containing the spins, values are ∈ {1, ..., q}
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

	function Potts(
		q         ::Int,
		L         ::Int,
		β         ::Real;
		n_sweeps  ::Int    =10^5,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

		params                = Parameters(L, β, n_sweeps, cutoff, start_type, seed, q)
		rng                   = MersenneTwister(seed)
		σ                     = initialize_σ_Potts(q, L, start_type, rng)
		observables           = Observables()
		observables.E_current = compute_E_Potts(σ)

		return new(σ, params, observables, rng)
	end
end


"""
	XY(
		L         ::Int,
		β         ::Real;
		n_sweeps  ::Int    =10^5,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

Type representing an XY model.

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
* A new instance of `XY`

Attributes
* `σ`          : Lattice containing the spins, values are 2D unit vectors
* `params`     : Parameters type, contains relevant parameters for the simulation
* `observables`: Observables type, contains relevant observables data for the simulation
* `rng`        : MersenneTwister type random number generator

# Examples
```julia-repl
julia> XY(32, 1)
XY
L          = 32
β          = 1.0
n_sweeps   = 100000
start_type = cold
seed       = 8
q          = 0
```
"""
struct XY <: AbstractMagnetismModel
	σ          ::Array{NTuple{2, Float64}, 2}
	params     ::Parameters
	observables::XYObservables
	rng        ::MersenneTwister

	function XY(
		L         ::Int,
		β         ::Real;
		n_sweeps  ::Int    =10^5,
		cutoff    ::Float64=0.2,
		start_type::Symbol =:cold,
		seed      ::Int    =8,
		)

		params      = Parameters(L, β, n_sweeps, cutoff, start_type, seed, 0)
		rng         = MersenneTwister(seed)
		σ           = initialize_σ_XY(L, start_type, rng)
		observables = XYObservables()
		observables.E_current, observables.Mx_current, observables.Mx_current = compute_E_M_XY(σ)

		return new(σ, params, observables, rng)
	end
end
