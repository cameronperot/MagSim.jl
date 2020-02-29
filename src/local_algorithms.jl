"""
	metropolis!(model::Ising)

Implementation of the famous Metropolis local update algorithm for the Ising model.
"""
function metropolis!(model::Ising)
	t₀   = floor(Int, model.params.cutoff * model.params.n_sweeps)
	ΔEs  = [-8, -4, 0, 4, 8]
	exps = exp(-model.params.β) .^ ΔEs

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			k = Int(model.σ[i, j] * sum_of_neighbors(model, i, j) / 2 + 3)

			if ΔEs[k] <= 0 || rand(model.rng) < exps[k]
				model.σ[i, j] *= -1
				model.observables.E_current += ΔEs[k]
				model.observables.M_current += 2 * model.σ[i, j]
			end
		end

		t > t₀ && update_observables!(model, t)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	metropolis!(model::Potts)

Implementation of the famous Metropolis local update algorithm for the Potts model.
"""
function metropolis!(model::Potts)
	t₀   = floor(Int, model.params.cutoff * model.params.n_sweeps)
	expᵦ = exp(-model.params.β)

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			counts   = compute_counts(model, i, j)
			old_spin = model.σ[i, j]
			new_spin = rand(model.rng, UnitRange{Int8}(1, model.params.q))
			ΔE       = counts[old_spin] - counts[new_spin]

			if ΔE <= 0 || rand(model.rng) < expᵦ^ΔE
				model.σ[i, j] = new_spin
				model.observables.E_current += ΔE
			end
		end

		t > t₀ && update_observables!(model, t)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	metropolis!(model::XY)

"""
function metropolis!(model::XY)
end


"""
	heat_bath!(model::Ising_2D)

"""
function heat_bath!(model::Ising_2D)
end


"""
	heat_bath!(model::Potts_2D)

"""
function heat_bath!(model::Potts_2D)
end


"""
	heat_bath!(model::XY_2D)

"""
function heat_bath!(model::XY_2D)
end


"""
	sum_of_neighbors(model::Ising, i::Int, j::Int)

Sums the values of the neighbors for a given site (i, j).
"""
function sum_of_neighbors(model::Ising, i::Int, j::Int)
	return (
		model.σ[plus(i, model.params.L), j] +
		model.σ[minus(i, model.params.L), j] +
		model.σ[i, plus(j, model.params.L)] +
		model.σ[i, minus(j, model.params.L)]
		)
end


"""
	compute_counts(q::Int, σ::Array{Int8, 1}, i::Int, j::Int)

"""
function compute_counts(model::Potts, i::Int, j::Int)
	counts = zeros(Int8, model.params.q)

	counts[model.σ[plus(i, model.params.L), j]]  += 1
	counts[model.σ[minus(i, model.params.L), j]] += 1
	counts[model.σ[i, plus(j, model.params.L)]]  += 1
	counts[model.σ[i, minus(j, model.params.L)]] += 1

	return counts
end
