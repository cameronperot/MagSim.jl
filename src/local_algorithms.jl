"""
	metropolis!(model::Ising)

Implementation of the Metropolis local update algorithm for the Ising model.
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

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	metropolis!(model::Potts)

Implementation of the Metropolis local update algorithm for the Potts model.
"""
function metropolis!(model::Potts)
	t₀ = floor(Int, model.params.cutoff * model.params.n_sweeps)
	eᵝ = exp(-model.params.β)

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			counts   = compute_counts(model, i, j)
			old_spin = model.σ[i, j]
			new_spin = rand(model.rng, UnitRange{Int8}(1, model.params.q))
			ΔE       = counts[old_spin] - counts[new_spin]

			if ΔE <= 0 || rand(model.rng) < eᵝ^ΔE
				model.σ[i, j] = new_spin
				model.observables.E_current += ΔE
			end
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	metropolis!(model::XY)

Implementation of the Metropolis local update algorithm for the XY model.
"""
function metropolis!(model::XY)
	t₀ = floor(Int, model.params.cutoff * model.params.n_sweeps)
	eᵝ = exp(-model.params.β)

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			old_spin = model.σ[i, j]
			new_spin = random_XYVector(model.rng)
			ΔE       = dot(old_spin .- new_spin, sum_of_neighbors(model, i, j))

			if ΔE <= 0 || rand(model.rng) < eᵝ^ΔE
				model.σ[i, j] = new_spin
				model.observables.E_current  += ΔE
				model.observables.Mx_current += new_spin[1] - old_spin[1]
				model.observables.My_current += new_spin[2] - old_spin[2]
			end
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	heat_bath!(model::Ising)

Implementation of the heat bath algorithm for the Ising model.
"""
function heat_bath!(model::Ising)
	t₀   = floor(Int, model.params.cutoff * model.params.n_sweeps)
	ΔEs  = [-8, -4, 0, 4, 8]
	exps = exp(-model.params.β) .^ ΔEs

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			k = Int(model.σ[i, j] * sum_of_neighbors(model, i, j) / 2 + 3)

			if rand(model.rng) < exps[k] / (exps[k] + 1)
				model.σ[i, j] *= -1
				model.observables.E_current += ΔEs[k]
				model.observables.M_current += 2 * model.σ[i, j]
			end
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	heat_bath!(model::Potts)

Implementation of the heat bath algorithm for the Potts model.
"""
function heat_bath!(model::Potts)
	t₀   = floor(Int, model.params.cutoff * model.params.n_sweeps)
	exps = exp.(model.params.β .* collect(0:4))

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			counts   = compute_counts(model, i, j)
			old_spin = model.σ[i, j]
			ps       = exp.(model.params.β .* counts); ps ./= sum(ps)
			R        = rand(model.rng)
			P        = 0.0

			for (new_spin, p) in enumerate(ps)
				P += p
				if R < P
					model.σ[i, j] = new_spin
					model.observables.E_current += counts[old_spin] - counts[new_spin]
					break
				end
			end
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	heat_bath!(model::XY)

Implementation of the Metropolis local update algorithm for the XY model.
"""
function heat_bath!(model::XY)
	t₀ = floor(Int, model.params.cutoff * model.params.n_sweeps)
	eᵝ = exp(-model.params.β)

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			old_spin = model.σ[i, j]
			new_spin = random_XYVector(model.rng)
			ΔE       = dot(old_spin .- new_spin, sum_of_neighbors(model, i, j))

			if rand(model.rng) < eᵝ^ΔE / (eᵝ^ΔE + 1)
				model.σ[i, j] = new_spin
				model.observables.E_current  += ΔE
				model.observables.Mx_current += new_spin[1] - old_spin[1]
				model.observables.My_current += new_spin[2] - old_spin[2]
			end
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
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
	compute_counts(model::Potts, i::Int, j::Int)

Computes the counts/frequency of the different spin vaues neighboring `σ[i, j]`.
"""
function compute_counts(model::Potts, i::Int, j::Int)
	counts = zeros(Int8, model.params.q)

	counts[model.σ[plus(i, model.params.L), j]]  += 1
	counts[model.σ[minus(i, model.params.L), j]] += 1
	counts[model.σ[i, plus(j, model.params.L)]]  += 1
	counts[model.σ[i, minus(j, model.params.L)]] += 1

	return counts
end


"""
	sum_of_neighbors(model::XY, i::Int, j::Int)

Sums the values of the neighbors for a given site (i, j).
"""
function sum_of_neighbors(model::XY, i::Int, j::Int)
	return (
		model.σ[plus(i, model.params.L), j] .+
		model.σ[minus(i, model.params.L), j] .+
		model.σ[i, plus(j, model.params.L)] .+
		model.σ[i, minus(j, model.params.L)]
		)
end
