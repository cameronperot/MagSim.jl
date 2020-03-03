"""
	metropolis!(model::Ising)

Implementation of the Metropolis algorithm for the Ising model.
"""
function metropolis!(model::Ising)
	t₀   = floor(Int, model.params.cutoff * model.params.n_sweeps)
	ΔEs  = [-8, -4, 0, 4, 8]
	exps = exp(-model.params.β) .^ ΔEs

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			k = Int(model.σ[i, j] * sum_of_neighbors(model, i, j) / 2 + 3)

			(ΔEs[k] <= 0 || rand(model.rng) < exps[k]) && (model.σ[i, j] *= -1)
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	metropolis!(model::Potts)

Implementation of the Metropolis algorithm for the Potts model.
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

			(ΔE <= 0 || rand(model.rng) < eᵝ^ΔE) && (model.σ[i, j] = new_spin)
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	metropolis!(model::XY)

Implementation of the Metropolis algorithm for the XY model.
"""
function metropolis!(model::XY)
	t₀ = floor(Int, model.params.cutoff * model.params.n_sweeps)
	eᵝ = exp(-model.params.β)

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			old_spin = model.σ[i, j]
			new_spin = random_XYVector(model.rng)
			ΔE       = dot(old_spin .- new_spin, sum_of_neighbors(model, i, j))

			(ΔE <= 0 || rand(model.rng) < eᵝ^ΔE) && (model.σ[i, j] = new_spin)
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

			(rand(model.rng) < exps[k] / (exps[k] + 1)) && (model.σ[i, j] *= -1)
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
				R < P && (model.σ[i, j] = new_spin; break)
			end
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end


"""
	heat_bath!(model::XY)

Implementation of the heat bath algorithm for the XY model.
"""
function heat_bath!(model::XY)
	t₀ = floor(Int, model.params.cutoff * model.params.n_sweeps)
	eᵝ = exp(-model.params.β)

	for t in 1:model.params.n_sweeps
		for j in 1:model.params.L, i in 1:model.params.L
			old_spin = model.σ[i, j]
			new_spin = random_XYVector(model.rng)
			ΔE       = dot(old_spin .- new_spin, sum_of_neighbors(model, i, j))

			(rand(model.rng) < eᵝ^ΔE / (eᵝ^ΔE + 1)) && (model.σ[i, j] = new_spin)
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
	return model
end
