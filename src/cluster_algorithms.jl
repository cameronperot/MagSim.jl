"""
	wolff!(model::Ising, avg_cluster_size::Int)

Implementation of the Wolff single cluster algorithm for the Ising model.

Arguments
* `model`           : Ising type
* `avg_cluster_size`: Average cluster size
"""
function wolff!(model::Ising, avg_cluster_size::Int)
	P          = 1 - exp(-2 * model.params.β)
	n_clusters = Int(ceil(model.params.L^2 / avg_cluster_size))
	t₀         = floor(Int, model.params.cutoff * model.params.n_sweeps)
	indices    = 1:length(model.σ)

	for t in 1:model.params.n_sweeps
		for n in 1:n_clusters
			idx           = CartesianIndices(model.σ)[rand(model.rng, indices)]
			i, j          = idx[1], idx[2]
			old_spin      = model.σ[i, j]
			new_spin      = -old_spin
			model.σ[i, j] = new_spin

			cluster = [(i, j)]
			for (i, j) in cluster
				for (k, l) in get_neighbor_indices(model, (i, j))
					if model.σ[k, l] == old_spin && rand(model.rng) < P
						model.σ[k, l] = new_spin
						push!(cluster, (k, l))
					end
				end
			end

			model.observables.avg_cluster_size += length(cluster)
		end

		t > t₀ && update_observables!(model)
	end

	model.observables.avg_cluster_size /= (model.params.n_sweeps * n_clusters)
	compute_observables_statistics!(model)
end


"""
	wolff!(model::Potts, avg_cluster_size::Int)

Implementation of the Wolff single cluster algorithm for the Potts model.
Arguments
* `model`           : Potts type
* `avg_cluster_size`: Average cluster size
"""
function wolff!(model::Potts, avg_cluster_size::Int)
	P          = 1 - exp(-model.params.β)
	n_clusters = Int(ceil(model.params.L^2 / avg_cluster_size))
	t₀         = floor(Int, model.params.cutoff * model.params.n_sweeps)
	indices    = 1:length(model.σ)

	for t in 1:model.params.n_sweeps
		for n in 1:n_clusters
			idx           = CartesianIndices(model.σ)[rand(model.rng, indices)]
			i, j          = idx[1], idx[2]
			old_spin      = model.σ[i, j]
			new_spin      = rand(model.rng, UnitRange{Int8}(1, model.params.q))
			model.σ[i, j] = new_spin

			cluster = [(i, j)]
			for (i, j) in cluster
				for (k, l) in get_neighbor_indices(model, (i, j))
					if model.σ[k, l] == old_spin && rand(model.rng) < P
						model.σ[k, l] = new_spin
						push!(cluster, (k, l))
					end
				end
			end

			model.observables.avg_cluster_size += length(cluster)
		end

		t > t₀ && update_observables!(model)
	end

	model.observables.avg_cluster_size /= (model.params.n_sweeps * n_clusters)
	compute_observables_statistics!(model)
end


"""
	wolff!(model::XY, avg_cluster_size::Int)

Implementation of the Wolff single cluster algorithm for the XY model.

Arguments
* `model`           : XY type
* `avg_cluster_size`: Average cluster size
"""
function wolff!(model::XY, avg_cluster_size::Int)
	P          = 1 - exp(-2 * model.params.β)
	n_clusters = Int(ceil(model.params.L^2 / avg_cluster_size))
	t₀         = floor(Int, model.params.cutoff * model.params.n_sweeps)
	indices    = 1:length(model.σ)

	for t in 1:model.params.n_sweeps
		for n in 1:n_clusters
			r       = random_XYVector(model.rng)
			idx     = CartesianIndices(model.σ)[rand(model.rng, indices)]
			i, j    = idx[1], idx[2]
			flip_spin!(model, (i, j), r)

			stack   = [(i, j)]
			cluster = Set(stack)
			while length(stack) > 0
				i, j = pop!(stack)
				spin = model.σ[i, j]

				for (k, l) in get_neighbor_indices(model, (i, j))
					if (k, l) ∉ cluster
						neighbor_spin = model.σ[k, l]

						if rand(model.rng) < 1 - exp(2 * model.params.β * dot(spin, r) * dot(neighbor_spin, r))
							push!(stack, (k, l))
							push!(cluster, (k, l))
							flip_spin!(model, (k, l), r)
						end
					end
				end
			end

			t > t₀ && (model.observables.avg_cluster_size += length(cluster))
		end

		t > t₀ && update_observables!(model)
	end

	model.observables.avg_cluster_size /= (model.params.n_sweeps * n_clusters)
	compute_observables_statistics!(model)
end


"""
	swendsen_wang!(model::Ising)

Implementation of the Swendsen-Wang multiple cluster algorithm for the Ising model.

Arguments
* `model`: Ising type
"""
function swendsen_wang!(model::Ising)
	P          = 1 - exp(-2 * model.params.β)
	t₀         = floor(Int, model.params.cutoff * model.params.n_sweeps)

	for t in 1:model.params.n_sweeps
		clustered = Set{Tuple{Int, Int}}()

		for j in 1:model.params.L, i in model.params.L
			(i, j) ∈ clustered && continue
			old_spin = model.σ[i, j]
			new_spin = rand(model.rng, Int8[-1, 1])
			model.σ[i, j] = new_spin
			push!(clustered, (i, j))

			stack = [(i, j)]
			while length(stack) > 0
				(i, j) = pop!(stack)
				for (k, l) in get_neighbor_indices(model, (i, j))
					if model.σ[k, l] == old_spin && (k, l) ∉ clustered && rand(model.rng) < P
						model.σ[k, l] = new_spin
						push!(stack, (k, l))
						push!(clustered, (k, l))
					end
				end
			end
		end

		t > t₀ && update_observables!(model)
	end

	compute_observables_statistics!(model)
end


"""
	swendsen_wang!(model::Potts)

"""
function swendsen_wang!(model::Potts)
end


"""
	swendsen_wang!(model::XY)

"""
function swendsen_wang!(model::XY)
end
