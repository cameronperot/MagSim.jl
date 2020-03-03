"""
	wolff!(model::Ising)

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
			model.observables.E_current -= 2 * new_spin * sum_of_neighbors(model, i, j)
			model.observables.M_current += 2 * new_spin

			cluster = Tuple{Int, Int}[(i, j)]
			for (i, j) in cluster
				neighbors = get_neighbor_indices(model, i, j)
				for (k, l) in neighbors
					if model.σ[k, l] == old_spin && rand(model.rng) < P
						model.σ[k, l] = new_spin
						push!(cluster, (k, l))
						model.observables.E_current -= 2 * new_spin * sum_of_neighbors(model, k, l)
						model.observables.M_current += 2 * new_spin
					end
				end
			end

			model.observables.avg_cluster_size += length(cluster)
		end

		t > t₀ && update_observables!(model)
	end

	model.observables.avg_cluster_size /= (model.params.n_sweeps * n_clusters)
	compute_observables_statistics!(model)
	return model
end


"""
	wolff!(model::Potts)

"""
function wolff!(model::Potts)
end


"""
	wolff!(model::XY)

"""
function wolff!(model::XY)
end


"""
	swendsen_wang!(model::Ising)

"""
function swendsen_wang!(model::Ising)
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
