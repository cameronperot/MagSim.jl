"""
	minus(L::Int, i::Int)

	Used for finding the (periodic boundary conditions) neighbors of a node
"""
minus(i::Int, L::Int) = i ≠ 1 ? i-1 : L


"""
	plus(L::Int, i::Int)

	Used for finding the (periodic boundary conditions) neighbors of a node
"""
plus(i::Int, L::Int) = i ≠ L ? i+1 : 1


"""
	sum_of_neighbors(I::Ising_2D, i::Int, j::Int)

"""
function sum_of_neighbors(I::Ising_2D, i::Int, j::Int)
	return (
		I.σ[plus(i, I.params.L), j] +
		I.σ[minus(i, I.params.L), j] +
		I.σ[i, plus(j, I.params.L)] +
		I.σ[i, minus(j, I.params.L)]
	)
end


"""
	update_observables!(I::AbstractMagnetismModel, t::Int)

"""
function update_observables!(model::AbstractMagnetismModel, t::Int)
	push!(model.observables.E, model.observables.E_current)
	push!(model.observables.M, model.observables.M_current)
end


"""
	finalize_observables!(I::AbstractMagnetismModel, t::Int)

"""
function compute_observables_statistics!(model::AbstractMagnetismModel)
	model.observables.statistics["e"] = mean(model.observables.E ./ length(model.σ))
	model.observables.statistics["m"] = mean(abs.(model.observables.M) ./ length(model.σ))
	model.observables.statistics["c"] = model.params.β^2 * var(model.observables.E) / length(model.σ)
	model.observables.statistics["χ"] = model.params.β * var(abs.(model.observables.M)) / length(model.σ)
	model.observables.statistics["U"] = 1 - mean(model.observables.M.^4) / 3 / mean(model.observables.M.^2)^2
	model.observables.statistics["τ_int"] = compute_autocorrelation_time(model, collect(1:400))
end


"""
	metropolis!(I::Ising_2D)

"""
function metropolis!(I::Ising_2D)
	t₀   = Int(floor(I.params.cutoff * I.params.n_sweeps))
	ΔEs  = [-8, -4, 0, 4, 8]
	exps = exp(-I.params.β) .^ ΔEs

	for t in 1:I.params.n_sweeps
		for j in 1:I.params.L
			for i in 1:I.params.L
				k = Int(I.σ[i, j] * sum_of_neighbors(I, i, j) / 2 + 3)

				if ΔEs[k] <= 0 || exps[k] > rand(I.rng)
					I.σ[i, j] *= -1
					I.observables.E_current += ΔEs[k]
					I.observables.M_current += 2 * I.σ[i, j]
				end
			end
		end

		t > t₀ ? update_observables!(I, t) : continue
	end

	compute_observables_statistics!(I)
	return I
end



function metropolis!(model::Potts_2D)
end


function metropolis!(model::XY_2D)
end


function heat_bath!(model::Ising_2D)
end


function heat_bath!(model::Potts_2D)
end


function heat_bath!(model::XY_2D)
end
