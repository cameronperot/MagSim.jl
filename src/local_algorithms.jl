"""
	metropolis!(I::Ising_2D)

Implementation of the famous Metropolis local update algorithm.
"""
function metropolis!(I::Ising_2D)
	t₀   = floor(Int, I.params.cutoff * I.params.n_sweeps)
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

		t > t₀ && update_observables!(I, t)
	end

	compute_observables_statistics!(I)
	return I
end


"""
	metropolis!(model::Potts_2D)

"""
function metropolis!(model::Potts_2D)
end


"""
	metropolis!(model::XY_2D)

"""
function metropolis!(model::XY_2D)
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
	sum_of_neighbors(I::Ising_2D, i::Int, j::Int)

Sums the values of the neighbors for a given site (i, j).
"""
function sum_of_neighbors(I::Ising_2D, i::Int, j::Int)
	return (
		I.σ[plus(i, I.params.L), j] +
		I.σ[minus(i, I.params.L), j] +
		I.σ[i, plus(j, I.params.L)] +
		I.σ[i, minus(j, I.params.L)]
		)
end
