"""
	get_neighbor_indices(model::Ising, i::Int, j::Int)

Returns the indices of the nearest neighbors for a given lattice site (i, j).
"""
function get_neighbor_indices(model::AbstractMagnetismModel, i::Int, j::Int)
	return (
		(plus(i, model.params.L), j),
		(minus(i, model.params.L), j),
		(i, plus(j, model.params.L)),
		(i, minus(j, model.params.L))
		)
end


"""
	sum_of_neighbors(model::Union{Ising, XY}, i::Int, j::Int)

Sums the values of the nearest neighbors for a given lattice site (i, j).
"""
function sum_of_neighbors(model::Union{Ising, XY}, i::Int, j::Int)
	return @. (
		model.σ[plus(i, model.params.L), j] +
		model.σ[minus(i, model.params.L), j] +
		model.σ[i, plus(j, model.params.L)] +
		model.σ[i, minus(j, model.params.L)]
		)
end


"""
	dot(a::NTuple{2, Float64}, b::NTuple{2, Float64})

Inner product of two 2D vectors.
"""
function dot(a::NTuple{2, Float64}, b::NTuple{2, Float64})
	return a[1] * b[1] + a[2] * b[2]
end



"""
	compute_counts(model::Potts, i::Int, j::Int)

Computes the counts/frequency of the different spin vaues neighboring `σ[i, j]`.
"""
function compute_counts(model::Potts, i::Int, j::Int)
	counts = zeros(Int, model.params.q)

	counts[model.σ[plus(i, model.params.L), j]]  += 1
	counts[model.σ[minus(i, model.params.L), j]] += 1
	counts[model.σ[i, plus(j, model.params.L)]]  += 1
	counts[model.σ[i, minus(j, model.params.L)]] += 1

	return counts
end
