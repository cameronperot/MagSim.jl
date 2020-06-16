"""
    get_neighbor_indices(model::Ising, (i, j)::NTuple{2, Int})

Returns the indices of the nearest neighbors for a given lattice site (i, j).

Arguments
* `model`  : AbstractMagnetismModel type
* `(i, j)` : Two-tuple containing the row and column indices of a lattice site
Returns
* Four-tuple of two-tuples representing indices of the neighbors of site (i, j)
"""
function get_neighbor_indices(model::AbstractMagnetismModel, (i, j)::NTuple{2,Int})
    return (
        (plus(i, model.params.L), j),
        (minus(i, model.params.L), j),
        (i, plus(j, model.params.L)),
        (i, minus(j, model.params.L)),
    )
end


"""
    sum_of_neighbors(model::Union{Ising, XY}, (i, j)::NTuple{2, Int})

Sums the values of the nearest neighbors for a given lattice site (i, j).

Arguments
* `model`  : Ising or XY type
* `(i, j)` : Two-tuple containing the row and column indices of a lattice site
Returns
* Sum of the spin values of the nearest-neighbors of site (i, j)
"""
function sum_of_neighbors(model::Union{Ising,XY}, (i, j)::NTuple{2,Int})
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

Arguments
* `a`: Two-tuple representing a 2D vector
* `b`: Two-tuple representing a 2D vector
Returns
* Dot product between `a` and `b`
"""
function dot(a::NTuple{2,Float64}, b::NTuple{2,Float64})
    return a[1] * b[1] + a[2] * b[2]
end



"""
    compute_counts(model::Potts, (i, j)::NTuple{2, Int})

Computes the counts/frequency of the different spin vaues neighboring `σ[i, j]`.

Arguments
* `model`  : Potts type
* `(i, j)` : Two-tuple containing the row and column indices of a lattice site
Returns
* `counts` : Array where `counts[i]` is the number of neighbors of (i, j) with spin value i
"""
function compute_counts(model::Potts, (i, j)::NTuple{2,Int})
    counts = zeros(Int, model.params.q)

    counts[model.σ[plus(i, model.params.L), j]] += 1
    counts[model.σ[minus(i, model.params.L), j]] += 1
    counts[model.σ[i, plus(j, model.params.L)]] += 1
    counts[model.σ[i, minus(j, model.params.L)]] += 1

    return counts
end


"""
    flip_spin(model::XY, (i, j)::NTuple{2, Int}, r::NTuple{2, Float64})

Flips an XY model spin w.r.t. a random vector `r`.

Arguments
* `model`  : XY type
* `(i, j)` : Two-tuple containing the row and column indices of a lattice site
* `r`      : The vector to flip the spin vector at (i, j) over
Returns
* The spin vector at site (i, j) flipped over the `r` vector
"""
function flip_spin(model::XY, (i, j)::NTuple{2,Int}, r::NTuple{2,Float64})
    return model.σ[i, j] .- 2 * dot(model.σ[i, j], r) .* r
end


"""
    flip_spin!(model::XY, (i, j)::NTuple{2, Int}, r::NTuple{2, Float64})

Flips an XY model spin w.r.t. a random vector `r`.

Arguments
* `model`  : XY type
* `(i, j)` : Two-tuple containing the row and column indices of a lattice site
* `r`      : The vector to flip the spin vector at (i, j) over
"""
function flip_spin!(model::XY, (i, j)::NTuple{2,Int}, r::NTuple{2,Float64})
    model.σ[i, j] = flip_spin(model, (i, j), r)
end


"""
    minus(i::Int, L::Int)

Returns the "minus" nearest neighbor index using periodic boundary conditions.

Arguments
* `i`: Linear index
* `L`: Maximum linear index
Returns
* The "minus" nearest neighbor index using periodic boundary conditions.
"""
minus(i::Int, L::Int) = i ≠ 1 ? i - 1 : L


"""
    plus(i::Int, L::Int)

Returns the "plus" nearest neighbor index using periodic boundary conditions.

Arguments
* `i`: Linear index
* `L`: Maximum linear index
Returns
* The "plus" nearest neighbor index using periodic boundary conditions.
"""
plus(i::Int, L::Int) = i ≠ L ? i + 1 : 1
