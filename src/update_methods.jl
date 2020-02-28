"""
	minus(i::Int, L::Int)

Helper function for computing nearest neighbor indices with periodic boundary conditions.
"""
minus(i::Int, L::Int) = i ≠ 1 ? i-1 : L


"""
	plus(i::Int, L::Int)

Helper function for computing nearest neighbor indices with periodic boundary conditions.
"""
plus(i::Int, L::Int) = i ≠ L ? i+1 : 1


"""
	update_observables!(model::AbstractMagnetismModel, t::Int)

Pushes the current energy and magnetization to their respective arrays.
"""
function update_observables!(model::AbstractMagnetismModel, t::Int)
	push!(model.observables.E, model.observables.E_current)
	push!(model.observables.M, model.observables.M_current)
end


"""
	compute_observables_statistics!(model::AbstractMagnetismModel)

Computes the average energy per spin, average magnetization per spin,
heat capacity per spin, magnetic susceptibility per spin, Binder parameter,
and the integrated autocorrelation time (computed using the binning method).
"""
function compute_observables_statistics!(model::AbstractMagnetismModel)
	model.observables.statistics[:e] = mean(model.observables.E ./ length(model.σ))
	model.observables.statistics[:m] = mean(abs.(model.observables.M) ./ length(model.σ))
	model.observables.statistics[:c] = model.params.β^2 * var(model.observables.E) / length(model.σ)
	model.observables.statistics[:χ] = model.params.β * var(abs.(model.observables.M)) / length(model.σ)
	model.observables.statistics[:U] = 1 - mean(model.observables.M.^4) / 3 / mean(model.observables.M.^2)^2
	model.observables.statistics[:τ_int] = compute_autocorrelation_time(model, collect(1:400))
end
