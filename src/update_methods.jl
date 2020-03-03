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
	update_observables!(model::Union{Ising, Potts})

Pushes the current energy and magnetization to their respective arrays.
"""
function update_observables!(model::Union{Ising, Potts})
	push!(model.observables.E, model.observables.E_current)
	push!(model.observables.M, model.observables.M_current)
end


"""
	update_observables!(model::XY)

Pushes the current energy and magnetization to their respective arrays.
"""
function update_observables!(model::XY)
	push!(model.observables.E, model.observables.E_current)
	push!(model.observables.Mx, model.observables.Mx_current)
	push!(model.observables.My, model.observables.My_current)
end



"""
	compute_observables_statistics!(model::Union{Ising, Potts})

Computes the average energy per spin, average magnetization per spin,
heat capacity per spin, magnetic susceptibility per spin, Binder parameter,
and the integrated autocorrelation time (computed using the binning method).
"""
function compute_observables_statistics!(model::Union{Ising, Potts})
	V = model.params.L^2
	model.observables.statistics[:e]  = mean(model.observables.E) / V
	model.observables.statistics[:e²] = mean(model.observables.E.^2) / V^2
	model.observables.statistics[:m]  = mean(abs.(model.observables.M)) / V
	model.observables.statistics[:m²] = mean(model.observables.M.^2) / V^2
	model.observables.statistics[:c]  = model.params.β^2 * var(model.observables.E) / V
	model.observables.statistics[:χ]  = model.params.β * var(abs.(model.observables.M)) / V
	model.observables.statistics[:U]  = 1 - mean(model.observables.M.^4) / 3 / mean(model.observables.M.^2)^2
	model.observables.statistics[:τ_int] = compute_autocorrelation_time(model, collect(1:400))
end


"""
	compute_observables_statistics!(model::XY)

Computes the average energy per spin, average magnetization per spin,
heat capacity per spin, magnetic susceptibility per spin,
and the integrated autocorrelation time (computed using the binning method).
"""
function compute_observables_statistics!(model::XY)
	V = model.params.L^2
	model.observables.statistics[:e]   = mean(model.observables.E) / V
	model.observables.statistics[:e²]  = mean(model.observables.E.^2) / V^2
	model.observables.statistics[:mx]  = mean(model.observables.Mx) / V
	model.observables.statistics[:mx²] = mean(model.observables.Mx.^2) / V^2
	model.observables.statistics[:my]  = mean(model.observables.My) / V
	model.observables.statistics[:my²] = mean(model.observables.My.^2) / V^2
	model.observables.statistics[:m²]  = mean(model.observables.Mx.^2 .+ model.observables.My.^2) / V^2
	model.observables.statistics[:c]   = model.params.β^2 * var(model.observables.E) / V
	model.observables.statistics[:χ]   = V * model.observables.statistics[:m²]
	# model.observables.statistics[:τ_int] = compute_autocorrelation_time(model, collect(1:400))
end
