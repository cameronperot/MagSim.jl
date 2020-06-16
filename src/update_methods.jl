"""
    compute_E_M_Ising(σ::Array{Int8, 2})

Computes the Ising model internal energy and magnetization for the input spin configuration `σ`.

Arguments
* `σ`     : A 2D array representing a square lattice in the Ising model
Returns
* `(E, M)`: The system's energy and magnetization
"""
function compute_E_M_Ising(σ::Array{Int8,2})
    E::Int = 0
    M::Int = 0
    L = size(σ, 1)
    for j = 1:L, i = 1:L
        E -= σ[i, j] * (σ[plus(i, L), j] + σ[i, plus(j, L)])
        M += σ[i, j]
    end

    return (E, M)
end


"""
    compute_E_Potts(σ::Array{Int8, 2})

Computes the Potts model internal energy for the input spin configuration `σ`.

Arguments
* `σ` : A 2D array representing a square lattice in the Potts model
Returns
* `E` : The system's energy
"""
function compute_E_Potts(σ::Array{Int8,2})
    E::Int = 0
    L = size(σ, 1)
    for j = 1:L, i = 1:L
        E -= Int(σ[i, j] == σ[plus(i, L), j]) + Int(σ[i, j] == σ[i, plus(j, L)])
    end

    return E
end


"""
    compute_E_M_XY(σ::Array{NTuple{2, Float64}, 2})

Computes the XY model internal energy and magnetization for the input spin configuration `σ`.

Arguments
* `σ`          : A 2D array representing a square lattice in the Ising model
Returns
* `(E, Mx, My)`: The system's energy, x-magnetization, and y-magnetization
"""
function compute_E_M_XY(σ::Array{NTuple{2,Float64},2})
    E = 0.0
    Mx = 0.0
    My = 0.0
    L = size(σ, 1)
    for j = 1:L, i = 1:L
        E -= dot(σ[i, j], σ[plus(i, L), j]) + dot(σ[i, j], σ[i, plus(j, L)])
        Mx += σ[i, j][1]
        My += σ[i, j][2]
    end

    return (E, Mx, My)
end


"""
    update_observables!(model::Ising)

Pushes the current energy and magnetization to their respective arrays.

Arguments
* `model`: Ising type
"""
function update_observables!(model::Ising)
    E, M = compute_E_M_Ising(model.σ)
    push!(model.observables.E, E)
    push!(model.observables.M, M)
end


"""
    update_observables!(model::Potts)

Pushes the current energy and magnetization to their respective arrays.

Arguments
* `model`: Potts type
"""
function update_observables!(model::Potts)
    E = compute_E_Potts(model.σ)
    push!(model.observables.E, E)
end


"""
    update_observables!(model::XY)

Pushes the current energy and magnetization to their respective arrays.

Arguments
* `model`: XY type
"""
function update_observables!(model::XY)
    E, Mx, My = compute_E_M_XY(model.σ)
    push!(model.observables.E, E)
    push!(model.observables.Mx, Mx)
    push!(model.observables.My, My)
end



"""
    compute_observables_statistics!(model::AbstractMagnetismModel)

Computes the average energy per spin, average magnetization per spin,
heat capacity per spin, magnetic susceptibility per spin, Binder parameter,
and the integrated autocorrelation times (computed using the binning method).
"""
function compute_observables_statistics!(model::AbstractMagnetismModel)
    V = model.params.L^2
    model.observables.statistics[:e] = mean(model.observables.E) / V
    model.observables.statistics[:e²] = mean(model.observables.E .^ 2) / V^2
    model.observables.statistics[:c] = model.params.β^2 * var(model.observables.E) / V

    if isa(model, Ising)
        model.observables.statistics[:m] = mean(abs.(model.observables.M)) / V
        model.observables.statistics[:m²] = mean(model.observables.M .^ 2) / V^2
        model.observables.statistics[:χ] =
            model.params.β * var(abs.(model.observables.M)) / V
        model.observables.statistics[:U] =
            1 - mean(model.observables.M .^ 4) / 3 / mean(model.observables.M .^ 2)^2
    elseif isa(model, XY)
        model.observables.statistics[:mx] = mean(model.observables.Mx) / V
        model.observables.statistics[:mx²] = mean(model.observables.Mx .^ 2) / V^2
        model.observables.statistics[:my] = mean(model.observables.My) / V
        model.observables.statistics[:my²] = mean(model.observables.My .^ 2) / V^2
        model.observables.statistics[:m²] =
            mean(model.observables.Mx .^ 2 .+ model.observables.My .^ 2) / V^2
        model.observables.statistics[:χ] = model.observables.statistics[:m²] * V
    end

    model.observables.statistics[:τ_int] =
        compute_autocorrelation_time(model, collect(1:400))
end
