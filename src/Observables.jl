"""
    Observables()

Type containing the observables data for a magnetism simulation.

Attributes
* `E`               : Array tracking the energy of the system over the updates
* `M`               : Array tracking the magnetization of the system over the updates
* `avg_cluster_size`: Average cluster size (applies only to cluster algorithms)
* `statistics`      : Dictionary containing the per spin values of average energy, average
    magnetization, heat capacity, magnetic susceptibility, Binder parameter, and the
    integrated autocorrelation time (computed using the binning method)
"""
mutable struct Observables{T<:Real}
    E::Array{T,1}
    M::Array{Int,1}
    Mx::Array{Float64,1}
    My::Array{Float64,1}
    avg_cluster_size::Float64
    statistics::Dict

    function Observables(T::Type{<:Real})
        return new{T}([], [], [], [], 0, Dict())
    end
end
