"""
	Observables()

Type containing the observables data for a magnetism simulation.

Attributes
* `E_current`       : The current energy of the system
* `E`               : Array tracking the energy of the system over the updates
* `M_current`       : The current magnetization of the system
* `M`               : Array tracking the magnetization of the system over the updates
* `avg_cluster_size`: Average cluster size (applies only to cluster algorithms)
* `statistics`      : Dictionary containing the per spin values of average energy, average
	magnetization, heat capacity, magnetic susceptibility, Binder parameter, and the
	integrated autocorrelation time (computed using the binning method)
"""
mutable struct Observables
	E_current        ::Int
	E                ::Array{Int, 1}
	M_current        ::Int
	M                ::Array{Int, 1}
	avg_cluster_size ::Float64
	statistics       ::Dict

	function Observables()
		new(0, [], 0, [], 0, Dict())
	end
end


"""
	XYObservables()

Type containing the observables data for an XY simulation.

Attributes
* `E_current`       : The current energy of the system
* `E`               : Array tracking the energy of the system over the updates
* `Mx_current`      : The current x magnetization of the system
* `Mx`              : Array tracking the x magnetization of the system over the updates
* `My_current`      : The current y magnetization of the system
* `My`              : Array tracking the y magnetization of the system over the updates
* `avg_cluster_size`: Average cluster size (applies only to cluster algorithms)
* `statistics`      : Dictionary containing the per spin values of average energy, average
	magnetization, heat capacity, magnetic susceptibility, Binder parameter, and the
	integrated autocorrelation time (computed using the binning method)
"""
mutable struct XYObservables
	E_current        ::Float64
	E                ::Array{Float64, 1}
	Mx_current       ::Float64
	Mx               ::Array{Float64, 1}
	My_current       ::Float64
	My               ::Array{Float64, 1}
	avg_cluster_size ::Float64
	statistics       ::Dict

	function XYObservables()
		new(0, [], 0, [], 0, [], 0, Dict())
	end
end
