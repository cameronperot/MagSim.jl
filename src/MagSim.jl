module MagSim

import Random: rand, MersenneTwister
import Statistics: mean, var

export

# AbstractModels
AbstractMagnetismModel,
Ising_2D,
Potts_2D,
XY_2D,

# analysis_methods
compute_autocorrelation_time,

# cluster_algorithms
wolff!,
swendsen_wang!,

# local_algorithms
metropolis!,
metropolis_optimized!,
heat_bath!,

# parallel_simulation_methods
simulate_βs_and_Ls,

# plotting_methods
plot_energy,
plot_heat_capacity,
plot_magnetization,
plot_susceptibility,
plot_2D_Ising_lattice

# includes
include("./AbstractModels.jl")
include("./analysis_methods.jl")
include("./cluster_algorithms.jl")
include("./local_algorithms.jl")
include("./parallel_simulation_methods.jl")
include("./plotting_methods.jl")

end # module
