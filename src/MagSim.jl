module MagSim

using CSV
using DataFrames
using Random: MersenneTwister, rand
using Statistics: mean, var

import Base.show

export

# AbstractModels
MagnetismModels,
Ising_2D,
Potts_2D,
XY_2D,

# cluster_algorithms
wolff!,
swendsen_wang!,

# local_algorithms
metropolis!,
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
include("./MagnetismModels.jl")
include("./analysis_methods.jl")
include("./cluster_algorithms.jl")
include("./initialization_methods.jl")
include("./local_algorithms.jl")
include("./parallel_simulation_methods.jl")
include("./plotting_methods.jl")
include("./update_methods.jl")

end # module
