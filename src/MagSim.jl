module MagSim

using CSV
using DataFrames
using Random
using Statistics

import Base: show

export

    # AbstractModels
    MagnetismModels,
    Ising,
    Potts,
    XY,

    # analysis_methods
    sim_dict_to_df,

    # cluster_algorithms
    wolff!,
    swendsen_wang!,

    # local_algorithms
    metropolis!,
    heat_bath!,

    # parallel_simulation_methods
    simulate_Ising_parallel,
    simulate_Potts_parallel,

    # plotting_methods
    plot_energy,
    plot_heat_capacity,
    plot_magnetization,
    plot_susceptibility,
    plot_Ising_lattice

# includes
include("./Observables.jl")
include("./Parameters.jl")
include("./MagnetismModels.jl")
include("./analysis_methods.jl")
include("./cluster_algorithms.jl")
include("./initialization_methods.jl")
include("./local_algorithms.jl")
include("./misc_methods.jl")
include("./parallel_simulation_methods.jl")
include("./plotting_methods.jl")
include("./update_methods.jl")

end # module
