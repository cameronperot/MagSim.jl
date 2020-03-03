"""
	simulate_Ising_parallel(
		algorithm ::Function,
		Ls        ::Array{Int, 1},
		βs        ::Union{Array{Float64}, AbstractRange{Float64}};
		n_sweeps  ::Int     =10^5,
		cutoff    ::Float64 =0.2,
		start_type::Symbol  =:cold,
		seed      ::Int     =8,
		)

Runs multiple Ising simulations in parallel.
"""
function simulate_Ising_parallel(
	algorithm ::Function,
	Ls        ::Array{Int, 1},
	βs        ::Union{Array{Float64}, AbstractRange{Float64}};
	n_sweeps  ::Int     =10^5,
	cutoff    ::Float64 =0.2,
	start_type::Symbol  =:cold,
	seed      ::Int     =8,
	)

	out = Dict{Tuple, Ising}()
	Threads.@threads for β in βs
		Threads.@threads for L in Ls
			out[(L, β)] =
				algorithm(
					Ising(
						L, β,
						n_sweeps=n_sweeps, cutoff=cutoff, start_type=start_type, seed=seed
						)
					)
		end
	end

	return out
end
