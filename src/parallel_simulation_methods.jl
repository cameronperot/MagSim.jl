"""
	simulate_βs_and_Ls(model::AbstractMagnetismModel, βs, Ls)

"""
function simulate_Ising_parallel(
	algorithm ::Function,
	Ls        ::Array{Int, 1},
	βs        ::Array{Float64, 1};
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
