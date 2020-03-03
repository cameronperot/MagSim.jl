"""
	compute_bin_averages(O::Array{Float64, 1}, n::Int)

Computes the bin averages, used in computing the autocorrelation time.
"""
function compute_bin_averages(O::Array{Float64, 1}, n::Int)
	N = length(O)
	N_B = Int(floor(N / n))
	bin_avgs = Array{Float64, 1}(undef, N_B)

	for k in 1:N_B
		a = (k - 1) * n + 1
		b = k * n
		bin_avgs[k] = mean(O[a:b])
	end

	return bin_avgs
end


"""
	compute_autocorrelation_time(model::Union{Ising, Potts}, ns::Array{Int, 1})

Computes the integrated autocorrelation time using the binning method. The autocorrelation
time is the value at which the output sequence (either energy `:e` or magnetization `:m` in
the return dictionary) asymptotes at.
"""
function compute_autocorrelation_time(model::Union{Ising, Potts}, ns::Array{Int, 1})
	V       = length(model.σ)
	e       = model.observables.E ./ V
	m       = abs.(model.observables.M) ./ V
	N       = length(e)
	σ²_e    = mean(e.^2) - mean(e)^2
	σ²_m    = mean(m.^2) - mean(m)^2
	τ_int_e = []
	τ_int_m = []

	for n in ns
		N_B        = floor(Int, N / n)
		e_bin_avgs = compute_bin_averages(e, n)
		m_bin_avgs = compute_bin_averages(m, n)
		ε²_e       = sum((e_bin_avgs .- mean(e_bin_avgs)).^2) / (N_B * (N_B - 1))
		ε²_m       = sum((m_bin_avgs .- mean(m_bin_avgs)).^2) / (N_B * (N_B - 1))
		push!(τ_int_e, ε²_e * N / 2σ²_e)
		push!(τ_int_m, ε²_m * N / 2σ²_m)
	end

	return Dict(:e => τ_int_e, :m => τ_int_m)
end


"""
	compute_autocorrelation_time(model::XY, ns::Array{Int, 1})

Computes the integrated autocorrelation time using the binning method. The autocorrelation
time is the value at which the output sequence (either energy `:e` or magnetization `:m` in
the return dictionary) asymptotes at.
"""
function compute_autocorrelation_time(model::XY, ns::Array{Int, 1})
	V        = length(model.σ)
	e        = model.observables.E ./ V
	mx       = abs.(model.observables.Mx) ./ V
	my       = abs.(model.observables.My) ./ V
	N        = length(e)
	σ²_e     = mean(e.^2) - mean(e)^2
	σ²_mx    = mean(mx.^2) - mean(mx)^2
	σ²_my    = mean(my.^2) - mean(my)^2
	τ_int_e  = []
	τ_int_mx = []
	τ_int_my = []

	for n in ns
		N_B         = floor(Int, N / n)
		e_bin_avgs  = compute_bin_averages(e, n)
		mx_bin_avgs = compute_bin_averages(mx, n)
		my_bin_avgs = compute_bin_averages(my, n)
		ε²_e        = sum((e_bin_avgs .- mean(e_bin_avgs)).^2) / (N_B * (N_B - 1))
		ε²_mx       = sum((mx_bin_avgs .- mean(mx_bin_avgs)).^2) / (N_B * (N_B - 1))
		ε²_my       = sum((my_bin_avgs .- mean(my_bin_avgs)).^2) / (N_B * (N_B - 1))
		push!(τ_int_e, ε²_e * N / 2σ²_e)
		push!(τ_int_mx, ε²_mx * N / 2σ²_mx)
		push!(τ_int_my, ε²_my * N / 2σ²_my)
	end

	return Dict(:e => τ_int_e, :mx => τ_int_mx, :my => τ_int_my)
end


"""
	sim_dict_to_df(sims::Dict)

Converts a dictionary where the values are models into a DataFrame object.
"""
function sim_dict_to_df(sims::Dict)
	df = DataFrame(L=Int[], β=Float64[],
		e=Float64[], c=Float64[],
		m=Float64[], χ=Float64[],
		U=Float64[],
		)

	for sim in values(sims)
		push!(df,
			[sim.params.L, sim.params.β,
			sim.observables.statistics[:e], sim.observables.statistics[:c],
			sim.observables.statistics[:m], sim.observables.statistics[:χ],
			sim.observables.statistics[:U]]
			)
	end
	sort!(df, [:L, :β])

	return df
end
