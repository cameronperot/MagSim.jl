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


function compute_autocorrelation_time(model::AbstractMagnetismModel, ns::Array{Int, 1})
	e       = model.observables.E ./ length(model.σ)
	m       = abs.(model.observables.M) ./ length(model.σ)
	N       = length(e)
	σ²_e    = mean(e.^2) - mean(e)^2
	σ²_m    = mean(m.^2) - mean(m)^2
	τ_int_e = []
	τ_int_m = []

	for n in ns
		N_B        = Int(floor(N / n))
		e_bin_avgs = compute_bin_averages(e, n)
		m_bin_avgs = compute_bin_averages(m, n)
		ε²_e       = sum((e_bin_avgs .- mean(e_bin_avgs)).^2) / (N_B * (N_B - 1))
		ε²_m       = sum((m_bin_avgs .- mean(m_bin_avgs)).^2) / (N_B * (N_B - 1))
		push!(τ_int_e, ε²_e * N / 2σ²_e)
		push!(τ_int_m, ε²_m * N / 2σ²_m)
	end

	return Dict("e" => τ_int_e, "m" => τ_int_m)
end
