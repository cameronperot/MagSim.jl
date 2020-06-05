# MagSim.jl

This is a Julia package for simulating ferromagnetic systems, e.g. simulation of the Ising model using the Metropolis algorithm.
More information on the underlying theory and motivation can be found in my [blog series on simulating magnetism](https://cameronperot.com/physics/2019/12/29/simulating-magnetism-part-1.html).

This is still a work in progress, more models/algorithms and documentation will be added in the future.

## Usage

```julia
using MagSim


L = 16;
β = log(1 + √2) / 2;
n_sweeps = 10^5;

model = Ising(L, β, n_sweeps=n_sweeps);
metropolis!(model);
println("Metropolis algorithm")
println("Energy per spin: ", model.observables.statistics[:e])
println("Magnetization per spin: ", model.observables.statistics[:m])

model = Ising(L, β, n_sweeps=n_sweeps);
heat_bath!(model);
println("Heat-bath algorithm")
println("Energy per spin: ", model.observables.statistics[:e])
println("Magnetization per spin: ", model.observables.statistics[:m])

model = Ising(L, β, n_sweeps=n_sweeps);
wolff!(model, 7);
println("Wolff single cluster algorithm")
println("Energy per spin: ", model.observables.statistics[:e])
println("Magnetization per spin: ", model.observables.statistics[:m])

model = Ising(L, β, n_sweeps=n_sweeps);
swendsen_wang!(model);
println("Swendsen-Wang multiple cluster algorithm")
println("Energy per spin: ", model.observables.statistics[:e])
println("Magnetization per spin: ", model.observables.statistics[:m])
```

Output:
```
Metropolis algorithm
Energy per spin: -1.45337421875
Magnetization per spin: 0.71599599609375
Heat-bath algorithm
Energy per spin: -1.450916015625
Magnetization per spin: 0.71159150390625
Wolff single cluster algorithm
Energy per spin: -1.4525763671875
Magnetization per spin: 0.71304951171875
Swendsen-Wang multiple cluster algorithm
Energy per spin: -1.45498671875
Magnetization per spin: 0.71501328125
```
