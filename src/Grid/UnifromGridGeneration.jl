"""
function optimal_cell_size(Δ_min, Δ_max, x_c; n_particles = 100, p_norm = Inf)

Uses particle swarm optimization to find a cell size in the range
[Δ_min,Δ_max] such that the nodes `x_c` are integer multiples of this cell size.


...
# Arguments
- `Δ_min`: minimum cell size
- `Δ_max`: maximum cell size
- `x_c`: constraint node coordinates
- `n_particles`: number of particles in the swarm (default: 100)
- `p_norm`: norm used to calculate the misfit objective function (default: Inf)
...

# Examples
```julia-repl
julia> Δ_max = 0.4;

julia> Δ_min = 0.1;

julia> x_c = [-.75,-.35,.6,.1]
4-element Vector{Float64}:
 -0.75
 -0.35
  0.6
  0.1

julia> optimal_cell_size(Δ_min,Δ_max,x_c)
[ Info: finding optimal spacing...
[ Info: optimal step found 0.12272727272727274 with loss = 0.17647058823529416
0.12272727272727274
```
...
# Details
The details of the algorithm are described in:
"Structured Mesh Generation: Open-source automatic nonuniform mesh generation for FDTD simulation."
link to the paper: https://ieeexplore.ieee.org/document/7458133
...
"""
function optimal_cell_size(Δ_min, Δ_max, x_c; n_particles = 100, p_norm = Inf,range = false)
    x0 = min(x_c...)
    xL = max(x_c...)
    Lx = xL - x0
    M_min = Lx / Δ_max
    M_max = Lx / Δ_min

    M_initial = M_max
    objective(M) = norm(abs.(M[1] / Lx * x_c - round.(M[1] / Lx * x_c)), p_norm)
    @info "finding optimal spacing..."
    op = optimize(objective, [M_initial], ParticleSwarm(; lower = [M_min], upper = [M_max], n_particles = n_particles))
    @info "optimal step found $(Lx/round.(Int,op.minimizer[1])) with loss = $(op.minimum)"
    return range ? LinRange(x0,xL,round.(Int, op.minimizer[1] + 1))  :  Lx / round.(Int, op.minimizer[1])
end