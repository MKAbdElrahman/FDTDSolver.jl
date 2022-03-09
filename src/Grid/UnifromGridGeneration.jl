"""
    optimal_cell_size(Δ_min, Δ_max, x_c; n_particles = 100, p_norm = Inf)

Uses particle swarm optimization to return a cell size in the range
[Δ_min,Δ_max] such that the nodes `x_c` are nearly integer multiples of the optimized cell size.

# Arguments
- `Δ_min::Float64`: minimum cell size
- `Δ_max::Float64`: maximum cell size
- `x_c::Vector{Float64}`: constraint node coordinates
- `n_particles::Int`: number of particles in the swarm (default: 100)
- `p_norm::Int`: norm used to calculate the misfit objective function (default: `Inf`)
- `range::Bool`: set `true` to return a `LinRange` instead  of cell spacing.

# Examples
```julia-repl
julia> optimal_cell_size(0.1,.3,[-1,-.5,1,1.25,2])
[ Info: finding optimal cell size...
[ Info: optimal cell size found is 0.125 with loss = 0.0
0.125
```

# References
The details of the algorithm are described in: 
[Structured Mesh Generation: Open-source automatic nonuniform mesh generation for FDTD simulation](https://ieeexplore.ieee.org/document/7458133)
"""
function optimal_cell_size(Δ_min, Δ_max, x_c; n_particles = 100, p_norm = Inf,range = false)
    x0 = min(x_c...)
    xL = max(x_c...)
    Lx = xL - x0
    M_min = Lx / Δ_max
    M_max = Lx / Δ_min

    M_initial = 0.5*(M_max + M_min)
    objective(M) = norm(M[1] / Lx * x_c - round.(M[1] / Lx * x_c), p_norm)
    @info "finding optimal cell size..."
    op = optimize(objective, [M_initial], ParticleSwarm(; lower = [M_min], upper = [M_max], n_particles = n_particles))
    @info "optimal cell size found is $(Lx/round.(Int,op.minimizer[1])) with loss = $(op.minimum)"
    return range ? LinRange(x0,xL,round.(Int, op.minimizer[1] + 1))  :  Lx / round.(Int, op.minimizer[1])
end