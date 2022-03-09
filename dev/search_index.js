var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FDTDSolver","category":"page"},{"location":"#FDTDSolver","page":"Home","title":"FDTDSolver","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FDTDSolver.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FDTDSolver]","category":"page"},{"location":"#FDTDSolver.optimal_cell_size-Tuple{Any, Any, Any}","page":"Home","title":"FDTDSolver.optimal_cell_size","text":"optimal_cell_size(Δ_min, Δ_max, x_c; n_particles = 100, p_norm = Inf)\n\nUses particle swarm optimization to return a cell size in the range [Δmin,Δmax] such that the nodes x_c are nearly integer multiples of the optimized cell size.\n\nArguments\n\nΔ_min::Float64: minimum cell size\nΔ_max::Float64: maximum cell size\nx_c::Vector{Float64}: constraint node coordinates\nn_particles::Int: number of particles in the swarm (default: 100)\np_norm::Int: norm used to calculate the misfit objective function (default: Inf)\nrange::Bool: set true to return a LinRange instead  of cell spacing.\n\nExamples\n\njulia> optimal_cell_size(0.1,.3,[-1,-.5,1,1.25,2])\n[ Info: finding optimal cell size...\n[ Info: optimal cell size found is 0.125 with loss = 0.0\n0.125\n\nReferences\n\nThe details of the algorithm are described in:  Structured Mesh Generation: Open-source automatic nonuniform mesh generation for FDTD simulation\n\n\n\n\n\n","category":"method"}]
}
