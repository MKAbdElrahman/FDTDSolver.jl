
export optimal_cell_size

export StructuredGrid, YeeGrid

export xnodes, ynodes, znodes,nodes
export cell_centroids, cell_sizes,cell_volumes

export Primal , Dual


include("./Grid/UnifromGridGeneration.jl")
include("./Grid/StructuredGrid.jl")
include("./Grid/YeeGrid.jl")
