mutable struct StructuredGrid{nD} <: AbstractArray{Real,nD}
    nodes::NTuple{nD,Vector{Float64}}
end

function StructuredGrid(nodes...)
    StructuredGrid{length(nodes)}(nodes)
end

Base.size(grid::StructuredGrid) = length.(grid.nodes)

nodes(sgrid::StructuredGrid) = sgrid.nodes
xnodes(grid::StructuredGrid) = grid.nodes[1]
ynodes(grid::StructuredGrid) = grid.nodes[2]
znodes(grid::StructuredGrid) = grid.nodes[3]

xnodes(grid::StructuredGrid, i::Int) = getindex(xnodes(grid), i)
ynodes(grid::StructuredGrid, j::Int) = getindex(ynodes(grid), j)
znodes(grid::StructuredGrid, k::Int) = getindex(znodes(grid), k)

Base.getindex(grid::StructuredGrid{1}, i::Int) = (xnodes(grid, i))
Base.getindex(grid::StructuredGrid{2}, i::Int, j::Int) = (xnodes(grid, i), ynodes(grid, j))
Base.getindex(grid::StructuredGrid{3}, i::Int, j::Int, k::Int) = (xnodes(grid, i), ynodes(grid, j), znodes(grid, k))

cell_sizes(sgrid::StructuredGrid) = StructuredGrid(diff.(nodes(sgrid)))
cell_volumes(sgrid::StructuredGrid) = map(prod, cell_sizes(sgrid))
_avg(nodes::AbstractVector) = 0.5(nodes[2:end] + nodes[1:end-1])
cell_centroids(sgrid::StructuredGrid) = StructuredGrid(_avg.(nodes(sgrid)))

abstract type AbstractNodeType end
struct Primal <: AbstractNodeType end
struct Dual <: AbstractNodeType end

function primal_nodes(ax::FieldComponent, reference_grid::StructuredGrid{Dim}) where {Dim}
    StructuredGrid(ntuple(i -> i == Int(ax) ? _avg(nodes(reference_grid)[i]) : nodes(reference_grid)[i], Dim))
end
function dual_nodes(ax::FieldComponent, reference_grid::StructuredGrid{Dim}) where {Dim}
    StructuredGrid(ntuple(i -> i == Int(ax) ? nodes(reference_grid)[i] : _avg(nodes(reference_grid)[i]), Dim))
end

function nodes(::Type{Primal}, ax::FieldComponent, reference_grid::StructuredGrid)
    primal_nodes(ax, reference_grid)
end
function nodes(::Type{Dual}, ax::FieldComponent, reference_grid::StructuredGrid)
    dual_nodes(ax, reference_grid)
end