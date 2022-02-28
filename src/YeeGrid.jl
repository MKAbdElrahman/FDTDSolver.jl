export StructuredGrid, YeeGrid
export xnodes, ynodes, znodes,nodes
export cell_centroids, cell_sizes,cell_volumes

mutable struct StructuredGrid{nD} <: AbstractArray{Real,nD}
	nodes::NTuple{nD,Vector{Float64}}
end

function StructuredGrid(nodes...)
        StructuredGrid{length(nodes)}(nodes)
end

Base.size(grid::StructuredGrid)	= length.(grid.nodes)

nodes(sgrid::StructuredGrid) = sgrid.nodes	
xnodes(grid::StructuredGrid)	= grid.nodes[1]
ynodes(grid::StructuredGrid)	= grid.nodes[2]
znodes(grid::StructuredGrid)	= grid.nodes[3]
	
xnodes(grid::StructuredGrid,i::Int)	= getindex(xnodes(grid),i)
ynodes(grid::StructuredGrid,j::Int)	=  getindex(ynodes(grid),j)
znodes(grid::StructuredGrid,k::Int)	= getindex(znodes(grid),k)
	
Base.getindex(grid::StructuredGrid{1},i::Int) = (xnodes(grid,i))
Base.getindex(grid::StructuredGrid{2},i::Int,j::Int) = (xnodes(grid,i),ynodes(grid,j))
Base.getindex(grid::StructuredGrid{3},i::Int,j::Int,k::Int) =  (xnodes(grid,i),ynodes(grid,j), znodes(grid,k))

cell_sizes(sgrid::StructuredGrid) = StructuredGrid(diff.(nodes(sgrid)))
cell_volumes(sgrid::StructuredGrid) = map(prod,cell_sizes(sgrid))	
_avg(nodes::AbstractVector) = 0.5(nodes[2:end] + nodes[1:end-1])	
cell_centroids(sgrid::StructuredGrid) = StructuredGrid(_avg.(nodes(sgrid)))
function primal_nodes(ax::FieldComponent, reference_grid::StructuredGrid{Dim}) where {Dim}
    StructuredGrid(ntuple(i -> i == Int(ax) ? _avg(nodes(reference_grid)[i]) : nodes(reference_grid)[i], Dim))
end
function dual_nodes(ax::FieldComponent, reference_grid::StructuredGrid{Dim}) where {Dim}
    StructuredGrid(ntuple(i -> i == Int(ax) ? nodes(reference_grid)[i] :  _avg(nodes(reference_grid)[i])  , Dim))
end
struct VectorFieldGrid{nD}
	xcomp::StructuredGrid{nD}
	ycomp::StructuredGrid{nD}
	zcomp::StructuredGrid{nD}
end

function primalYeeGrid(refgrid::StructuredGrid)
	VectorFieldGrid(primal_nodes(XComponent,refgrid),
		primal_nodes(YComponent,refgrid),
		primal_nodes(ZComponent,refgrid))
end
function dualYeeGrid(refgrid::StructuredGrid)
	VectorFieldGrid(dual_nodes(XComponent,refgrid),
		dual_nodes(YComponent,refgrid),
		dual_nodes(ZComponent,refgrid))
end		


function prrimalYeeGrid(nodes...)
	refgrid = StructuredGrid(nodes...)
	primalYeeGrid(refgrid)
end
function dualYeeGrid(nodes...)
	refgrid = StructuredGrid(nodes...)
	dualYeeGrid(refgrid)
end



struct YeeGrid{nD}
	ref_grid::StructuredGrid{nD}
	primal::VectorFieldGrid{nD}
	dual::VectorFieldGrid{nD}
end

function YeeGrid(refgrid::StructuredGrid)
	YeeGrid(refgrid,primalYeeGrid(refgrid),dualYeeGrid(refgrid))
end

function YeeGrid(nodes...)
	refgrid = StructuredGrid(nodes...)
	YeeGrid(refgrid,primalYeeGrid(refgrid),dualYeeGrid(refgrid))
end