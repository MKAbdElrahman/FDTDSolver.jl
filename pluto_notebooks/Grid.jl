### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 456df07b-1be8-4adb-9302-9290314311e3
begin
	@enum  CartesianDirection XDirection = 1 YDirection ZDirection 
	@enum  CartesianAxis XAxis = 1 YAxis ZAxis
	@enum  FieldComponent XComponent = 1 YComponent ZComponent
end

# ╔═╡ 4f7a7193-e75e-4756-9423-7d81e2f7e4db
begin

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

abstract type AbstractNodeType end
struct 	Primal <: AbstractNodeType  end
struct 	Dual <: AbstractNodeType  end
	
function primal_nodes(ax::FieldComponent, reference_grid::StructuredGrid{Dim}) where {Dim}
    StructuredGrid(ntuple(i -> i == Int(ax) ? _avg(nodes(reference_grid)[i]) : nodes(reference_grid)[i], Dim))
end
function dual_nodes(ax::FieldComponent, reference_grid::StructuredGrid{Dim}) where {Dim}
    StructuredGrid(ntuple(i -> i == Int(ax) ? nodes(reference_grid)[i] :  _avg(nodes(reference_grid)[i])  , Dim))
end

function nodes(::Type{Primal},ax::FieldComponent, reference_grid::StructuredGrid)
primal_nodes(ax,reference_grid)
end
function nodes(::Type{Dual},ax::FieldComponent, reference_grid::StructuredGrid)
dual_nodes(ax,reference_grid)
end

end

# ╔═╡ 1bc97871-d861-47d4-bde9-e2373b1e903a
begin	
struct YeeGrid{T,nD} 
	x::StructuredGrid{nD}
	y::StructuredGrid{nD}
	z::StructuredGrid{nD}
end
	
function YeeGrid{T}(refgrid::StructuredGrid{nD}) where {nD,T}
	YeeGrid{T,nD}(
		nodes(T,XComponent,refgrid),
		nodes(T,YComponent,refgrid),
		nodes(T,ZComponent,refgrid))
end
function YeeGrid{T}(nodes...) where T
	refgrid = StructuredGrid(nodes...)
	YeeGrid{T}(refgrid)
end

	
end

# ╔═╡ 4b3b120f-52b0-47ae-bec2-207ea12ddd6a
begin
	grid_3d_primal = YeeGrid{Primal}(LinRange(-2,2,10),LinRange(-5,5,100),LinRange(-3,3,40))
	grid_3d_dual = YeeGrid{Dual}(LinRange(-2,2,10),LinRange(-5,5,100),LinRange(-3,3,40))
	# nonuniform grid
	grid_2d_primal = YeeGrid{Primal}([1,4,7,8,9,10,15,20],LinRange(-1,1,15))
	grid_2d_dual = YeeGrid{Primal}([1,4,7,8,9,10,15,20],LinRange(-1,1,15))

	grid_1d_dual = YeeGrid{Primal}([1,4,7,8,9,10,15,20])
end

# ╔═╡ Cell order:
# ╠═456df07b-1be8-4adb-9302-9290314311e3
# ╠═4f7a7193-e75e-4756-9423-7d81e2f7e4db
# ╠═1bc97871-d861-47d4-bde9-e2373b1e903a
# ╠═4b3b120f-52b0-47ae-bec2-207ea12ddd6a
