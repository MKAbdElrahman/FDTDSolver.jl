### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ f32b184a-ba70-4e80-8f8d-e12def369128
using Parameters

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
Base.size(yeegrid::YeeGrid) = (size(yeegrid.x),size(yeegrid.y),size(yeegrid.z))
	
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

# ╔═╡ 6ad42c82-25c5-41c3-83f1-0ff002edbafe
begin
struct VectorField{PD,T,N} 
	x::Array{T,N}
	y::Array{T,N}
	z::Array{T,N}
end
function VectorField{T}(yeegrid::YeeGrid{PD,N}) where {T,PD,N}
	VectorField{PD,T,N}(zeros.(T,size(yeegrid))...)
end
struct TensorField{PD,T,N} 
	xx::Array{T,N}
	xy::Array{T,N}
	xz::Array{T,N}
	yx::Array{T,N}
	yy::Array{T,N}
	yz::Array{T,N}
	zx::Array{T,N}
	zy::Array{T,N}
	zz::Array{T,N}
end
function TensorField{T}(yeegrid::YeeGrid{PD,N}) where {T,PD,N}
	TensorField{PD,T,N}(zeros.(T,size(yeegrid))...,zeros.(T,size(yeegrid))...,zeros.(T,size(yeegrid))...)
end	
end

# ╔═╡ 4fd9f13e-8279-41e2-b7f5-e9865ed869be
function FaradayStep!(b,e,m,cb;ref_grid)

(nx,ny,nz) = size(ref_grid)	
	
for k = 1:nz-1
	for j = 1:ny-1
		for i = 1:nx-1
b.x[i,j,k] = 
	 b.x[i,j,k] + cb.x[i,j,k] * (e.z[i,j,k]-e.z[i,j+1,k]+e.y[i,j,k+1]-e.y[i,j,k]) -m.x[i,j,k]
b.y[i,j,k] = 
	 b.y[i,j,k] + cb.y[i,j,k] * (e.x[i,j,k]-e.x[i,j,k+1]+e.z[i+1,j,k]-e.z[i,j,k]) - m.y[i,j,k]
b.z[i,j,k] = 
	 b.z[i,j,k] + cb.z[i,j,k] * (e.y[i,j,k]-e.y[i+1,j,k]+e.x[i,j+1,k]-e.x[i,j,k]) - 
	 m.z[i,j,k]
		end
	end
end
	
for k = 1:nz-1
	for j = 1:ny-1
		for i = nx:nx
b.x[i,j,k] = 
	 b.x[i,j,k] + cb.x[i,j,k] * (e.z[i,j,k]-e.z[i,j+1,k]+e.y[i,j,k+1]-e.y[i,j,k]) -m.x[i,j,k]
		end
	end
end
for k = 1:nz-1
	for j = ny:ny
		for i = 1:nx-1
b.y[i,j,k] = 
	 b.y[i,j,k] + cb.y[i,j,k] * (e.x[i,j,k]-e.x[i,j,k+1]+e.z[i+1,j,k]-e.z[i,j,k]) - m.y[i,j,k]
		end
	end
end	
for k = nz:nz
	for j = 1:ny-1
		for i = 1:nx-1
b.z[i,j,k] = 
	 b.z[i,j,k] + cb.z[i,j,k] * (e.y[i,j,k]-e.y[i+1,j,k]+e.x[i,j+1,k]-e.x[i,j,k]) - 
	 m.z[i,j,k]
		end
	end
end

end

# ╔═╡ fcb5f651-3e20-4181-9e72-9310f823adf0
function AmpereStep!(d,h,J,cb;ref_grid)

(nx,ny,nz) = size(ref_grid)	
	
for k = 2:nz-1
	for j = 2:ny-1
		for i = 2:nx-1
d.x[i,j,k] = d.x[i,j,k] +
cb.x[i,j,k] * (h.z[i,j,k]-h.z[i,j-1,k]-h.y[i,j,k]+h.y[i,j,k-1]) -J.x[i,j,k]
d.y[i,j,k] = d.y[i,j,k] +
cb.y[i,j,k] * (h.x[i,j,k]-h.x[i,j,k-1]-h.z[i,j,k]+h.z[i-1,j,k]) -J.y[i,j,k]
d.z[i,j,k] = d.z[i,j,k] +
cb.z[i,j,k] * (h.y[i,j,k]-h.y[i-1,j,k]-h.x[i,j,k]+h.z[i,j-1,k]) - J.z[i,j,k]
		end
	end
end
	
for k = 2:nz-1
	for j = 2:ny-1
		for i = 1:1
d.x[i,j,k] = d.x[i,j,k] + 
cb.x[i,j,k] * (h.z[i,j,k]-h.z[i,j-1,k]-h.y[i,j,k]+h.y[i,j,k-1]) -J.x[i,j,k]		end
	end
end	
for k = 2:nz-1
	for j = 1:1
		for i = 2:nx-1
d.y[i,j,k] = d.y[i,j,k] +
cb.x[i,j,k] * (h.x[i,j,k]-h.x[i,j,k-1]-h.z[i,j,k]+h.z[i-1,j,k]) -J.y[i,j,k]		end
	end
end
for k = 1:1
	for j = 2:ny-1
		for i = 2:nx-1
d.z[i,j,k] = d.z[i,j,k] +
cb.x[i,j,k] * (h.y[i,j,k]-h.y[i-1,j,k]-h.x[i,j,k]+h.z[i,j-1,k]) - J.z[i,j,k]
		end
	end
end	
end

# ╔═╡ 7f1e505f-7223-4ae8-a6d4-9ff455553256
function _update_primal_constitutive!(e,d,C,nx,ny,nz)
for k = 2:nz-1
	for j = 2:ny-1
		for i = 2:nx-1
dy =  0.25 * (d.y[i,j,k] + d.y[i+1,j,k] + d.y[i+1,j-1,k] + d.y[i,j-1,k] )
dz =  0.25 * (d.z[i,j,k] + d.z[i+1,j,k] + d.z[i+1,j,k-1] + d.z[i,j,k-1] ) 		
e.x[i,j,k]  = C.xx[i,j,k] * d.x[i,j,k]  + C.xy[i,j,k] * dy  + C.xz[i,j,k] * dz  


dx = 0.25*(d.x[i,j,k]+d.x[i,j+1,k]+d.x[i-1,j,k]+d.x[i-1,j+1,k])
dz = 0.25*(d.z[i,j,k] + d.z[i+1,j,k]+ d.z[i+1,j,k-1] + d.z[i,j,k-1])			
			
e.y[i,j,k]  = C.yx[i,j,k] * dx  + C.yy[i,j,k] * d.y[i,j,k]  + C.yz[i,j,k] * dz  

dx = 0.25 * (d.x[i,j,k] + d.x[i,j,k+1] + d.x[i-1,j,k+1]+d.x[i-1,j,k])
dy = 0.25 * (d.y[i,j,k] + d.y[i,j,k+1] + d.y[i,j-1,k+1]+d.y[i,j-1,k])

e.z[i,j,k]  = C.zx[i,j,k] * dx  + C.zy[i,j,k] * dy  + C.yz[i,j,k] *   d.z[i,j,k]  
			
		end
	end
end

for k = 2:nz-1
	for j = 2:ny-1
		for i = 1:1

dy =  0.25 * (d.y[i,j,k] + d.y[i+1,j,k] + d.y[i+1,j-1,k] + d.y[i,j-1,k] )
dz =  0.25 * (d.z[i,j,k] + d.z[i+1,j,k] + d.z[i+1,j,k-1] + d.z[i,j,k-1] ) 		
e.x[i,j,k]  = C.xx[i,j,k] * d.x[i,j,k]  + C.xy[i,j,k] * dy  + C.xz[i,j,k] * dz  

	end
	end
end	
for k = 2:nz-1
	for j = 1:1
		for i = 2:nx-1

dx = 0.25*(d.x[i,j,k]+d.x[i,j+1,k]+d.x[i-1,j,k]+d.x[i-1,j+1,k])
dz = 0.25*(d.z[i,j,k] + d.z[i,j+1,k]+ d.z[i,j+1,k-1] + d.z[i,j,k-1])			
			
e.y[i,j,k]  = C.yx[i,j,k] * dx  + C.yy[i,j,k] * d.y[i,j,k]  + C.yz[i,j,k] * dz  
		end
	end
end
for k = 1:1
	for j = 2:ny-1
		for i = 2:nx-1
dx = 0.25 * (d.x[i,j,k] + d.x[i,j,k+1] + d.x[i-1,j,k+1]+d.x[i-1,j,k])
dy = 0.25 * (d.y[i,j,k] + d.y[i,j,k+1] + d.y[i,j-1,k+1]+d.y[i,j-1,k])

e.z[i,j,k]  = C.zx[i,j,k] * dx  + C.zy[i,j,k] * dy  + C.yz[i,j,k] *   d.z[i,j,k]  	
		end
	end
end	
	
	
end
	

# ╔═╡ 8176e6dd-0768-4d9b-a2b7-957236e4c55b
function _update_dual_constitutive!(h,b,C,nx,ny,nz)
	for k = 2:nz-1
	for j = 2:ny-1
		for i = 2:nx-1
by = 0.25 * (b.y[i,j,k] + b.y[i-1,j,k] + b.y[i-1,j+1,k] + b.y[i,j+1,k])
bz = 0.25 * (b.z[i,j,k] + b.z[i,j,k+1] + b.z[i-1,j,k+1] + b.z[i-1,j,k])
h.x[i,j,k]  = C.xx[i,j,k] * b.x[i,j,k]  + C.xy[i,j,k] * by  + C.xz[i,j,k] * bz  


bx = 0.25*(b.x[i,j,k]+b.x[i,j-1,k]+b.x[i+1,j-1,k]+b.x[i+1,j,k])
bz = 0.25*(b.z[i,j,k] + b.z[i,j-1,k]+ b.z[i,j-1,k+1] + b.z[i,j,k+1])			
			
h.y[i,j,k]  = C.yx[i,j,k] * bx  + C.yy[i,j,k] * b.y[i,j,k]  + C.yz[i,j,k] * bz  

bx = 0.25*(b.x[i,j,k]+b.x[i+1,j,k]+b.x[i+1,j,k-1]+b.x[i,j,k-1])
by = 0.25*(b.y[i,j,k]+b.y[i,j+1,k]+b.y[i,j+1,k-1]+b.y[i,j,k-1])

h.z[i,j,k]  = C.zx[i,j,k] * bx  + C.zy[i,j,k] * by  + C.zz[i,j,k] * b.z[i,j,k]  
			 			
		end
	end
end
	
for k = 1:1
	for j = 1:1
		for i = 2:nx-1
by = 0.25 * (b.y[i,j,k] + b.y[i-1,j,k] + b.y[i-1,j+1,k] + b.y[i,j+1,k])
bz = 0.25 * (b.z[i,j,k] + b.z[i,j,k+1] + b.z[i-1,j,k+1] + b.z[i-1,j,k])
h.x[i,j,k]  = C.xx[i,j,k] * b.x[i,j,k]  + C.xy[i,j,k] * by  + C.xz[i,j,k] * bz  
		end
	end
end
for k = 1:1
	for j = 2:ny-1
		for i = 1:1
bx = 0.25*(b.x[i,j,k]+b.x[i,j-1,k]+b.x[i+1,j-1,k]+b.x[i+1,j,k])
bz = 0.25*(b.z[i,j,k] + b.z[i,j-1,k]+ b.z[i,j-1,k+1] + b.z[i,j,k+1])				
h.y[i,j,k]  = C.yx[i,j,k] * bx  + C.yy[i,j,k] * b.y[i,j,k]  + C.yz[i,j,k] * bz 
		end
	end
end	
for k = 2:ny-1
	for j = 1:1
		for i = 1:1
bx = 0.25*(b.x[i,j,k]+b.x[i+1,j,k]+b.x[i+1,j,k-1]+b.x[i,j,k-1])
by = 0.25*(b.y[i,j,k]+b.y[i,j+1,k]+b.y[i,j+1,k-1]+b.y[i,j,k-1])

h.z[i,j,k]  = C.zx[i,j,k] * bx  + C.zy[i,j,k] * by  + C.zz[i,j,k] * b.z[i,j,k] 
		end
	end
end

end

# ╔═╡ d87ef4f9-e4e9-4b1d-907d-054505d3ff54
begin
	ref_grid = StructuredGrid(LinRange(-5,15,30),LinRange(-5,17,10),LinRange(-5,1,101))

	e = VectorField{Float64}(YeeGrid{Primal}(ref_grid))
	d = VectorField{Float64}(YeeGrid{Primal}(ref_grid))
	ϵ = TensorField{Float64}(YeeGrid{Primal}(ref_grid))

	b = VectorField{Float64}(YeeGrid{Dual}(ref_grid))
	h = VectorField{Float64}(YeeGrid{Dual}(ref_grid))
	μ = TensorField{Float64}(YeeGrid{Dual}(ref_grid))

	j =  VectorField{Float64}(YeeGrid{Primal}(ref_grid)) 
	σₑ = TensorField{Float64}(YeeGrid{Primal}(ref_grid))	

	m =  VectorField{Float64}(YeeGrid{Dual}(ref_grid))
	σₘ = TensorField{Float64}(YeeGrid{Dual}(ref_grid))	

	cf = VectorField{Float64}(YeeGrid{Dual}(ref_grid))
	ca = VectorField{Float64}(YeeGrid{Primal}(ref_grid))

end

# ╔═╡ 03ac2bca-403e-478e-b947-208427fdf579
begin
	FaradayStep!(b,e,m,cf, ref_grid = ref_grid)
	_update_dual_constitutive!(h,b,μ,size(ref_grid)...)
	_update_dual_constitutive!(m,h,σₘ,size(ref_grid)...)
	
	
	AmpereStep!(d,h,j,ca, ref_grid = ref_grid)
	_update_primal_constitutive!(e,d,ϵ,size(ref_grid)...)
	_update_primal_constitutive!(j,e,σₑ,size(ref_grid)...)

end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"

[compat]
Parameters = "~0.12.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"
"""

# ╔═╡ Cell order:
# ╠═456df07b-1be8-4adb-9302-9290314311e3
# ╠═4f7a7193-e75e-4756-9423-7d81e2f7e4db
# ╠═1bc97871-d861-47d4-bde9-e2373b1e903a
# ╠═4b3b120f-52b0-47ae-bec2-207ea12ddd6a
# ╠═6ad42c82-25c5-41c3-83f1-0ff002edbafe
# ╠═f32b184a-ba70-4e80-8f8d-e12def369128
# ╠═4fd9f13e-8279-41e2-b7f5-e9865ed869be
# ╠═fcb5f651-3e20-4181-9e72-9310f823adf0
# ╠═7f1e505f-7223-4ae8-a6d4-9ff455553256
# ╠═8176e6dd-0768-4d9b-a2b7-957236e4c55b
# ╠═d87ef4f9-e4e9-4b1d-907d-054505d3ff54
# ╠═03ac2bca-403e-478e-b947-208427fdf579
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
