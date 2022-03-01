

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
