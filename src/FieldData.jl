export VectorField, TensorField

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