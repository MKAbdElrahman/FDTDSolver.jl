ref_grid = StructuredGrid(LinRange(-5, 5, 30), LinRange(-5, 7, 10), LinRange(-5, 1, 100))

e = VectorField{Float64}(YeeGrid{Primal}(ref_grid))
d = VectorField{Float64}(YeeGrid{Primal}(ref_grid))
ϵ = TensorField{Float64}(YeeGrid{Primal}(ref_grid))

b = VectorField{Float64}(YeeGrid{Dual}(ref_grid))
h = VectorField{Float64}(YeeGrid{Dual}(ref_grid))
μ = TensorField{Float64}(YeeGrid{Dual}(ref_grid))

j = VectorField{Float64}(YeeGrid{Primal}(ref_grid))
σₑ = TensorField{Float64}(YeeGrid{Primal}(ref_grid))

m = VectorField{Float64}(YeeGrid{Dual}(ref_grid))
σₘ = TensorField{Float64}(YeeGrid{Dual}(ref_grid))

cf = VectorField{Float64}(YeeGrid{Dual}(ref_grid))
ca = VectorField{Float64}(YeeGrid{Primal}(ref_grid))


FaradayStep!(b, e, m, cf, ref_grid = ref_grid)

_update_dual_constitutive!(h, b, μ, size(ref_grid)...)
_update_dual_constitutive!(m, h, σₘ, size(ref_grid)...)


AmpereStep!(d, h, j, ca, ref_grid = ref_grid)

_update_primal_constitutive!(e, d, ϵ, size(ref_grid)...)
_update_primal_constitutive!(j, e, σₑ, size(ref_grid)...)