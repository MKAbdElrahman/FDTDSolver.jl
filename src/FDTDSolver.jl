module FDTDSolver

using LoopVectorization
using Parameters
using RecipesBase
using FFTW
include("EnumeTypes.jl")
include("Grid.jl")
include("Waveforms.jl")
include("Fields.jl")
include("FullAnistropicEngine/FullAnistropicEngine.jl")

end
