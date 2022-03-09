module FDTDSolver

import LoopVectorization: @turbo
import Parameters: @with_kw
import RecipesBase: @recipe
import FFTW: fft,fftfreq
import Optim: optimize, ParticleSwarm
import LinearAlgebra: norm
include("EnumeTypes.jl")
include("Grid.jl")
include("Waveforms.jl")
include("Fields.jl")
include("FullAnistropicEngine/FullAnistropicEngine.jl")

end
