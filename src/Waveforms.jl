export ContinuousWave
export GaussianPulse, DGaussianPulse
export BlackmanHarrisPulse, DBlackmanHarrisPulse


function fourier(t,signal)
	δt = step(t) ; N = length(t)
	signal_fft_mag = abs.(fft(signal))* δt 
	freqs = fftfreq(N,1/δt)
	return freqs, signal_fft_mag
end

abstract type AbstractWaveform end

include("./Waveforms/ContinuousWaveForm.jl")
include("./Waveforms/GaussianWaveForm.jl")
include("./Waveforms/BlackmanHarrisWaveForm.jl")




