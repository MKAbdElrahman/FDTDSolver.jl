@with_kw struct BlackmanHarrisPulse <: AbstractWaveform
	half_bandwidth::Float64
	fc::Float64 = 0.0
end
function (pulse::BlackmanHarrisPulse)(t)
	a0 = 0.353222222
	a1 = −0.488
	a2 = 0.145
	a3 = −0.010222222
	T = 1 / pulse.half_bandwidth
	 if 0 ≤ t ≤ T
	 return (a0 + a1 * cos((2π*t)/T) + a2 * cos((2π*2*t)/T) +  a3 * cos((2π*3*t)/T)) * cos(2π*pulse.fc*(t)) 
	 else
		 return 0
	 end
end
@with_kw struct DBlackmanHarrisPulse <: AbstractWaveform
	half_bandwidth::Float64
	fc::Float64 = 0.0
end	
function (pulse::DBlackmanHarrisPulse)(t)
	a0 = 0.353222222
	a1 = −0.488
	a2 = 0.145
	a3 = −0.010222222
	T = 1 / pulse.half_bandwidth
	 if 0 ≤ t ≤ T
	 return -( a1 * sin((2π*t)/T) + a2 * sin((2π*2*t)/T) +  a3 * sin((2π*3*t)/T)) * cos(2π*pulse.fc*(t)) 
	 else
		 return 0
	 end
end	

const BlackmanHarrisWaveForm = Union{BlackmanHarrisPulse,DBlackmanHarrisPulse}	
	
@recipe function f(s::BlackmanHarrisWaveForm;nsamples = 20,spectrum = true)
	linewidth := 2
	label := "BlackmanHarris Pulse"
	T = 1/s.half_bandwidth
	if spectrum 
		xlabel := "frequency"
		xrange := [-s.fc- 10*s.half_bandwidth,s.fc+ 10*s.half_bandwidth]
	else
		xlabel := "time"
		xrange := [0,2T]
	end
	δt = 1/(s.fc+ s.half_bandwidth)/nsamples
	t_samples = range(0,20T,step = δt)
 if spectrum return fourier(t_samples,s.(t_samples)) else return (t_samples,s.(t_samples)) end
end
