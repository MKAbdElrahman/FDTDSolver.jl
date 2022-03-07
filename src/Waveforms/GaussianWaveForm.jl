@with_kw struct DGaussianPulse <: AbstractWaveform
    fc::Float64 = 0.0
    half_bandwidth::Float64
    half_width::Float64 = 1 / (2 * sqrt(log(2))) * 1 / (π * half_bandwidth)
    time_delay::Float64 = 6 * half_width
end
(s::DGaussianPulse)(t) = -(t - s.time_delay) / s.half_width * exp(-(t - s.time_delay)^2 / s.half_width^2)
@with_kw struct GaussianPulse <: AbstractWaveform
    fc::Float64 = 0.0
    half_bandwidth::Float64
    half_width::Float64 = 1 / (2 * sqrt(log(2))) * 1 / (π * half_bandwidth)
    time_delay::Float64 = 6 * half_width
end

(s::GaussianPulse)(t) = cos(2π * s.fc * (t - s.time_delay)) * exp(-(t - s.time_delay)^2 / s.half_width^2)

const GaussianWaveForm = Union{GaussianPulse,DGaussianPulse}

@recipe function f(s::GaussianWaveForm; nsamples = 20, spectrum = true)
    if spectrum
        xlabel := "frequency"
        xrange := [-s.fc - 10 * s.half_bandwidth, s.fc + 10 * s.half_bandwidth]
    else
        xrange := [s.time_delay - 3s.half_width, s.time_delay + 3s.half_width]
    end
    δt = 1 / (s.fc + s.half_bandwidth) / nsamples
    t_samples = range(s.time_delay - 20s.half_width, s.time_delay + 20s.half_width, step = δt)
    if spectrum
        return fourier(t_samples, s.(t_samples))
    else
        return (t_samples, s.(t_samples))
    end
end
