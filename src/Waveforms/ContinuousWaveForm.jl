@with_kw struct ContinuousWave <: AbstractWaveform
    fc::Float64
    t0::Float64 = 6 / fc
end

function (s::ContinuousWave)(t)
    g(t) = sin(π / 2 * (t / s.t0))
    0 ≤ t ≤ s.t0 ? g(t) * sin(2π * s.fc * (t - s.t0)) : sin(2π * s.fc * (t - s.t0))
end

@recipe function f(s::ContinuousWave; nsamples = 30, spectrum = true)
    linewidth := 2
    label := "Sinusoidal"
    δt = 1 / s.fc / nsamples
    t_samples = range(0, s.t0 + 10 * 1 / s.fc, step = δt)
    if spectrum
        return fourier(t_samples, s.(t_samples))
    else
        return (t_samples, s.(t_samples))
    end
end

