module Signal
    include("functions.jl")
    include("transformations.jl")
    include("field.jl")
    include("lines.jl")
    include("sparks.jl")

    using LinearAlgebra

    function generate_signal(psr)
        # line of sight points at the polar cap
        los_points = []

        for line in psr.line_of_sight
            push!(los_points, [line[1][end], line[2][end], line[3][end]])
        end
        # pulsar radio signal
        signal_number = size(psr.sparks_locations)[1]
        bin_number = size(los_points)[1]
        psr.signal = zeros(signal_number, bin_number)
        sigma = psr.spark_radius / 3.72 # 2.355->FWHM, 3.03->1%, 3.72->0.1%   
        for (i, p) in enumerate(los_points)
            for (j, sparks) in enumerate(psr.sparks_locations)
                for s in sparks
                    dist = norm(p - s)
                    psr.signal[j, i] += exp(-dist^2 / (2 * sigma^2)) 
                end
            end
        end

    end
    function generate_pulses(psr)

        signal_number, bin_number = size(psr.signal)

        pulse_number = 100

        psr.pulses = zeros(pulse_number, bin_number)

        step_skip = 10 # TODO work here
        pulse_idx = 1

        for i in 1:signal_number
            if i % step_skip == 0
                psr.pulses[pulse_idx, :] = psr.signal[i, :]
                pulse_idx += 1
            end
            if pulse_idx > pulse_number
                break
            end
        end

        #println("sig. $signal_number  bin $bin_number puls. $pulse_number")


        
    end

end # module end