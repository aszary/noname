 module Signal
    using LinearAlgebra

    function generate(psr)
        # line of sight points at the polar cap
        los_points = []

        for line in psr.los_lines
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


end # module end