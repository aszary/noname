 module Signal
    using LinearAlgebra

    function generate(psr)
        # line of sight points at the polar cap
        los_points = []

        for line in psr.los_lines
            push!(los_points, [line[1][end], line[2][end], line[3][end]])
        end
        # pulsar radio signal
        pulse_number = size(psr.sparks_locations)[1]
        bin_number = size(los_points)[1]
        psr.signal = zeros(pulse_number, bin_number)
        sigma = psr.spark_radius / 3 # TODO play with this
        for (i, p) in enumerate(los_points)
            for (j, sparks) in enumerate(psr.sparks_locations)
                for s in sparks
                    dist = norm(p - s)
                    # dist = norm(p[1:2] - s[1:2]) # 2D does not help
                    # TODO somthing wrong in animation, check this
                    psr.signal[j, i] += exp(-dist^2 / (2 * sigma^2)) 
                end
            end
        end

        println(size(los_points))

    end


end # module end