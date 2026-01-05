module Signal


    function generate_signal(psr)
        los_points = [] 
            for line in psr.los_lines
            push!(los_points, [line[1][end], line[2][end], line[3][end]])
        end
        psr.signal = zeros(size(los_points[1]))
        println(size(los_points))

    end


end # module end