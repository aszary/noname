module Tools


    """
    Calculates average profile
    """
    function average_profile(data; norm=true)
        pulses, bins = size(data)
        ave = zeros(bins)
        for i in 1:pulses
            for j in 1:bins
                ave[j] += data[i,j]
            end
        end
        if norm == true
            ma = maximum(ave)
            ave = ave / ma
        elseif length(norm) == 2
            me = mean(ave[norm[1]:norm[2]])
            ave .-= me
            ma = maximum(ave)
            ave = ave ./ ma
        end
        return ave
    end


    """
    Calculates intensity of single pulses
    """
    function pulses_intensity(data)
        pulse_num, bins = size(data)
        intensity = zeros(pulse_num)
        for i in 1:pulse_num
            intensity[i] = sum(data[i, :])
        end
        pulses = collect(1:pulse_num)
        mi = minimum(intensity)
        if mi < 0
            intensity .+= -mi
        end
        intensity /= maximum(intensity)
        return (intensity, pulses)
    end


end # module end