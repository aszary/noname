 module Signal
    using LinearAlgebra

    include("functions.jl")

    function generate_signal(psr; noise_level=0.1)
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

        # adding gaussian noise
        noise = noise_level * randn(size(psr.signal))
        psr.signal .+= noise

    end

    function generate_pulses(psr; pulse_max=100)
        # use skip_steps in simulate_sparks to have single pulses

        signal_number, bin_number = size(psr.signal)
        psr.pulses = zeros(pulse_max, bin_number)

        for i in 1:signal_number
            psr.pulses[i, :] = psr.signal[i, :]
            if i == pulse_max
                break
            end
        end

    end

    """
    Generates longitudes based on pulse width  
    """
    function generate_longitudes(psr)
        # TODO equal spacing check if it is OK

        signal_number, bin_number = size(psr.signal)

        # calculate longitudes based on profile width for given pulsar geometry
        x,y,z = psr.los_lines[end][1][1], psr.los_lines[end][2][1], psr.los_lines[end][3][1] # boundry line
        sph = Functions.cartesian2spherical([x,y,z])
        rho = rad2deg(rho_from_theta(sph[2]))
        w = Signal.pulse_width_deg(psr.alpha, psr.beta, rho)
        psr.longitudes = range(-w/2, w/2, length=bin_number) # seem hacky
    end

    """
    Generates longitudes based on pulse width for each line
    """
    function generate_longitudes1(psr)

        signal_number, bin_number = size(psr.signal)
        psr.longitudes = zeros(bin_number)
        
        println(size(psr.los_lines))

        num = size(psr.los_lines)[1]
        mid = num ÷ 2  # half

        for i in 1:num
            # calculate longitudes based on profile width for given pulsar geometry
            x,y,z = psr.los_lines[i][1][1], psr.los_lines[i][2][1], psr.los_lines[i][3][1] # boundry line
            sph = Functions.cartesian2spherical([x,y,z])
            rho = rad2deg(rho_from_theta(sph[2]))
            w = Signal.pulse_width_deg(psr.alpha, psr.beta, rho)
            if i <= mid
                psr.longitudes[i] = w / 2
            else
                psr.longitudes[i] = -w / 2
            end
            println(w)
        end
    end



"""
    pulse_width(α, β, ρ)

Oblicza szerokość profilu pulsara W na podstawie geometrii wiązki.

# Argumenty
- `α`: kąt inklinacji między osią rotacji a osią magnetyczną [radiany]
- `β`: parametr uderzenia - najbliższe przejście linii widzenia od osi magnetycznej [radiany]
- `ρ`: kąt otwarcia wiązki emisyjnej [radiany]

# Zwraca
- `W`: szerokość profilu w radianach (lub `NaN` jeśli wiązka nie jest przecięta)

# Uwagi
Używa równania: sin²(W/4) = [sin²(ρ/2) - sin²(β/2)] / [sin(α)·sin(α+β)]
Warunek |β| ≤ ρ musi być spełniony, aby wiązka była widoczna.
"""
function pulse_width(α, β, ρ)
    # Sprawdź czy wiązka jest przecięta
    if abs(β) > ρ
        return NaN
    end

    # align rotator
    if abs(α) < 1e-10
        return 2π  # 360°
    end
    
    numerator = sin(ρ/2)^2 - sin(β/2)^2
    denominator = sin(α) * sin(α + β)
    
    # Sprawdź poprawność wartości
    if denominator ≤ 0 || numerator < 0
        return NaN
    end
    
    sin2_W4 = numerator / denominator
    
    # sin² nie może przekroczyć 1
    if sin2_W4 > 1
        return 2π
    end
    
    W = 4 * asin(sqrt(sin2_W4))
    return W
end

# Wersja z kątami w stopniach
function pulse_width_deg(α_deg, β_deg, ρ_deg)
    W_rad = pulse_width(deg2rad(α_deg), deg2rad(β_deg), deg2rad(ρ_deg))
    return rad2deg(W_rad)
end


"""
    rho_from_theta(θ)

Numeric solution for ρ from θ 
"""
function rho_from_theta(θ; tol=1e-12, max_iter=50)

    ρ = 1.5 * θ
    for _ in 1:max_iter
        x = 3 / (2 * tan(ρ))
        θ_calc = atan(-x + sqrt(2 + x^2))
        Δ = θ - θ_calc
        
        if abs(Δ) < tol
            return ρ
        end
        
        ρ += 0.5 * Δ * 1.5
    end
    
    return ρ

end


end # module end