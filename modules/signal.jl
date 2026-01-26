module Signal
using LinearAlgebra
include("functions.jl")
include("field.jl")
include("sparks.jl")
include("transformations.jl")
include("lines.jl")


function signal_from_sparks(psr; sigma = psr.spark_radius / 3.72)

    # --- sanity checks ---
    if isempty(psr.sparks_locations)
        error("psr.sparks_locations is empty – run spark simulation first")
    end

    los_lines = psr.fields.magnetic_lines_from_los
    if isempty(los_lines)
        error("LOS magnetic lines not initialized")
    end

    # --- LOS points on polar cap ---
    los_points = Vector{Vector{Float64}}()
    for line in los_lines
        lx, ly, lz = line
        push!(los_points, [lx[end], ly[end], lz[end]])
    end

    signal_number = length(psr.sparks_locations)  # time steps
    bin_number    = length(los_points)            # LOS bins

    psr.signal = zeros(signal_number, bin_number)

    # --- main loop ---
    for (i, p) in enumerate(los_points)              # LOS bin
        for (j, sparks) in enumerate(psr.sparks_locations)  # time step
            for s in sparks                          # each spark
                dist = norm(p .- s)
                psr.signal[j, i] += exp(-dist^2 / (2 * sigma^2))
            end
        end
    end

    return psr.signal
end

function generate_pulses(psr; pulse_number = 100, step_skip = 10)

    if isempty(psr.signal)
        error("Signal not generated yet")
    end

    signal_number, bin_number = size(psr.signal)
    psr.pulses = zeros(pulse_number, bin_number)

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
end





end # module
