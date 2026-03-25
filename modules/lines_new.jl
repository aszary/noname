module LinesNew


    function init_line_of_sight(psr; num=100)

         # calculates line of sight
        psr.line_of_sight = []
        psr.longitudes = zeros(num)

        α = deg2rad(psr.alpha)
        β = deg2rad(psr.beta)

        φ_s = Geometry.generate_uniform_phase_array(num, α ,β , psr.r_em, psr.p)
        #θ_array, ψ_array = Geometry.emission_points_from_phase(φ_s, α, β, psr.r_em, psr.p)
        θ_array, ψ_array = Geometry.emission_points_with_ar(φ_s, α, β, psr.r_em, psr.p)

        for i in eachindex(θ_array)
            ψ = ψ_array[i]
            θ = θ_array[i]
            push!(psr.line_of_sight, [psr.r_em, θ, ψ]) # spherical coordinates
            psr.longitudes[i] = rad2deg(φ_s[i])
        end

    end





end # module end