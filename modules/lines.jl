module Lines
    using LinearAlgebra

    include("field.jl")
    include("functions.jl")

    function calculate_polarcaps!(psr; phi_num=100)
        theta = Functions.theta_max(1, psr)
        phis = range(0, 2*pi, length=phi_num)
        x = Array{Float64}(undef, phi_num)
        y = Array{Float64}(undef, phi_num)
        z = Array{Float64}(undef, phi_num)
        x2 = Array{Float64}(undef, phi_num)
        y2 = Array{Float64}(undef, phi_num)
        z2 = Array{Float64}(undef, phi_num)
    
        for (i,ph) in enumerate(phis)
            ca = Functions.spherical2cartesian([psr.r, theta, ph])
            x[i] = ca[1]
            y[i] = ca[2]
            z[i] = ca[3]
            ca2 = Functions.spherical2cartesian([psr.r, pi - theta, ph]) # south pole
            x2[i] = ca2[1]
            y2[i] = ca2[2]
            z2[i] = ca2[3]
        end
        psr.polar_caps = [[x, y, z], [x2, y2, z2]]
        psr.pc = [x, y, z]
    end

    function generate_open!(psr, step=10, stepsnum=2000)
        fv = psr.fields

        # two polar caps
        for (i,pc) in enumerate(psr.polar_caps)
            # points at the polar cap
            xs = pc[1]
            ys = pc[2]
            zs = pc[3]
            # goint other direction at south pole
            if i == 2
                step = -step
            end
            # all points
            for (j,x) in enumerate(xs) 
                pos = [xs[j], ys[j], zs[j]]
                pos_sph = Functions.cartesian2spherical(pos)
                b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                push!(psr.open_lines[i], [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                ml = psr.open_lines[i][j] # magnetic line add values to ml[1] is x ml[2] is y ml[3] is z 
                for k in 1:stepsnum
                    pos_sph = Functions.cartesian2spherical(pos)
                    b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                    b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                    st = b / norm(b) * step
                    pos += st # new position for magnetic line
                    push!(ml[1], pos[1])
                    push!(ml[2], pos[2])
                    push!(ml[3], pos[3])
                end
            end
        end
    end

    function init_line_of_sight(psr)
        omega = Functions.spherical2cartesian(psr.rotation_axis)
        mu = Functions.spherical2cartesian(psr.magnetic_axis)
        beta = deg2rad(0.4)           # impact parameter
        theta = Functions.theta_max(1, psr) # TODO TODO 
        rho = 3 / 2 * theta           # beam half-opening angle
        println("beta $beta")
        println("theta $theta")
        println("rho $rho")
        #rho = deg2rad(2.0)
        #println(rho)
        r = psr.r
        #psr.line_of_sight = observer_line_in_beam(omega, mu, beta, rho, r, n_points=50)
        
        psr.line_of_sight = footpoints_in_beam(omega, mu, beta, rho, r, 500_00 ,n_points=50)

    end
    

function observer_line_in_beam(Ω::Vector{T}, μ::Vector{T}, β::Real, ρ::Real, r::Real; 
                                n_points::Int=100) where T<:Real
    # Normalize input vectors
    Ω_hat = Ω / norm(Ω)
    μ_hat = μ / norm(μ)

    # Magnetic inclination angle (angle between rotation and magnetic axes)
    cos_α = dot(Ω_hat, μ_hat)
    α = acos(clamp(cos_α, -1, 1))
    
    # Observer angle (angle between rotation axis and line of sight)
    ζ = α + β
    
    # Check if observer intersects the beam at all
    # Minimum angle between observer and magnetic axis is |beta|
    if abs(β) > ρ
        return Vector{T}[]  # observer does not see the beam
    end
    
    # Calculate phase range where observer is inside the beam
    # From: cos(rho) = cos(alpha)*cos(zeta) + sin(alpha)*sin(zeta)*cos(phi)
    cos_phi_max = (cos(ρ) - cos(α) * cos(ζ)) / (sin(α) * sin(ζ))
    cos_phi_max = clamp(cos_phi_max, -1, 1)
    phi_max = acos(cos_phi_max)
    
    # First perpendicular vector: in the (Omega, mu) plane, orthogonal to Omega
    e1 = μ_hat - cos_α * Ω_hat
    e1 = e1 / norm(e1)
    
    # Second perpendicular vector: completes the right-handed system
    e2 = cross(Ω_hat, e1)
    
    # Generate points along the arc inside the beam
    phi_range = range(-phi_max, phi_max, length=n_points)
    
    # Points at distance r from the center
    points = [r * (Ω_hat * cos(ζ) + sin(ζ) * (e1 * cos(φ) + e2 * sin(φ))) 
              for φ in phi_range]
    
    return points
end

    

function footpoints_in_beam(Ω::Vector{T}, μ::Vector{T}, β::Real, ρ::Real,
                            R_NS::Real, h::Real; n_points::Int=100) where T<:Real
    Ω_hat, μ_hat = Ω / norm(Ω), μ / norm(μ)
    
    cos_α = dot(Ω_hat, μ_hat)
    α = acos(clamp(cos_α, -1, 1))
    ζ = α + β
    
    abs(β) > ρ && return Vector{T}[]
    
    cos_ϕ_max = clamp((cos(ρ) - cos(α) * cos(ζ)) / (sin(α) * sin(ζ)), -1, 1)
    ϕ_range = range(-acos(cos_ϕ_max), acos(cos_ϕ_max), length=n_points)
    
    e1 = normalize(μ_hat - cos_α * Ω_hat)
    e2 = cross(Ω_hat, e1)
    
    r_em = R_NS + h
    
    map(ϕ_range) do ϕ
        n_hat = Ω_hat * cos(ζ) + sin(ζ) * (e1 * cos(ϕ) + e2 * sin(ϕ))
        cos_θ = clamp(dot(n_hat, μ_hat), -1, 1)
        θ_surf = asin(clamp(sqrt(R_NS * sin(acos(cos_θ))^2 / r_em), 0, 1))
        n_perp = normalize(n_hat - cos_θ * μ_hat)
        R_NS * (μ_hat * cos(θ_surf) + n_perp * sin(θ_surf))
    end
end



end # module end