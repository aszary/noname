module Lines
    using LinearAlgebra

    include("field.jl")
    include("functions.jl")
    include("transformations.jl")

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

        theta_max = Functions.theta_max(psr.r_em/psr.r, psr) # TODO TODO 
        rho_approx = 3 / 2 * theta_max           # beam half-opening angle
        rho = Functions.rho(theta_max)

        println("beta $(psr.beta)")
        println("theta_max $theta_max")
        println("rho_approx $rho_approx")
        println("rho $rho")

        println()
        psr.line_of_sight = []
        phis = range(0, 2π, length=1000)

        for phi in phis
            vec = Transformations.beaming(Functions.spherical2cartesian(psr.rotation_axis), deg2rad(psr.alpha+psr.beta), phi)
            push!(psr.line_of_sight, Functions.spherical2cartesian(vec))
        end

        #rho = deg2rad(2.0)
        #println(rho)
        r = psr.r
        #psr.line_of_sight = observer_line_in_beam(omega, mu, beta, rho, r, n_points=50)
        
        #psr.line_of_sight = footpoints_in_beam(omega, mu, beta, theta, r, 500_00 ,n_points=50)
        #res = line_of_sight_coords(psr.alpha, psr.beta, rad2deg(rho))
        #println(res)
    end

function line_of_sight_coords(alpha_deg, beta_deg, rho_deg)
    # konwersja na radiany
    α   = deg2rad(alpha_deg)
    β   = deg2rad(beta_deg)
    ρ   = deg2rad(rho_deg)
    ζ   = α + β

    sinα = sin(α)
    sinζ = sin(ζ)

    denom = sinα * sinζ
    num   = cos(ρ) - cos(α)*cos(ζ)

    if abs(denom) < 1e-12
        error("Przypadek osobliwy: sin(alpha)*sin(zeta) ≈ 0 (osiowa geometria).")
    end

    arg = num / denom

    if arg < -1 || arg > 1
        return (
            visible = false,
            reason  = "No real intersection: |arg| > 1",
            arg     = arg
        )
    end

    φ0 = acos(arg)

    # wektory kierunku obserwatora dla faz ±φ0
    s_plus  = (
        sin(ζ)*cos(φ0),
        sin(ζ)*sin(φ0),
        cos(ζ)
    )
    s_minus = (
        sin(ζ)*cos(φ0),
       -sin(ζ)*sin(φ0),
        cos(ζ)
    )

    return (
        visible          = true,
        phi0_rad         = φ0,
        phi0_deg         = rad2deg(φ0),
        pulse_width_deg  = 2*rad2deg(φ0),
        s_plus           = s_plus,
        s_minus          = s_minus,
        arg              = arg
    )
end





end # module end