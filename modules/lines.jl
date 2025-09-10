module Lines
    using LinearAlgebra

    include("field.jl")
    include("functions.jl")

    function calculate_polarcaps!(psr; phi_num=10)
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
                #println(i, " ", j, " ", pos, " ", pos_sph)
                # TODO start here..
                push!(psr.open_lines[i], [[pos[1]], [pos[2]], [pos[3]]]) # adding initial position
                ml = psr.open_lines[i][j] # magnetic line add values to ml[1] is x ml[2] is y ml[3] is z 
                #println("\t$(pos[1]) ", size(psr.open_lines[i][j]), ml, " ", ml[2])
                for k in 1:stepsnum
                    pos_sph = Functions.cartesian2spherical(pos)
                    b_sph = Field.bvac(pos_sph, psr.r, fv.beq)
                    b = Functions.vec_spherical2cartesian(pos_sph, b_sph)
                    st = b / norm(b) * step
                    pos += st # new position for magnetic line
                    push!(ml[1], pos[1])
                    push!(ml[2], pos[2])
                    push!(ml[3], pos[3])
                    #println(k)

                end
                #println(ml)

                #return

            end
            
        end
    end


end # module end