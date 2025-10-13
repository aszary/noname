module Sparks
    using LinearAlgebra
    using PyCall
    include("functions.jl")
    include("field.jl")
     """

    Random sparks at the polar cap (grids will be generated later!)

    # Arguments

    - min_dist: minimum distant in meters

    """
    function random_sparks!(psr; min_dist=30, trials=10)
        sp = []

        # maximum theta
        thm = asin((psr.r_pc - min_dist) / psr.r) # plus distance from pc boundry condition
        for i in 1:trials
            phi = rand() * 2pi
            theta = rand() * thm
            car = Functions.spherical2cartesian([psr.r, theta, phi])
            md = 2 * min_dist
            for s in sp
                dist = norm([s[1], s[2], s[3]] - [car[1], car[2], car[3]])
                if dist < md
                    md = dist
                end
            end
            if md > min_dist
                push!(sp, [car[1], car[2], car[3]])
            end
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end


    """
    simple square grid, but different data size => applicable for gradient calculation

    # try this https://stackoverflow.com/questions/40338386/calculating-a-3d-gradient-with-unevenly-spaced-points next time?
    """
    function create_grid!(psr; size=50)
        r = psr.r # stellar radius in meters

        x_min, x_max = extrema(psr.pc[1])
        y_min, y_max = extrema(psr.pc[2])
        z_min = psr.pc[3][1]

        xs = range(1.1*x_min, 1.1*x_max, length=size)
        ys = range(1.1*y_min, 1.1*y_max, length=size)
        gr_x = xs
        gr_y = ys
        gr_z = zeros((size, size))
        for (i, xx) in enumerate(xs)
            for (j, yy) in enumerate(ys)
                zz = sqrt(r^2 - xx^2 - yy^2)
                if zz >= z_min
                    gr_z[i, j] = zz
                else
                    gr_z[i, j] = 0  # zero used as magic number if point is outsied of polar cap
                end
            end
        end
        psr.grid = [gr_x, gr_y, gr_z]
    end

    """

    Random sparks for new grid shape (x[size], y[size], z[size, size])
    sparks only by index [i, j] -> location is in grid
    # Arguments

    - min_dist: minimum distant in meters

    """

        function random_sparks_grid!(psr; min_dist=30, trials=10)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = []
        for i in 1:trials
            ii = rand(1:grid_size)
            jj = rand(1:grid_size)
            x = gr[1][ii]
            y = gr[2][jj]
            z = gr[3][ii, jj]
            if !([ii, jj] in sp) && (z !=0) # not in sparks # remember to skip z=0
                md = 2 * min_dist
                # check distance between sparks
                for si in sp
                    sx = gr[1][si[1]]
                    sy = gr[2][si[2]]
                    sz = gr[3][si[1], si[2]]
                    dist = norm([sx, sy, sz] - [x, y, z])
                    #println(dist)
                    if dist < md
                        md = dist
                    end
                end
                for i in 1:size(psr.pc[1])[1]
                    dist = norm([x, y, z] - [psr.pc[1][i], psr.pc[2][i], psr.pc[3][i]])
                    if dist < md
                        md = dist
                    end
                end
                if md > min_dist
                    push!(sp, [ii, jj])
                end
            end
        end
        #psr.sparks = convert(Array{Float64,1}, sp)
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end

    """
    Electric potential [Filaments]

# Arguments

-r: distance from the spark forming region
    """
    function v(r; a=1)
        return a * log(r)
    end



    """
    Calculates electric potential, electric field and drift velocity for sparks defined by grid indexes...
    """
    function calculate_potential!(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = psr.sparks
        spark_num = size(sp)[1]

        vs = Array{Float64}(undef, grid_size, grid_size)

        v_min =  1e50
        v_max =  -1e50

        for i in 1:grid_size
            for j in 1:grid_size
                vv = 0
                for k in 1:spark_num
                    (ii, jj) = sp[k]
                    if (gr[3][i, j]!=0) # (ii != i) && (jj !=j) && (gr[3][i, j]!=0) # why?
                        sx = gr[1][ii]
                        sy = gr[2][jj]
                        sz = gr[3][ii, jj]
                        dist = norm([gr[1][i], gr[2][j], gr[3][i, j]] - [sx, sy, sz])
                        #vv += v(dist) # nice looking dots (Inf) in the plot, but no
                        if dist != 0
                            vv += v(dist)
                        end
                    end
                end
                if gr[3][i, j] != 0
                    vs[i, j] = vv
                else
                    vs[i, j] = 0
                end
            end
        end

        # set 0 potential (beyound the polar cap) to NaN
        for i in 1:grid_size
            for j in 1:grid_size
                if vs[i, j] == 0
                    vs[i, j] = NaN # kind of works, set to 0 for checks
                end
            end
        end

        # calculate electric field
        ex = Array{Float64}(undef, grid_size, grid_size)
        ey = Array{Float64}(undef, grid_size, grid_size)

        # julia solution
        grad_v2 = Functions.numpy_gradient_2d(vs)
        grad_vx = grad_v2[1]
        grad_vy = grad_v2[2]
        ex = - grad_vx
        ey =  - grad_vy
        # calculate drift velocity
        vdx = Array{Float64}(undef, grid_size, grid_size)
        vdy = Array{Float64}(undef, grid_size, grid_size)
        for i in 1:grid_size
            for j in 1:grid_size
                B = Field.bd(gr[1][i], gr[1][j], psr) #* 10 # * X to have longer arrows!
                E = [ex[i, j], ey[i, j], 0]
                #println(B)
                #println(E)
                v = cross(E, B)
                vdx[i, j] = v[1]
                vdy[i, j] = v[2]
            end
        end

        #println(typeof(grad_v2))
        psr.potential = vs
        psr.electric_field = [ex, ey]
        psr.drift_velocity = [vdx, vdy]
    end


    

    """
    Electric potential [Filaments]

    # Arguments

    -r: distance from the spark forming region
    """
    function init_sparks1!(psr; rfs=[0.295, 0.5], num=3, center=true)
        sp = []

        # at the center
        if center == true
            car = Functions.spherical2cartesian([psr.r, 0, 0])
            push!(sp, [car[1], car[2], car[3]])
        end

        c1 = nothing
        for (i, rf) in enumerate(rfs)
            r = psr.r_pc * rf
            thm = asin(r / psr.r)
            if i == 1
                c1 = 2 * pi * r
                for phi in range(0, 2pi, length=num+1)[1:num] # really?
                    car = Functions.spherical2cartesian([psr.r, thm, phi])
                    push!(sp, [car[1], car[2], car[3]])
                    #println(phi)
                end
            else
                # calculate track cicumference to adjust track radius
                ci = 2 * pi * r
                num_new = convert(Int, ceil(ci / c1) * num)
                ci_new = ceil(ci / c1) * c1
                r_new = ci_new / (2 * pi)
                println("Radius of track no. $i adjusted to $(r_new/psr.r_pc)")
                thm = asin(r_new / psr.r)
                for phi in range(0, 2pi, length=num_new+1)[1:num_new] # really?
                    car = Functions.spherical2cartesian([psr.r, thm, phi])
                    push!(sp, [car[1], car[2], car[3]])
                    #println(phi)
                end
            end
            #println("$r $c1")
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end


    """

    Initiate sparks at the polar cap (grids will be generated later!)
    Same number of sparks at all tracks

    # Arguments

    - rs: list of tracks radius
    - num: number of sparks at the inner track and all the rest

    """
    function init_sparks2!(psr; rfs=[0.235, 0.71], num=3, offset=true, center=true)
        sp = []
        # at the center
        if center == true
            car = Functions.spherical2cartesian([psr.r, 0, 0])
            push!(sp, [car[1], car[2], car[3]])
        end

        dphi = 2pi / num
        for (i, rf) in enumerate(rfs)
            r = psr.r_pc * rf
            thm = asin(r / psr.r)
            for phi in range(0, 2pi, length=num+1)[1:num] # really?
                if offset == true
                    phi += (i % 2) * dphi / 2
                end
                car = Functions.spherical2cartesian([psr.r, thm, phi])
                push!(sp, [car[1], car[2], car[3]])
                #println(phi)
            end
            #println("$r $c1")
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end



    """
    Initiate sparks at the polar cap (grids will be generated later!)
    "Fibonacci" distribution

    # Arguments

    - num: number of sparks at the whole polar cap
    - rfmax: fraction of the polar cap radius covered by sparks

    """
    function init_sparks3!(psr; num=10, rfmax=0.7)
        sp = []
        # TODO start here
        r = psr.r
        theta_max = rfmax * Functions.theta_max(1, psr)
        rp = r * cos(theta_max) # calculation only at the polar cap
        phi = pi * (3. - sqrt(5.))  # golden angle in radians
        for i in 1:num
            zz = r - ((i-1) / (num-1)) * (r-rp) #  # 2 * r  # z goes from r to -r
            #println(zz)
            radius = sqrt(r^2 - zz * zz)  # radius at z
            theta = phi * i  # golden angle increment
            xx = cos(theta) * radius
            yy = sin(theta) * radius
            push!(sp, [xx, yy, zz])
        end
        psr.sparks = sp
        println("Number of sparks added: ", size(sp)[1])
    end
    function create_grids!(psr, prec=0.3, grid_size=5)

        if psr.sparks == nothing
            println("Run init_sparks! frirst..")
            return
        end

        sp = psr.sparks
        gr = []

        for s in sp
            gr_x = range(s[1] - prec*2, s[1]+prec*2, length=grid_size)
            gr_y = range(s[2] - prec*2, s[2]+prec*2, length=grid_size)
            gr_z = zeros((grid_size, grid_size))
            for (i, x) in enumerate(gr_x)
                for (j, y) in enumerate(gr_y)
                    z = sqrt(psr.r^2 - x^2 - y^2)
                    gr_z[i, j] = z
                end
            end
            push!(gr, [gr_x, gr_y, gr_z])
        end
        psr.grid = gr
    end


end # module