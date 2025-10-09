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
    function random_sparks!(psr; min_dist=30, trials=1000)
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

        function random_sparks_grid!(psr; min_dist=30, trials=1000)
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
                    vs[i, j] = NaN
                end
            end
        end

        # calculate electric field
        ex = Array{Float64}(undef, grid_size, grid_size)
        ey = Array{Float64}(undef, grid_size, grid_size)

        # julia using diff (bad results!)
        """
        grad_vx = diff(vs, dims=1) # 199x200
        grad_vy = diff(vs, dims=2) # 200x199

        for i in 1:grid_size-1
            for j in 1:grid_size
                ex[i, j] = -grad_vx[i, j]
            end
        end
        # dirty hacks below
        for j in 1:grid_size
            ex[grid_size, j] = -grad_vx[grid_size-1, j]
        end

        for i in 1:grid_size
            for j in 1:grid_size-1
                ey[i, j] = -grad_vy[i, j]
            end
        end
        for j in 1:grid_size
            ey[j, grid_size] = -grad_vy[j, grid_size-1]
        end
        """

        # python gradient calculation
        # TODO find julia solution
        np = pyimport("numpy")
        grad_v2 = np.gradient(vs)
        grad_vx = grad_v2[1]
        grad_vy = grad_v2[2]
        ex = - grad_vx
        ey =  - grad_vy

        # calculate drift velocity
        vdx = Array{Float64}(undef, grid_size, grid_size)
        vdy = Array{Float64}(undef, grid_size, grid_size)
        for i in 1:grid_size
            for j in 1:grid_size
                B = Field.bd(gr[1][i], gr[1][j], psr)
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
    Calculates electric potential, electric field and drift velocity for custom sparks defined by (sx, sy, sz)...
    """
    function calculate_potential_custom!(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = psr.sparks
        spark_num = size(sp)[1]

        vs = Array{Float64}(undef, grid_size, grid_size)

        v_min =  1e50
        v_max =  -1e50

        # TODO start here
        for i in 1:grid_size
            for j in 1:grid_size
                vv = 0
                for k in 1:spark_num
                    sx = sp[k][1]
                    sy = sp[k][2]
                    sz = sp[k][3]
                    dist = norm([gr[1][i], gr[2][j], gr[3][i, j]] - [sx, sy, sz])
                    #vv += v(dist) # nice looking dots (Inf) in the plot, but no
                    if dist != 0
                        vv += v(dist)
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
                    vs[i, j] = NaN
                end
            end
        end

        # calculate electric field
        ex = Array{Float64}(undef, grid_size, grid_size)
        ey = Array{Float64}(undef, grid_size, grid_size)

        # python gradient calculation
        # TODO find julia solution check calculate_potential!() for inspiration
        np = pyimport("numpy")
        grad_v2 = np.gradient(vs)
        grad_vx = grad_v2[1]
        grad_vy = grad_v2[2]
        ex = - grad_vx
        ey =  - grad_vy

        # calculate drift velocity
        vdx = Array{Float64}(undef, grid_size, grid_size)
        vdy = Array{Float64}(undef, grid_size, grid_size)
        for i in 1:grid_size
            for j in 1:grid_size
                B = Field.bd(gr[1][i], gr[1][j], psr)
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


end # module