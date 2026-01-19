module Sparks
    using LinearAlgebra
    using PyCall
    using JLD2
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
            
            # Extract coordinates immediately
            x = gr[1][ii]
            y = gr[2][jj]
            z = gr[3][ii, jj]
            
            # Skip points outside polar cap (z=0)
            if z != 0
                md = 2 * min_dist
                
                # Check distance against existing sparks
                for si in sp
                    dist = norm(si - [x, y, z])
                    if dist < md
                        md = dist
                    end
                end

                # Check distance against Polar Cap boundary
                if isdefined(psr, :pc) && psr.pc !== nothing
                    for k in 1:length(psr.pc[1])
                        pc_p = [psr.pc[1][k], psr.pc[2][k], psr.pc[3][k]]
                        dist = norm([x, y, z] - pc_p)
                        if dist < md; md = dist; end
                    end
                end
                
                if md > min_dist
                    push!(sp, [x, y, z]) # Storing COORDINATES
                end
            end
        end
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

        for i in 1:grid_size
            for j in 1:grid_size
                vv = 0
                for k in 1:spark_num
                    # FIX: Read coordinates directly (no (ii, jj) unpacking)
                    sx = sp[k][1]
                    sy = sp[k][2]
                    sz = sp[k][3]
                    
                    dist = norm([gr[1][i], gr[2][j], gr[3][i, j]] - [sx, sy, sz])
                    
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

        # Handle potentials beyond polar cap
        for i in 1:grid_size
            for j in 1:grid_size
                if vs[i, j] == 0
                    vs[i, j] = NaN 
                end
            end
        end

        # Calculate E-field and Drift Velocity
        grad_v2 = Functions.numpy_gradient_2d(vs)
        ex = -grad_v2[1]
        ey = -grad_v2[2]

        vdx = Array{Float64}(undef, grid_size, grid_size)
        vdy = Array{Float64}(undef, grid_size, grid_size)
        
        for i in 1:grid_size
            for j in 1:grid_size
                B = Field.bd(gr[1][i], gr[1][j], psr)
                E = [ex[i, j], ey[i, j], 0]
                v = cross(E, B)
                vdx[i, j] = v[1]
                vdy[i, j] = v[2]
            end
        end

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
"""
    Calculates electric potential, electric field and drift velocity for grids around sparks
    """
    function calculate_potentials!(psr; save=false)
        grids = psr.grid
        grids_num = size(grids)[1]
        sp = psr.sparks
        spark_num = size(sp)[1]

        potentials = []
        electric_fields = []
        drift_velocities = []
        sparks_velocities = []

        psr.pot_minmax = [1e50, -1e50]

        for ii in 1:grids_num
            gr = grids[ii]
            grid_size = size(gr[1])[1]

            vs = Array{Float64}(undef, grid_size, grid_size)

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
                    if vv < psr.pot_minmax[1]
                        psr.pot_minmax[1] = vv
                    end
                    if vv > psr.pot_minmax[2]
                        psr.pot_minmax[2] = vv
                    end
                    vs[i, j] = vv
                    #println(ii, " ", i, " ", j, " ", vs[i, j])
                end
            end
            push!(potentials, vs)

            # calculate electric field
            ex = Array{Float64}(undef, grid_size, grid_size)
            ey = Array{Float64}(undef, grid_size, grid_size)

            # julia gradient calculation
            grad_v = Functions.numpy_gradient_2d(vs)
            grad_vx = grad_v[1]
            grad_vy = grad_v[2]
            ex = - grad_vx
            ey = - grad_vy

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
            ind = convert(Int, ceil(grid_size / 2)) # works for odd grid sizes
            push!(sparks_velocities, [vdx[ind, ind], vdy[ind, ind]])
            push!(electric_fields, [ex, ey])
            push!(drift_velocities, [vdx, vdy])
        end
        #println(typeof(grad_v2))
        psr.potential = potentials
        psr.electric_field = electric_fields
        psr.drift_velocity = drift_velocities # for full grids
        psr.sparks_velocity = sparks_velocities # center velocity in grids
        if save == true
            push!(psr.sparks_locations, deepcopy(psr.sparks))
        end
    end
    function step(psr; speedup=1)

        sv = psr.sparks_velocity
        for (i,s) in enumerate(psr.sparks)
            s[1] = s[1] + sv[i][1] * speedup
            s[2] = s[2] + sv[i][2] * speedup
            # at the setellar surface
            s[3] = sqrt(psr.r^2 - s[1]^2 - s[2]^2)
        end


    end
"""
    Calculates electric potential, electric field and drift velocity for sparks locations in psr.sparls ..
    """
    function calculate_potential_sparks!(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]
        sp = psr.sparks
        spark_num = size(sp)[1]

        vs = Array{Float64}(undef, grid_size, grid_size)

        for i in 1:grid_size
            for j in 1:grid_size
                vv = 0
                for k in 1:spark_num
                    # FIX: Read coordinates directly
                    sx = sp[k][1]
                    sy = sp[k][2]
                    sz = sp[k][3]
                    
                    dist = norm([gr[1][i], gr[2][j], gr[3][i, j]] - [sx, sy, sz])
                    if dist != 0
                        vv += v(dist)
                    end
                end
                if gr[3][i, j] != 0; vs[i, j] = vv; else; vs[i, j] = 0; end
            end
        end

        # Set 0 potential to NaN
        for i in 1:grid_size
            for j in 1:grid_size
                if vs[i, j] == 0; vs[i, j] = NaN; end
            end
        end

        # Calculate fields (Same logic as before)
        grad_v2 = Functions.numpy_gradient_2d(vs)
        ex = -grad_v2[1]; ey = -grad_v2[2]

        vdx = Array{Float64}(undef, grid_size, grid_size)
        vdy = Array{Float64}(undef, grid_size, grid_size)
        for i in 1:grid_size
            for j in 1:grid_size
                B = Field.bd(gr[1][i], gr[1][j], psr)
                E = [ex[i, j], ey[i, j], 0]
                v = cross(E, B)
                vdx[i, j] = v[1]; vdy[i, j] = v[2]
            end
        end

        psr.potential = vs
        push!(psr.potential_simulation, deepcopy(vs))
        psr.electric_field = [ex, ey]
        psr.drift_velocity = [vdx, vdy]
    end
    """
    Runs sparks simulation, periodiclly saving results (for better accuracy) and performing full grid calculation
    """
    function simulate_sparks!(psr; n_steps=2000, skip_steps=10, speedup=10)
        for i in 1:n_steps
            save = (i % skip_steps == 0)
            # small grids around sparks
            Sparks.create_grids!(psr) 
            # potentials around sparks (drift direction, etc.)
            Sparks.calculate_potentials!(psr; save=save) 
            # one step based on sparks_velocity
            Sparks.step(psr; speedup=speedup) 
            if save == true
                println("\n--- Simulation step $i/$n_steps ---")
                # full grid for potential calculation
                Sparks.create_grid!(psr) 
                # potential at the polar cap
                Sparks.calculate_potential_sparks!(psr) 
            end
        end

    end

    """
    Electric potential [Filaments]

    # Arguments

    -r: distance from the spark forming region
    """
    function v(r; a=1)
        #println("R ", r)
        #println("v ", a * log(r / 10))
        return a * log(r)  # 2017 paper 
        #return a * log(r/20)
        #return a * log(r) * exp(-(r-1)/40) # wykładnicze wygaszanie dla 150
        #return r < 150 ? - log(r) * exp(-(r-1)/40) : 0.0 # zero po 150
        #return cos(π*r/(2*150)) # some tests
        #return a * (r/150)^2 # solid-body like ?
    end
     """
    Saves sparks positions 
    """
    function save_sparks(psr; num=1)
        dir = "output/$num"
        mkpath(dir) # check if exists if not creates one
        sl = psr.sparks_locations
        @save "$dir/sparks_locations.jld2" sl
    end

    """
    Loads sparks positions
    """
    function load_sparks(psr; num=1) 
        @load "output/$num/sparks_locations.jld2" sl
        psr.sparks_locations = sl
    end

end # module