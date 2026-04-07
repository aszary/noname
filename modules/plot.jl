module Plot
    #using CairoMakie
    using LsqFit
    using GLMakie
    using GeometryBasics
    using Glob
    using Statistics
    using LinearAlgebra
    include("functions.jl")
    include("sparks.jl")
    include("tools.jl")



    function pulsar(psr)

        #fig, ax, p = mesh(Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true) # better camera control (Scene), but zlims does not work

        # better accuracy for the sphere 
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true) # better camera control (Scene), but zlims does not work

        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        #draw test points
        r = 1.4e4
        theta = 30 # rad2deg(Functions.theta_max(1, psr)) # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)
        scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)

        # line of sight points
        if !isnothing(psr.line_of_sight)
            for l in psr.line_of_sight
                scatter!(ax, l[1], l[2], l[3], markersize=10, color=:red)
            end
        end

        # line of sight lines
        for line in psr.los_lines
            lines!(ax, line[1], line[2], line[3], color=:red, linewidth=1)
        end

        magnetic_field!(ax, psr)

        magnetic_lines!(ax, psr)

        # draw polar caps
        for pc in psr.polar_caps
            lines!(ax, pc[1], pc[2], pc[3], color=:red, linewidth=1)
        end

        # draw open lines
        for ml in psr.open_lines
            lines!(ax, ml[1], ml[2], ml[3], color=:green)
        end


        # plot sparks for random_sparks! function, not very useful
        #=
        if psr.sparks !== nothing
            for sp in psr.sparks
                scatter!(ax, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            end
        end
        =#

        # plot grid
        if psr.grid !== nothing
            scatter!(ax, psr.grid[1], psr.grid[2], psr.grid[3], marker=:diamond, color=:blue)
        end       

        # plot sparks for random_sparks_grid
        gr = psr.grid
        if psr.sparks !== nothing
            for (i, j) in psr.sparks
                scatter!(ax, gr[1][round(Int,i)], gr[2][round(Int,j)], gr[3][round(Int,i), round(Int,j)], marker=:xcross, color=:red)
            end
        end

        # calculate electric potential
        gr = psr.grid
        grid_size = size(gr[1])[1]

        # data for potential plotting
        x = Array{Float64}(undef, grid_size * grid_size)
        y = Array{Float64}(undef, grid_size * grid_size)
        z = Array{Float64}(undef, grid_size * grid_size)
        v = Array{Float64}(undef, grid_size * grid_size)
        ex = Array{Float64}(undef, grid_size * grid_size)
        ey = Array{Float64}(undef, grid_size * grid_size)

        ind = 0
        for i in 1:grid_size
            for j in 1:grid_size
                ind += 1
                x[ind] = gr[1][i]
                y[ind] = gr[2][j]
                z[ind] = gr[3][i,j]
                v[ind] = psr.potential[i, j]
                ex[ind] = psr.electric_field[1][i, j]
                ey[ind] = psr.electric_field[2][i, j]
            end
        end

        # TODO plot potential

        #mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        cam = cam3d!(ax.scene, eyeposition=[10000, 10000, 20000], lookat =[0, 0, 10000], upvector=[0,0,1], center = false)

        # Try accessing the scene's camera directly
        println("Scene camera type: ", typeof(cam))
        println("Scene camera fields: ", fieldnames(typeof(cam)))
        println("Scene camera properties: ", propertynames(cam))

        # Add a button to print camera state
        #button = Button(f[7, 1], label = "Print Camera State")
        button = Button(fig[1, 1], label = "Print", 
               width = 80, height = 25,
               halign = :right, valign = :top,
               tellwidth = false, tellheight = false)


        on(button.clicks) do n
            println("\n--- Camera State (Click $n) ---")
            println("Camera type: ", typeof(cam))
            println("Eye position: ", cam.eyeposition[])
            println("View direction: ", cam.lookat[])
            println("Up vector: ", cam.upvector[])
        end
        
        display(fig)



    end


    function small_grids(psr)
        #fig = Figure()
        #ax = Axis3(f[1, 1], aspect = :equal)
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true) # better camera control (Scene), but zlims does not work

        # Draw a sphere centered at (0,0,0) with radius r
        #mesh!(ax, Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
        
        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        #draw test points
        r = 1.4e4
        theta = 30 #rad2deg(Functions.theta_max(1, psr)) # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)
        scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)

        magnetic_field!(ax, psr)

        magnetic_lines!(ax, psr)

        # draw polar caps
        for pc in psr.polar_caps
            lines!(ax, pc[1], pc[2], pc[3], color=:red, linewidth=1)
        end

        # draw open lines
        for ml in psr.open_lines
            lines!(ax, ml[1], ml[2], ml[3], color=:green)
        end


        # plot sparks for random_sparks! function, not very useful
        if psr.sparks !== nothing
            for sp in psr.sparks
                scatter!(ax, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            end
        end

        # plot grids
        if psr.grid !== nothing
            for gr in psr.grid
                scatter!(ax, gr[1], gr[2], gr[3], marker=:diamond, color=:black)
            #scatter!(ax, psr.grid[1], psr.grid[2], psr.grid[3], marker=:diamond, color=:black)
            end
        end


        #mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        cam = cam3d!(ax.scene, eyeposition=[10000, 10000, 20000], lookat =[0, 0, 10000], upvector=[0,0,1], center = false)

        # Try accessing the scene's camera directly
        println("Scene camera type: ", typeof(cam))
        println("Scene camera fields: ", fieldnames(typeof(cam)))
        println("Scene camera properties: ", propertynames(cam))

        # Add a button to print camera state
        #button = Button(f[7, 1], label = "Print Camera State")
        button = Button(fig[1, 1], label = "Print", 
               width = 80, height = 25,
               halign = :right, valign = :top,
               tellwidth = false, tellheight = false)

        on(button.clicks) do n
            println("\n--- Camera State (Click $n) ---")
            println("Camera type: ", typeof(cam))
            println("Eye position: ", cam.eyeposition[])
            println("View direction: ", cam.lookat[])
            println("Up vector: ", cam.upvector[])
        end
        
        display(fig)

    end


    function magnetic_field!(ax, psr)
        fv = psr.fields

        positions = fv.locations
        magnitudes = [Functions.norm(B) for B in fv.magnetic_fields]

        xs = [p[1] for p in positions]
        ys = [p[2] for p in positions]
        zs = [p[3] for p in positions]

        scatter!(ax, xs, ys, zs, color=magnitudes, colormap=:viridis, markersize=5)
    end 


    function magnetic_lines!(ax, psr)
        fv = psr.fields
        for line in fv.magnetic_lines
            xs, ys, zs = line[1], line[2], line[3]
            lines!(ax, xs, ys, zs, color=:blue, linewidth=1)
        end
    end


    function potential2D(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]

        # data for potential plotting
        x = Array{Float64}(undef, grid_size * grid_size)
        y = Array{Float64}(undef, grid_size * grid_size)
        z = Array{Float64}(undef, grid_size * grid_size)
        v = Array{Float64}(undef, grid_size * grid_size)
        ex = Array{Float64}(undef, grid_size * grid_size)
        ey = Array{Float64}(undef, grid_size * grid_size)
        vx = Array{Float64}(undef, grid_size * grid_size)
        vy = Array{Float64}(undef, grid_size * grid_size)

        ind = 0
        for i in 1:grid_size
            for j in 1:grid_size
                ind += 1
                x[ind] = gr[1][i]
                y[ind] = gr[2][j]
                z[ind] = gr[3][i,j]
                v[ind] = psr.potential[i, j]
                ex[ind] = psr.electric_field[1][i, j]
                ey[ind] = psr.electric_field[2][i, j]
                vx[ind] = psr.drift_velocity[1][i, j]
                vy[ind] = psr.drift_velocity[2][i, j]
            end
        end

        fig, ax1, p = heatmap(x, y, v, interpolate=false) #, colorrange=[-155, -135])
        #hm = meshscatter!(ax1, x, y, ze; markersize=1.25, color=v, transparency=false)
        #arrows!(x, y, ex, ey, color=:white)
        arrows2d!(x, y, vx, vy, color=:white)

        display(fig)
    end


    function potential2Dv2(psr)
        gr = psr.grid
        grid_size = size(gr[1])[1]

        # data for potential plotting
        x = Array{Float64}(undef, grid_size * grid_size)
        y = Array{Float64}(undef, grid_size * grid_size)
        z = Array{Float64}(undef, grid_size * grid_size)
        v = Array{Float64}(undef, grid_size * grid_size)
        ex = Array{Float64}(undef, grid_size * grid_size)
        ey = Array{Float64}(undef, grid_size * grid_size)
        vx = Array{Float64}(undef, grid_size * grid_size)
        vy = Array{Float64}(undef, grid_size * grid_size)

        ind = 0
        for i in 1:grid_size
            for j in 1:grid_size
                ind += 1
                x[ind] = gr[1][i]
                y[ind] = gr[2][j]
                z[ind] = gr[3][i,j]
                v[ind] = psr.potential[i, j]
                ex[ind] = psr.electric_field[1][i, j]
                ey[ind] = psr.electric_field[2][i, j]
                vx[ind] = psr.drift_velocity[1][i, j]
                vy[ind] = psr.drift_velocity[2][i, j]
            end
        end

        #fig = Figure(; resolution=(600, 600))
        #ax = Axis(fig[1, 1]; aspect=(1,1))
        #heatmap!(fig, x, y, v, interpolate=false) #, colorrange=[-155, -135])

        fig, ax1, p = heatmap(x, y, v, interpolate=false) #, colorrange=[-155, -135])
        #resize!(fig, (700, 700)) # changes resolution
        resize_to_layout!(fig)
        #hm = meshscatter!(ax1, x, y, ze; markersize=1.25, color=v, transparency=false)

        #arrows!(x, y, ex, ey, color=:white)
        #arrows!(x, y, vx, vy, color=:white)

        # get last file
        filename = get_newfilename("output", "potential2D_", "png")
        println("Filename: $filename")

        save(filename, fig)
        #save("output/potential2D.svg",fig) # does not work
        #save("output/potential2D.pdf",fig) # does not work
        display(fig)
    end
 

    """
    Returns new file name incremented by +1
    """
    function get_newfilename(dir, filestart, ext="mp4")
        # get last file
        files = glob("$filestart*.$ext", "$dir")
        if length(files) == 0
                return "$dir/$(filestart)1.$ext"
        end
        files = replace.(files, "$dir/$filestart"=>"", ".$ext"=>"")
        nums = parse.(Int, files)
        num = maximum(nums) + 1
        return "$dir/$filestart$num.$ext"
    end



    function steps(psr; n_steps=500, skip_steps=10, speedup=10, delay=0.01)
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true)
        
        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)
        arrows3d!(ax, [0, 0], [0, 0], [0, 0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])
        
        # draw polar caps
        for pc in psr.polar_caps
            lines!(ax, pc[1], pc[2], pc[3], color=:red, linewidth=1)
        end
        
        spark_plots = []
        # plot sparks
        for sp in psr.sparks
            sp_plot = scatter!(ax, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            push!(spark_plots, sp_plot)
        end
        
        cam = cam3d!(ax.scene, eyeposition=[629.1865281011553, 799.8239786011065, 10488.971033158241], lookat=[16.260756148965456, 120.13202591109358, 9893.546008743007], upvector=[0.02591359646710347, 0.028736371988157157, 0.9992510727755556], center = false)
        
        display(fig)
        
        # Automatyczna animacja
        for iteration in 1:n_steps
            println("\n--- Step $iteration/$n_steps ---")
            
            for i in 1:skip_steps
                Sparks.step(psr; speedup=speedup)
                Sparks.create_grids!(psr)
                Sparks.calculate_potentials!(psr)
            end
            
            # Aktualizuj pozycje sparków na wykresie
            for (i, sp) in enumerate(psr.sparks)
                spark_plots[i][1] = sp[1]
                spark_plots[i][2] = sp[2]
                spark_plots[i][3] = sp[3]
            end
            
            sleep(delay)  # Opóźnienie między krokami (w sekundach)
        end
        
        println("\n--- Animation complete ---")
        return fig
    end    


    function steps2D(psr; delay=0.1)

        # calculate electric potential
        gr = psr.grid
        grid_size = size(gr[1])[1]

        # data for potential plotting
        x = Array{Float64}(undef, grid_size * grid_size)
        y = Array{Float64}(undef, grid_size * grid_size)
        z = Array{Float64}(undef, grid_size * grid_size)
        v = Array{Float64}(undef, grid_size * grid_size)
        ex = Array{Float64}(undef, grid_size * grid_size)
        ey = Array{Float64}(undef, grid_size * grid_size)

        vs = [] # potential at the whole polar cap
        vxs = [] # x-cut potential (bottom panel)
        vys = [] # y-cut potential (right panel)

        mid_index = div(grid_size, 2)
        
        for ii in 1:length(psr.potential_simulation)
            ind = 0
            for i in 1:grid_size
                for j in 1:grid_size
                    ind += 1
                    x[ind] = gr[1][i]
                    y[ind] = gr[2][j]
                    z[ind] = gr[3][i,j]
                    #v[ind] = psr.potential[i, j]
                    v[ind] = psr.potential_simulation[ii][i, j]
                    ex[ind] = psr.electric_field[1][i, j]
                    ey[ind] = psr.electric_field[2][i, j]
                end
            end

            push!(vs, deepcopy(v))
    
            # Horizontal cross-section (constant y, middle row)
            vx = psr.potential_simulation[ii][:, mid_index]
            push!(vxs, deepcopy(vx))
    
            # Vertical cross-section (constant x, middle column)
            vy = psr.potential_simulation[ii][mid_index, :]
            push!(vys, deepcopy(vy))

        end



        # PLOTTING
        GLMakie.activate!()
        #CairoMakie.activate!()

        # Figure size
        size_inches = (17/2.54, 11/2.54) # 17cm x 11cm
        size_pt = 72 .* size_inches
        #println(size_pt)
        fig = Figure(size=size_pt, fontsize=8, figure_padding=(1, 2, 3, 3)) # left, right, bottom, top
        ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", xminorticksvisible=true, yminorticksvisible=true, xaxisposition=:bottom)

        # plot polar cap
        lines!(ax, psr.pc[1], psr.pc[2], psr.pc[3])

        v_observable = Observable(vs[1])
        heatmap!(ax, x, y, v_observable, interpolate=false)

        # spark positions Observable
        spark_positions_obs = Observable(psr.sparks_locations[1])
        # circles Obs.
        circles_obs = Observable([Circle(Point2f(p), psr.spark_radius) for p in psr.sparks_locations[1]])
        # plotting sparks
        poly!(ax, circles_obs, color=(:red, 0), strokewidth=2, strokecolor=:red)

        # plot grid
        #scatter!(ax, x, y, marker=:diamond, color=:blue) # simple, works
        #= # use simple scatter..
        grid_size = size(gr[1])[1]
        if psr.grid !== nothing
            for i in 1:grid_size
                for j in 1:grid_size
                    scatter!(ax, psr.grid[1][i], psr.grid[2][j], marker=:diamond, color=:blue)
                end
            end
        end
        =#       


        button1 = Button(fig[2, 2], label = "Speed up", 
               width = 80, height = 25,
               halign = :right, valign = :top,
               tellwidth = false, tellheight = false)
        on(button1.clicks) do n
            delay = delay / 2
            println("delay: $delay")
        end


        button2 = Button(fig[2, 2], label = "Speed down", 
               width = 80, height = 25,
               halign = :right, valign = :bottom,
               tellwidth = false, tellheight = false)
        on(button2.clicks) do n
            delay = delay * 2
            println("delay: $delay")
        end

        ax_right = Axis(fig[1, 2], xlabel="Potential (V)", ylabel="y (m)")
    
        ax_bottom = Axis(fig[2, 1], xlabel="x (m)", ylabel="Potential (V)")

        # setting panel sizes
        colsize!(fig.layout, 2, Relative(0.25))  # lewy panel 25% szerokości
        rowsize!(fig.layout, 2, Relative(0.25))  # dolny panel 25% wysokości

        # cross-section plots
        # Create observables for cross-section plots
        vx_observable = Observable(vxs[1])
        vy_observable = Observable(vys[1])
        # Grid coordinates for x and y axes
        x_coords = gr[1]
        y_coords = gr[2]
        # Plot horizontal cross-section (bottom panel)
        lines!(ax_bottom, x_coords, vx_observable, color=:blue, linewidth=2)
        # Plot vertical cross-section (right panel)
        lines!(ax_right, vy_observable, y_coords, color=:red, linewidth=2)

        # adding solid-body like potential
        # Find center of polar cap (assumed to be potential minimum)
        x_center = mean(psr.pc[1])
        y_center = mean(psr.pc[2])

        # NaN x-data cleaning
        xx = []
        rx = []
        vx = []
        for (i,v) in enumerate(vxs[1])
            if !isnan(v)
                push!(xx, x_coords[i])
                push!(vx, v)
                push!(rx, sqrt((x_coords[i] - x_center)^2))
                #println(v)
            end
        end
        # Funkcja modelu: V(r) = V_offset + V0 * (r/r_scale)^2
        @. model(r, p) = p[1] + p[2] * (r/p[3])^2
        # Usuń punkt centralny (r=0) aby uniknąć problemów numerycznych
        mask_x = rx .> 0.1
        r_x_fit = rx[mask_x]
        vxs_fit = vx[mask_x]
        # Początkowe wartości parametrów [V_offset, V0, r_scale]
        p0_x = [median(vxs_fit), 1.0, 35.0]
        # Dopasowanie dla przekroju x
        fit_x = curve_fit(model, r_x_fit, vxs_fit, p0_x)
        params_x = fit_x.param
        V_theory_x = model(rx, params_x)
        lines!(ax_bottom, xx, V_theory_x, color=:black, linewidth=2, linestyle=:dash, label="Theory (solid-body)")

        # NaN y-data cleaning
        yy = []
        ry = []
        vy = []
        for (i,v) in enumerate(vys[1])
            if !isnan(v)
                push!(yy, y_coords[i])
                push!(vy, v)
                push!(ry, sqrt((y_coords[i] - y_center)^2))
                #println(v)
            end
        end
        # Przygotuj dane dla przekroju y (at x = x_center)
        mask_y = ry .> 0.1
        r_y_fit = ry[mask_y]
        vys_fit = vy[mask_y]
        # Początkowe wartości parametrów
        p0_y = [median(vys_fit), 1.0, 38.0]
        # Dopasowanie dla przekroju y
        fit_y = curve_fit(model, r_y_fit, vys_fit, p0_y)
        params_y = fit_y.param
        V_theory_y = model(ry, params_y)
        lines!(ax_right, V_theory_y, yy, color=:black, linewidth=2, linestyle=:dash, label="Theory (solid-body)")       

        display(fig)

        # plot all steps
        n_steps = length(psr.sparks_locations)

        # Automatyczna animacja
        i = 1
        while (i < n_steps)
            println("\n--- Animation step $i/$n_steps ---")
            # changing potential
            v_observable[] = vs[i]
           
            # update spark positions 
            spark_positions_obs[] = psr.sparks_locations[i]
            circles_obs[] = [Circle(Point2f(p), psr.spark_radius) for p in psr.sparks_locations[i]]

            # update cross-sections
            vx_observable[] = vxs[i]
            vy_observable[] = vys[i]

            sleep(delay)  # Opóźnienie między krokami (w sekundach)
            
            i = i+1
            if i == n_steps -1 # infinite loop
                i = 1
            end
        end

        # save to file # use CairoMakie
        #=
        filename = "output/steps_2D.pdf"
        println(filename)
        save(filename, fig, pt_per_unit = 1)
        =#

        #display(fig)


    end


    function sparks(psr;)

        # better accuracy for the sphere 
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true) # better camera control (Scene), but zlims does not work

        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        # plot polar cap
        lines!(ax, psr.pc[1], psr.pc[2], psr.pc[3])

        # line of sight end points
        for line in psr.los_lines
            scatter!(ax, line[1][end], line[2][end], line[3][end], color=:blue, marker=:xcross)
        end

        if psr.sparks !== nothing
            for sp in psr.sparks
                scatter!(ax, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            end
        end

        cam = cam3d!(ax.scene, eyeposition=[902.365098608735, 388.66374763125975, 10660.389838857573], lookat =[-90.40642962540288, 22.67516168954977, 10092.052717582405], upvector=[0.11471181283596832, 0.042288898277857076, 0.9924982866878566], center = false)

        display(fig)
    end


    function signal(psr; delay=0.1)

        # better accuracy for the sphere 
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true) # better camera control (Scene), but zlims does not work

        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        # plot polar cap
        lines!(ax, psr.pc[1], psr.pc[2], psr.pc[3])

        # line of sight end points
        for line in psr.los_lines
            scatter!(ax, line[1][end], line[2][end], line[3][end], color=:blue, marker=:xcross)
        end

        # spark positions Observable
        spark_positions_obs = Observable(Point3f.(psr.sparks_locations[1]))
        # plotting sparks as spheres if the same radius
        if isnothing(psr.spark_radii)
            meshscatter!(ax, spark_positions_obs, markersize=psr.spark_radius, color=:red)
        else
            spark_radii_obs = Observable(psr.spark_radii[1])
            meshscatter!(ax, spark_positions_obs, markersize=spark_radii_obs, color=:red)
        end
        #scatter!(ax, spark_positions_obs, markersize=psr.spark_radius, color=:red)


        # PULSAR SIGNAL BELOW
        ax_bottom = Axis(fig[2, 1], xlabel="bins", ylabel="Intensity")

        # setting panel sizes
        #colsize!(fig.layout, 2, Relative(0.25))  # lewy panel 25% szerokości
        rowsize!(fig.layout, 2, Relative(0.25))  # dolny panel 25% wysokości
        
        # set y range based on maximum signal
        signal_max = maximum(psr.signal)
        ylims!(ax_bottom, -0.1*signal_max, signal_max * 1.1) 

        signal_obs = Observable(psr.signal[1, :])
        lines!(ax_bottom, signal_obs , color=:black, linewidth=2)


        # Camera aimed at polar cap based on its location and size
        pc_xs, pc_ys, pc_zs = psr.pc[1], psr.pc[2], psr.pc[3]
        pc_center = [mean(pc_xs), mean(pc_ys), mean(pc_zs)]
        pc_size = maximum(sqrt.((pc_xs .- pc_center[1]).^2 .+ (pc_ys .- pc_center[2]).^2 .+ (pc_zs .- pc_center[3]).^2))
        pc_dir = pc_center / norm(pc_center)
        cam_dist = psr.r * 0.1 + 2 * pc_size
        eyepos = pc_center + pc_dir * cam_dist
        view_dir_norm = -pc_dir
        up_candidate = abs(dot(view_dir_norm, [0.0, 0.0, 1.0])) < 0.9 ? [0.0, 0.0, 1.0] : [1.0, 0.0, 0.0]
        upvec = up_candidate - dot(up_candidate, view_dir_norm) * view_dir_norm
        upvec = upvec / norm(upvec)
        cam = cam3d!(ax.scene, eyeposition=eyepos, lookat=pc_center, upvector=upvec, center=false)

        # Try accessing the scene's camera directly
        println("Scene camera type: ", typeof(cam))
        println("Scene camera fields: ", fieldnames(typeof(cam)))
        println("Scene camera properties: ", propertynames(cam))

        # Add a button to print camera state
        #button = Button(f[7, 1], label = "Print Camera State")
        button = Button(fig[1, 1], label = "Print", 
               width = 80, height = 25,
               halign = :right, valign = :top,
               tellwidth = false, tellheight = false)

        on(button.clicks) do n
            println("\n--- Camera State (Click $n) ---")
            println("Camera type: ", typeof(cam))
            println("Eye position: ", cam.eyeposition[])
            println("View direction: ", cam.lookat[])
            println("Up vector: ", cam.upvector[])
        end
       
        # more buttons
        button1 = Button(fig[1, 1], label = "Speed up", 
               width = 80, height = 25,
               halign = :right, valign = :center,
               tellwidth = false, tellheight = false)
        on(button1.clicks) do n
            delay = delay / 2
            println("delay: $delay")
        end

        button2 = Button(fig[1, 1], label = "Speed down", 
               width = 80, height = 25,
               halign = :right, valign = :bottom,
               tellwidth = false, tellheight = false)
        on(button2.clicks) do n
            delay = delay * 2
            println("delay: $delay")
        end

       
        display(fig)

        # plot all steps
        n_steps = length(psr.sparks_locations)

        # Automatyczna animacja
        i = 1
        while (i < n_steps)
            println("\n--- Animation step $i/$n_steps ---")
           
            # update spark positions
            spark_positions_obs[] = psr.sparks_locations[i]
            if !isnothing(psr.spark_radii)
                spark_radii_obs[] = psr.spark_radii[i]
            end

            signal_obs[] = psr.signal[i, :]

            sleep(delay)  # Opóźnienie między krokami (w sekundach)
            
            i = i+1
            if i == n_steps -1 # infinite loop
                i = 1
            end
        end


    end


    function pulses0(psr; start=1, number=100, norm=3.0, name_mod="PSR_NAME")

        data = psr.pulses

        num, bins = size(data)
        if isnothing(number)
            number = num - start  # missing one?
        end

        # Figure size
        size_inches = (8 / 2.54, 11 / 2.54) # 8cm x 11cm
        dpi = 150
        size_pt = dpi .* size_inches

        fig = Figure(size=size_pt, fontsize=8)

        gl = fig[1, 1] = GridLayout()

        ax = Axis(gl[1, 1], xlabel=L"bin number $$", ylabel=L"Pulse number $$", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)
        hideydecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)

        for i in start:1:start+number-1
            da = data[i, :] .* norm .+ i
            #da = da[bin_st:bin_end]
            lines!(ax, psr.longitudes, da, color=:grey, linewidth=0.7)
            #band!(ax, bin_numbers,  ones(length(da)) * i,  da, color=:white) # work on this one day
        end

        ax_profile = Axis(gl[2, 1], xlabel=L"longitude ($^\circ$)", ylabel=L"Intensity $$", xminorticksvisible=true, yminorticksvisible=true)
        mean_profile = vec(mean(data[start:start+number-1, :], dims=1))
        lines!(ax_profile, psr.longitudes, mean_profile, color=:black, linewidth=1.0)
        #lines!(ax_profile, mean_profile, color=:black, linewidth=1.0)
        xlims!(ax_profile, [psr.longitudes[1], psr.longitudes[end]])

        rowsize!(gl, 1, Relative(0.8))
        rowsize!(gl, 2, Relative(0.2))

        display(fig)
        readline(stdin; keep=false)

        #filename = "$outdir/$(name_mod)_single0.pdf"
        #println(filename)
        #save(filename, fig, pt_per_unit=1)

    end


    function pulses1(psr; start=1, number=100, norm=3.0, name_mod="PSR_NAME")

        data = psr.pulses

        num, bins = size(data)
        if isnothing(number)
            number = num - start  # missing one?
        end

        # Figure size
        size_inches = (8 / 2.54, 11 / 2.54) # 8cm x 11cm
        dpi = 150 #72
        size_pt = dpi .* size_inches

        fig = Figure(size=size_pt, fontsize=8)

        gl = fig[1, 1] = GridLayout()

        ax = Axis(gl[1, 1], xlabel="bin number", ylabel=L"Pulse number $$", xminorticksvisible=true, yminorticksvisible=true)
        hidexdecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)
        hideydecorations!(ax, label=false, ticklabels=false, ticks=false, grid=true, minorgrid=false, minorticks=false)

        heatmap!(ax, transpose(data))

        ax_profile = Axis(gl[2, 1], xlabel=L"longitude ($^\circ$)", ylabel=L"Intensity $$", xminorticksvisible=true, yminorticksvisible=true)
        mean_profile = vec(mean(data[start:start+number-1, :], dims=1))
        lines!(ax_profile, psr.longitudes, mean_profile, color=:black, linewidth=1.0)
        #xlims!(ax_profile, [psr.longitudes[1], psr.longitudes[end]])

        rowsize!(gl, 1, Relative(0.8))
        rowsize!(gl, 2, Relative(0.2))

        display(fig)
        readline(stdin; keep=false)

        #filename = "$outdir/$(name_mod)_single0.pdf"
        #println(filename)
        #save(filename, fig, pt_per_unit=1)

    end



    function pulses(psr; start=1, number=100, times=1, cmap="viridis", darkness=0.5, name_mod="PSR_NAME", show_=false)

        data = psr.pulses

        # PREPARE DATA
        num, bins = size(data)
        if number === nothing
            number = num - start  # missing one?
        end

        da = data[start:start+number-1, :]
        da = repeat(da, times) # repeat data X times
        average = Tools.average_profile(da)
        intensity, pulses = Tools.pulses_intensity(da)
        intensity .-= minimum(intensity)
        intensity ./= maximum(intensity)

        pulses .+= start - 1  # julia

        # CREATE FIGURE
        fig, p = triple_panels()
        p.left.ylabel = L"Pulse number $$"
        p.left.xlabel = L"intensity $$"
        p.bottom.xlabel = L"longitude ($^\circ$)"

        # PLOTTING DATA
        lines!(p.left, intensity, pulses, color=:grey, linewidth=0.5)
        #xlims!(left, [0.01, 1.01])
        ylims!(p.left, [pulses[1] - 0.5, pulses[end] + 0.5])

        heatmap!(p.center, transpose(da))

        lines!(p.bottom, psr.longitudes, average, color=:grey, linewidth=0.5)
        xlims!(p.bottom, [psr.longitudes[1], psr.longitudes[end]])

        screen = display(fig)
        readline(stdin; keep=false)
        #resize!(screen, 500, 800)
        
        #filename = "$outdir/$(name_mod)_single.pdf"
        #println(filename)
        #save(filename, fig, pt_per_unit=1)        
        
    end


    struct Panels
        left
        right
        top
        bottom
        center
    end


    function triple_panels()

        # Figure size
        size_inches = (8 / 2.54, 11 / 2.54) # 8cm x 11cm
        dpi = 150 #72
        size_pt = dpi .* size_inches
        #println(size_pt)

        fig = Figure(size=size_pt, fontsize=8)

        left = Axis(fig[1:6, 1], xminorticksvisible=true, yminorticksvisible=true, xticks=[0.5]) # 2:6, 1
        #top = Axis(fig[1, 2:3], xaxisposition=:top, yaxisposition = :right)
        center = Axis(fig[1:6, 2:3]) # 2:6, 2:3
        bottom = Axis(fig[7, 2:3], yaxisposition=:left, xminorticksvisible=true, yminorticksvisible=true, yticklabelsvisible=false)

        left.xreversed = true

        hidedecorations!.(center)
        hidedecorations!.(left, grid=true, ticks=false, ticklabels=false, label=false, minorticks=false)
        hidedecorations!.(bottom, grid=true, ticks=false, ticklabels=false, label=false, minorticks=false)

        colgap!(fig.layout, 0)
        rowgap!(fig.layout, 0)

        return fig, Panels(left, nothing, nothing, bottom, center)

    end



    function anomalies(psr; delay=0.1)

        # better accuracy for the sphere 
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true) # better camera control (Scene), but zlims does not work

        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        # plot polar cap
        #lines!(ax, psr.pc[1], psr.pc[2], psr.pc[3])

        # magnetic field lines from line of sight
        for line in psr.los_lines
            lines!(ax, line[1], line[2], line[3], color=:blue, linewidth=1)
            scatter!(ax, line[1][end], line[2][end], line[3][end], color=:blue, marker=:xcross)
        end

        # line of sight emission points (already in cartesian)
        if !isnothing(psr.line_of_sight)
            scatter!(ax, [p[1] for p in psr.line_of_sight], [p[2] for p in psr.line_of_sight], [p[3] for p in psr.line_of_sight], color=:green, markersize=8)
        end

        # draw open lines
        for ml in psr.open_lines
            lines!(ax, ml[1], ml[2], ml[3], color=:green)
        end

        # anomalies
        for a in psr.nsfield.anomalies
            pos = Functions.spherical2cartesian([a.r * psr.r, a.theta_r, a.phi_r])
            dir = Functions.spherical2cartesian([a.m * psr.r, a.theta_m, a.phi_m])
            arrows3d!(ax, [pos[1]], [pos[2]], [pos[3]], [dir[1]], [dir[2]], [dir[3]], color=:orange)
        end

        cam = cam3d!(ax.scene, eyeposition=[902.365098608735, 388.66374763125975, 10660.389838857573], lookat =[-90.40642962540288, 22.67516168954977, 10092.052717582405], upvector=[0.11471181283596832, 0.042288898277857076, 0.9924982866878566], center = false)

        # Try accessing the scene's camera directly
        println("Scene camera type: ", typeof(cam))
        println("Scene camera fields: ", fieldnames(typeof(cam)))
        println("Scene camera properties: ", propertynames(cam))

        # Add a button to print camera state
        #button = Button(f[7, 1], label = "Print Camera State")
        button = Button(fig[1, 1], label = "Print", 
               width = 80, height = 25,
               halign = :right, valign = :top,
               tellwidth = false, tellheight = false)

        on(button.clicks) do n
            println("\n--- Camera State (Click $n) ---")
            println("Camera type: ", typeof(cam))
            println("Eye position: ", cam.eyeposition[])
            println("View direction: ", cam.lookat[])
            println("Up vector: ", cam.upvector[])
        end
 



        display(fig)

    end




    function anomalies2D(psr)

        fig = Figure(size = (1200, 600))

        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)
        theta_range = range(0, 2π, length=200)
        clip_r = 5 * psr.r

        # legend elements (shared)
        rot_elem = [LineElement(color = :red, linewidth = 2, points = Point2f[(0.0, 0.5), (0.7, 0.5)]),
                    MarkerElement(color = :red, marker = :rtriangle, markersize = 10, points = Point2f[(0.9, 0.5)])]
        mag_elem = [LineElement(color = :blue, linewidth = 2, points = Point2f[(0.0, 0.5), (0.7, 0.5)]),
                    MarkerElement(color = :blue, marker = :rtriangle, markersize = 10, points = Point2f[(0.9, 0.5)])]
        los_elem = [LineElement(color = :red, linewidth = 1),
                    MarkerElement(color = :red, marker = :xcross, markersize = 8)]

        # --- left panel: x-z plane ---
        ax1 = Axis(fig[1, 1], aspect = DataAspect(), xlabel = "x [m]", ylabel = "z [m]", title = "Anomalies (x-z plane)")
        lines!(ax1, psr.r .* cos.(theta_range), psr.r .* sin.(theta_range), color = :black)
        arrows2d!(ax1, [0.0], [0.0], [rot_vec[1]], [rot_vec[3]], color = :red)
        arrows2d!(ax1, [0.0], [0.0], [mag_vec[1]], [mag_vec[3]], color = :blue)
        for a in psr.nsfield.anomalies
            pos = Functions.spherical2cartesian([a.r * psr.r, a.theta_r, a.phi_r])
            dir = Functions.spherical2cartesian([a.m * psr.r, a.theta_m, a.phi_m])
            arrows2d!(ax1, [pos[1]], [pos[3]], [dir[1]], [dir[3]], color = :orange)
            scatter!(ax1, pos[1], pos[3], color = :orange, marker = :circle)
        end
        for line in psr.los_lines
            mask = [sqrt(line[1][i]^2 + line[2][i]^2 + line[3][i]^2) <= clip_r for i in eachindex(line[1])]
            xs = line[1][mask]
            zs = line[3][mask]
            isempty(xs) && continue
            lines!(ax1, xs, zs, color = :red, linewidth = 1)
            scatter!(ax1, line[1][end], line[3][end], color = :red, marker = :xcross)
        end
        #=
        for ml in psr.open_lines
            mask = [sqrt(ml[1][i]^2 + ml[2][i]^2 + ml[3][i]^2) <= clip_r for i in eachindex(ml[1])]
            xs = ml[1][mask]
            zs = ml[3][mask]
            isempty(xs) && continue
            lines!(ax1, xs, zs, color=:green, linewidth=1)
        end
        =#
        xlims!(ax1, -2.5 * psr.r, 2.5 * psr.r)
        ylims!(ax1, -0.5 * psr.r, 3.5 * psr.r)

        # --- right panel: y-z plane ---
        ax2 = Axis(fig[1, 2], aspect = DataAspect(), xlabel = "y [m]", ylabel = "z [m]", title = "Anomalies (y-z plane)")
        lines!(ax2, psr.r .* cos.(theta_range), psr.r .* sin.(theta_range), color = :black)
        arrows2d!(ax2, [0.0], [0.0], [rot_vec[2]], [rot_vec[3]], color = :red)
        arrows2d!(ax2, [0.0], [0.0], [mag_vec[2]], [mag_vec[3]], color = :blue)
        for a in psr.nsfield.anomalies
            pos = Functions.spherical2cartesian([a.r * psr.r, a.theta_r, a.phi_r])
            dir = Functions.spherical2cartesian([a.m * psr.r, a.theta_m, a.phi_m])
            arrows2d!(ax2, [pos[2]], [pos[3]], [dir[2]], [dir[3]], color = :orange)
            scatter!(ax2, pos[2], pos[3], color = :orange, marker = :circle)
        end
        for line in psr.los_lines
            mask = [sqrt(line[1][i]^2 + line[2][i]^2 + line[3][i]^2) <= clip_r for i in eachindex(line[1])]
            ys = line[2][mask]
            zs = line[3][mask]
            isempty(ys) && continue
            lines!(ax2, ys, zs, color = :red, linewidth = 1)
            scatter!(ax2, line[2][end], line[3][end], color = :red, marker = :xcross)
        end
        #=
        for ml in psr.open_lines
            mask = [sqrt(ml[1][i]^2 + ml[2][i]^2 + ml[3][i]^2) <= clip_r for i in eachindex(ml[1])]
            ys = ml[2][mask]
            zs = ml[3][mask]
            isempty(ys) && continue
            lines!(ax2, ys, zs, color=:green, linewidth=1)
        end
        =#
        xlims!(ax2, -2.5 * psr.r, 2.5 * psr.r)
        ylims!(ax2, -0.5 * psr.r, 3.5 * psr.r)

        axislegend(ax2, [rot_elem, mag_elem, los_elem], ["rotation axis", "magnetic axis", "los field lines"])
        display(fig)

    end



    """
    polar_cap2D(psr)

    Top-down 2D view of the polar cap (looking along the magnetic axis / z-axis).
    Plots the footprints of open field lines on the stellar surface, and marks
    the magnetic axis, rotation axis, and anomaly positions.
    """
    function polar_cap2D(psr)

        fig = Figure(size = (700, 700))
        ax = Axis(fig[1, 1], aspect = DataAspect(),
                  xlabel = "x [m]", ylabel = "y [m]",
                  title = "Polar cap – top view (along magnetic axis)")

        # --- open line footprints (done first to get scale) ---
        view_cx, view_cy, arrow_len, ef_a = 0.0, 0.0, 1.5 * psr.r_pc, psr.r_pc
        has_ef = false
        if !isempty(psr.open_lines)
            xs = [ml[1][end] for ml in psr.open_lines]
            ys = [ml[2][end] for ml in psr.open_lines]
            scatter!(ax, xs, ys, color = :green, markersize = 6)
        end

        # --- ellipse fit from psr.ellipse_fit ---
        if !isnothing(psr.ellipse_fit)
            ef = psr.ellipse_fit
            has_ef = true

            # draw ellipse
            cx, cy = ef.center_local
            el_xs = Float64[]
            el_ys = Float64[]
            for t in range(0, 2π, length=300)
                u = cx + ef.a * cos(t) * cos(ef.θ) - ef.b * sin(t) * sin(ef.θ)
                v = cy + ef.a * cos(t) * sin(ef.θ) + ef.b * sin(t) * cos(ef.θ)
                p3 = ef.centroid + u * ef.x_hat + v * ef.y_hat
                push!(el_xs, p3[1])
                push!(el_ys, p3[2])
            end
            lines!(ax, el_xs, el_ys, color = :green, linewidth = 2)
            text!(ax, ef.centroid[1], ef.centroid[2],
                  text = "a=$(round(ef.a, digits=1)) b=$(round(ef.b, digits=1)) θ=$(round(rad2deg(ef.θ), digits=1))°",
                  fontsize = 11, color = :darkgreen, align = (:center, :bottom))

            # view centred on polar cap; arrows scaled to 2× semi-major axis
            view_cx   = ef.centroid[1]
            view_cy   = ef.centroid[2]
            arrow_len = 2.0 * ef.a
            ef_a      = ef.a
        end

        # --- magnetic axis: projects to origin in x-y top-down view ---
        scatter!(ax, [0.0], [0.0], color = :blue, marker = :circle, markersize = 14)

        # --- rotation axis: normalised arrow from origin ---
        rot_vec  = Functions.spherical2cartesian(psr.rotation_axis)
        rot_xy   = [rot_vec[1], rot_vec[2]]
        rot_norm = sqrt(rot_xy[1]^2 + rot_xy[2]^2)
        if rot_norm > 0
            rot_dir = rot_xy ./ rot_norm .* arrow_len
            arrows2d!(ax, [0.0], [0.0], [rot_dir[1]], [rot_dir[2]], color = :red, shaftwidth = 2)
        end

        # --- anomalies: x-y projection of position; moment arrow normalised ---
        for a in psr.nsfield.anomalies
            pos    = Functions.spherical2cartesian([a.r * psr.r, a.theta_r, a.phi_r])
            dir    = Functions.spherical2cartesian([a.m * psr.r, a.theta_m, a.phi_m])
            dir_xy = [dir[1], dir[2]]
            d_norm = sqrt(dir_xy[1]^2 + dir_xy[2]^2)
            if d_norm > 0
                dir_scaled = dir_xy ./ d_norm .* arrow_len
                arrows2d!(ax, [pos[1]], [pos[2]], [dir_scaled[1]], [dir_scaled[2]],
                        color = :orange, shaftwidth = 2)
            end
            scatter!(ax, [pos[1]], [pos[2]], color = :orange, marker = :star5, markersize = 12)
        end

        # --- axis limits: centred on polar cap, 5× semi-major axis padding ---
        pad = has_ef ? 5.0 * ef_a : 3.0 * psr.r_pc
        xlims!(ax, view_cx - pad, view_cx + pad)
        ylims!(ax, view_cy - pad, view_cy + pad)

        # --- manual legend ---
        legend_elems  = Any[MarkerElement(color = :green,  marker = :circle,    markersize = 8),
                            LineElement(  color = :green,  linewidth = 2),
                            MarkerElement(color = :blue,   marker = :circle,    markersize = 12),
                            LineElement(  color = :red,    linewidth = 2)]
        legend_labels = ["open line footprints", "ellipse fit", "magnetic axis", "rotation axis"]
        if !isempty(psr.nsfield.anomalies)
            push!(legend_elems,  MarkerElement(color = :orange, marker = :star5, markersize = 10))
            push!(legend_labels, "anomaly")
        end
        Legend(fig[1, 2], legend_elems, legend_labels)

        display(fig)
    end


end # module end