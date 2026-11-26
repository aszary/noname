module Plot
    #using CairoMakie
    using LsqFit
    using GLMakie
    using GeometryBasics
    using Glob
    using Statistics
    using LinearAlgebra
    using FFTW
    
    import ..Functions
    import ..Sparks
    import ..Tools




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
            lines!(ax, xs, ys, zs, color=:blue, linewidth=1,fxaa=true)
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
            number = num - start + 1
        end

        da = data[start:start+number-1, :]
        da = repeat(da, times) # repeat data X times

        # exact pulse numbers for the Y axis
        total_pulses = size(da, 1)
        pulse_indices = collect(start : start+total_pulses-1)

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

        # reverse the Y axis so that pulse number (time) increases top to
        # bottom, the standard convention for subpulse drift diagrams
        p.left.yreversed = true
        p.center.yreversed = true

        # PLOTTING DATA
        lines!(p.left, intensity, pulses, color=:grey, linewidth=0.5)
        #xlims!(left, [0.01, 1.01])
        ylims!(p.left, [pulse_indices[1] - 0.5, pulse_indices[end] + 0.5])

        # pass longitude/pulse axes explicitly so the heatmap lines up with
        # the average profile and pulse-number panels
        heatmap!(p.center, psr.longitudes, pulse_indices, transpose(da), colormap=cmap)

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


    function average_stokes(psr)
        avg_I = vec(mean(psr.pulses, dims=1))
        avg_Q = vec(mean(psr.stokes_q[1:psr.npulse, :], dims=1))
        avg_U = vec(mean(psr.stokes_u[1:psr.npulse, :], dims=1))
        avg_V = vec(mean(psr.stokes_v[1:psr.npulse, :], dims=1))
        avg_L = sqrt.(avg_Q .^ 2 .+ avg_U .^ 2)

        size_inches = (12 / 2.54, 16 / 2.54)
        size_pt = 72 .* size_inches
        fig = Figure(size=size_pt, fontsize=8, figure_padding=(4, 8, 4, 4))

        ax_pa = Axis(fig[1, 1], ylabel="PA [deg]",
                     xticklabelsvisible=false, xminorticksvisible=true, yminorticksvisible=true)
        ax_flux = Axis(fig[2, 1], xlabel=L"Longitude ($^\circ$)", ylabel="Flux Density",
                       xminorticksvisible=true, yminorticksvisible=true)

        scatter!(ax_pa, psr.longitudes, psr.pa, color=:black, markersize=3)

        # RVM overlay + slope comparison at actual inflection point
        α = deg2rad(psr.alpha)
        β = deg2rad(psr.beta)
        φ = deg2rad.(psr.longitudes)
        n = length(psr.longitudes)
        dφ_deg = psr.longitudes[2] - psr.longitudes[1]
        # numerical dPA/dφ at every bin (central difference, skip edges)
        dpa = [(psr.pa[i+1] - psr.pa[i-1]) / (2 * dφ_deg) for i in 2:n-1]
        # inflection point = steepest slope (max |dPA/dφ|)
        infl = argmax(abs.(dpa)) + 1  # +1 because dpa is indexed from bin 2
        pa0 = psr.pa[infl]
        φ0  = φ[infl]
        slope_numerical  = dpa[infl - 1]  # dpa is offset by 1
        rvm_sign = sign(slope_numerical) != sign(sin(α) / sin(β)) ? -1.0 : 1.0
        slope_analytical = rvm_sign * sin(α) / sin(β)
        β_eff = rad2deg(asin(clamp(abs(sin(α) / slope_numerical), 0.0, 1.0)))
        ζ_eff = α + sign(β) * deg2rad(β_eff)
        rvm_raw = rad2deg.(atan.(sin(α) .* sin.(φ .- φ0),
                                 sin(ζ_eff) .* cos(α) .- cos(ζ_eff) .* sin(α) .* cos.(φ .- φ0)))
        # wrap to (-90°, 90°] to match the PA data normalisation
        pa_rvm = rad2deg.(atan.(tan.(deg2rad.(pa0 .+ rvm_sign .* rvm_raw))))
        lines!(ax_pa, psr.longitudes, pa_rvm, color=:orange, linewidth=1.5, label="RVM (β_eff)")
        vlines!(ax_pa, [psr.longitudes[infl]], color=:blue, linewidth=1, linestyle=:dash)
        scatter!(ax_pa, [psr.longitudes[infl]], [pa0], color=:blue, markersize=7, marker=:diamond, label="inflection")
        axislegend(ax_pa, position=:rt, framevisible=false, labelsize=8)

        println("PA slope at inflection point (lon=$(round(psr.longitudes[infl], digits=2))°):")
        println("  analytical  sin(α)/sin(β) = $(round(slope_analytical, digits=4)) deg/deg")
        println("  numerical   dPA/dφ        = $(round(slope_numerical,  digits=4)) deg/deg")
        println("  |difference|              = $(round(abs(slope_numerical - slope_analytical), digits=4)) deg/deg")
        println("  effective β at r_em       = $(round(β_eff, digits=4))°  (input β = $(psr.beta)°)")

        lines!(ax_flux, psr.longitudes, avg_I, color=:black,  linewidth=1.5, label="I")
        lines!(ax_flux, psr.longitudes, avg_L, color=:red,    linewidth=1.5, label="L")
        lines!(ax_flux, psr.longitudes, avg_V, color=:blue,   linewidth=1.5, label="V")

        xlims!(ax_pa,   psr.longitudes[1], psr.longitudes[end])
        xlims!(ax_flux, psr.longitudes[1], psr.longitudes[end])

        axislegend(ax_flux, position=:rt, framevisible=false, labelsize=8)

        rowsize!(fig.layout, 1, Relative(0.3))
        rowgap!(fig.layout, 3)

        display(fig)
    end


    function polarization_vector_study(psr; delay=0.05)
        fig = Figure(size = (1400, 800))
        
        sphere_mesh_left = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 64))
        sphere_mesh_right = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 64))
        
        axis_scale = 6.0
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis) .* axis_scale
        mag_vec_static = Functions.spherical2cartesian(psr.magnetic_axis) .* axis_scale

        mid_idx = div(length(psr.line_of_sight), 2)
        p_los_center = psr.line_of_sight[mid_idx]

        norm_los = norm(p_los_center)
        dir_los = [p_los_center[1]/norm_los, p_los_center[2]/norm_los, p_los_center[3]/norm_los]
        p_obs_right_static = Point3f(dir_los[1] * psr.r * 4.5, dir_los[2] * psr.r * 4.5, dir_los[3] * psr.r * 4.5)

        vec_len_left = psr.r * 10.0  
        vec_len_right = psr.r * 2.0  
        
        # ==========================================
        # PANEL 1 (Top Left): Rotating Frame
        # ==========================================
        ax1 = Axis3(fig[1, 1], title = "Rotating Frame (Telescope Orbits)")
        mesh!(ax1, sphere_mesh_left, color = (:teal, 0.7), transparency = true)
        
        for line in psr.los_lines
            lines!(ax1, line[1], line[2], line[3], color = (:grey, 0.15), linewidth = 1)
        end

        # Fixed arrow sizes for slim proportions
        arrows3d!(ax1, [Point3f(0,0,0)], [Vec3f(rot_vec...)], color = :red, shaftradius=0.025, tipradius=0.075, tiplength=0.25)
        arrows3d!(ax1, [Point3f(0,0,0)], [Vec3f(mag_vec_static...)], color = :blue, shaftradius=0.025, tipradius=0.075, tiplength=0.25)
        
        active_line_obs_left = Observable(Point3f.(psr.los_lines[1][1], psr.los_lines[1][2], psr.los_lines[1][3]))
        los_obs_left = Observable(Point3f(psr.los_lines[1][1][1], psr.los_lines[1][2][1], psr.los_lines[1][3][1]))
        
        pol_pos_left = Observable([Point3f(los_obs_left[])])
        pol_dir_left = Observable([Vec3f(0,0,0)])
        
        lines!(ax1, active_line_obs_left, color = :blue, linewidth = 3)
        scatter!(ax1, los_obs_left, color=:red, marker=:diamond, markersize=12)
        
        # Yellow polarization vector
        arrows3d!(ax1, pol_pos_left, pol_dir_left, color = :yellow, shaftradius=0.025, tipradius=0.075, tiplength=0.25)

        # ==========================================
        # PANEL 2 (Top Right): Observer Frame
        # ==========================================
        ax2 = Axis3(fig[1, 2], aspect = :data, title = "Observer Frame (Star Rotates)")
        mesh!(ax2, sphere_mesh_right, color = (:teal, 0.7), transparency = true)
        
        arrows3d!(ax2, [Point3f(0,0,0)], [Vec3f(rot_vec...)], color = :red, shaftradius=0.01, tipradius=0.03, tiplength=0.1)
        arrows3d!(ax2, [Point3f(0,0,0)], [Vec3f(p_obs_right_static...)], color = :green, shaftradius=0.01, tipradius=0.03, tiplength=0.1)
        scatter!(ax2, [p_obs_right_static], color=:red, marker=:diamond, markersize=12)
        
        mag_pos_right = Observable([Point3f(0,0,0)])
        mag_dir_right = Observable([Vec3f(mag_vec_static...)])
        
        pol_pos_right = Observable([p_obs_right_static])
        pol_dir_right = Observable([Vec3f(0,0,0)])
        
        arrows3d!(ax2, mag_pos_right, mag_dir_right, color = :blue, shaftradius=0.01, tipradius=0.03, tiplength=0.1)
        arrows3d!(ax2, pol_pos_right, pol_dir_right, color = :yellow, shaftradius=0.025, tipradius=0.075, tiplength=0.25)

        lim = psr.r * 6.0
        limits!(ax2, -lim, lim, -lim, lim, -lim, lim)

        # ==========================================
        # PANEL 3 (Bottom): RVM Curve
        # ==========================================
        ax_pa = Axis(fig[2, 1:2], xlabel = "Longitude [deg]", ylabel = "PA [deg]", title = "RVM Polarization Curve")
        scatter!(ax_pa, psr.longitudes, psr.pa, color = (:grey, 0.5), markersize = 4)
        
        current_pa_obs = Observable(Point2f[Point2f(psr.longitudes[1], psr.pa[1])])
        scatter!(ax_pa, current_pa_obs, color = :yellow, markersize = 12, strokewidth = 1, strokecolor = :black)

        rowsize!(fig.layout, 1, Relative(0.65))
        display(fig)

        # ==========================================
        # ANIMATION LOOP
        # ==========================================
        n_bins = length(psr.longitudes)
        
        while events(fig).window_open[]
            for i in 1:n_bins
                if !events(fig).window_open[] break end
                
                phi_rad = deg2rad(psr.longitudes[i])
                pa_rad = deg2rad(psr.pa[i])
                
                # --- Update Left Panel ---
                line_data = psr.los_lines[i]
                active_line_obs_left[] = Point3f.(line_data[1], line_data[2], line_data[3])
                
                p_los_left = Point3f(line_data[1][1], line_data[2][1], line_data[3][1])
                los_obs_left[] = p_los_left
                
                dx_left = vec_len_left * cos(pa_rad)
                dy_left = vec_len_left * sin(pa_rad)
                
                pol_pos_left[] = [p_los_left]
                pol_dir_left[] = [Vec3f(dx_left, dy_left, 0.0)]
                
                # --- Update Right Panel ---
                mag_vec_rotated = Functions.rotate_vector(mag_vec_static, rot_vec, -phi_rad)
                mag_dir_right[] = [Vec3f(mag_vec_rotated...)]
                
                dx_right = vec_len_right * cos(pa_rad)
                dy_right = vec_len_right * sin(pa_rad)
                
                pol_dir_right[] = [Vec3f(dx_right, dy_right, 0.0)]
                
                # --- Update RVM Plot ---
                current_pa_obs[] = [Point2f(psr.longitudes[i], psr.pa[i])]
                
                sleep(delay)
            end
        end
    end





    """
    closed_lines(psr)

    Visualizes the pulsar magnetosphere with:
    - Stellar surface (teal sphere) 
    - Rotation axis (red) and Magnetic axis (blue)
    - Anomalous dipole moments (orange arrows)
    - General/Closed anomaly-aware field lines (blue) 
    - Open field lines (green) 
    - Line-of-sight paths (red) 
    """
    function closed_lines(psr)
        # 1. Draw the stellar surface
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true)

        rot_vec = Functions.spherical2cartesian(psr.rotation_axis) 
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis) 
        

        rot_scaled = rot_vec .* 0.99
        mag_scaled = mag_vec .* 0.99

        arrows3d!(ax, [0, 0], [0, 0], [0, 0],
                [rot_scaled[1], mag_scaled[1]],
                [rot_scaled[2], mag_scaled[2]],
                [rot_scaled[3], mag_scaled[3]],
                color = [:red, :blue])
        #drawing anomalies
        if hasproperty(psr, :nsfield) && hasproperty(psr.nsfield, :anomalies)
            for a in psr.nsfield.anomalies
                pos = Functions.spherical2cartesian([a.r * psr.r, a.theta_r, a.phi_r])
                dir = Functions.spherical2cartesian([a.m * psr.r, a.theta_m, a.phi_m])

                dir_mag = norm(dir)
                if dir_mag > 0
                    dir_scaled = dir ./ dir_mag .* (psr.r * a.m * 10)
                    arrows3d!(ax, [pos[1]], [pos[2]], [pos[3]],
                            [dir_scaled[1]], [dir_scaled[2]], [dir_scaled[3]],
                            color=:orange)
                end
            end
        end

        # 4. Draw closed field lines (blue)
        for line in psr.nsfield.magnetic_lines
            lines!(ax, line[1], line[2], line[3], color=:blue, linewidth=1, fxaa=true )
        end

        # 5. Draw open field lines (green)
        if !isnothing(psr.open_lines)
            for ml in psr.open_lines 
                lines!(ax, ml[1], ml[2], ml[3], color=:green, linewidth=2,fxaa=true) 
            end # <--- End of open lines loop
        end # <--- End of if block

        # 6. Draw line-of-sight magnetic field lines (red)
        if !isnothing(psr.los_lines) 
            for line in psr.los_lines
                lines!(ax, line[1], line[2], line[3], color=:red, linewidth=1.5,fxaa=true) 
                # Mark the specific emission point/cross at the end of the trace
                scatter!(ax, line[1][end], line[2][end], line[3][end], color=:red, marker=:xcross)
            end # <--- End of LOS lines loop
        end # <--- End of if block

        # Adjust camera for a close 3D view of the surface and anomalies
        cam3d!(ax.scene, eyeposition=[psr.r*2.5, psr.r*2.5, psr.r*2.5], 
               lookat=[0, 0, 0], upvector=[0,0,1], center=false)

        display(fig) 
    end


    function lrfs(psr; start=1, number=nothing, darkness=1.0)
        data = psr.pulses
        
        num_total, bins = size(data)
        if isnothing(number)
            number = num_total - start + 1
        end
        
        # Extract the relevant block of single pulses
        da = data[start:start+number-1, :]
        n_pulses, n_bins = size(da)
        
        # Step 1: Calculate the average profile for the bottom panel
        avg_profile = vec(mean(da, dims=1))
        
        # Step 2: Subtract the mean from each phase bin 
        # This removes the dominant DC component (at frequency f=0)
        da_mean_sub = da .- avg_profile'
        
        # Step 3: Fourier Transform along the time (pulse) axis using FFTW
        # We use rfft since the input data is strictly real.
        # Dimension 1 means we perform FFT for each column (longitude bin) separately.
        fft_complex = rfft(da_mean_sub, 1)
        
        # Calculate the power spectrum (squared amplitude)
        lrfs_power = abs.(fft_complex).^2
        
        # Step 4: Integrated fluctuation spectrum (sum along the longitude axis)
        int_spectrum = vec(sum(lrfs_power, dims=2))
        
        # Create the frequency axis (Y-axis) - Nyquist limit is always 0.5 c/P
        n_freq = size(lrfs_power, 1)
        freqs = range(0.0, 0.5, length=n_freq)
        
        # -----------------------------------------------------
        # CALCULATIONS FOR THE TOP PANEL (FLUCTUATION PHASE)
        # -----------------------------------------------------
        # Find the index of the dominant frequency 
        # (ignoring f=0, which is index 1, so we search from index 2 onwards)
        peak_idx = argmax(int_spectrum[2:end]) + 1
        
        # Extract complex numbers for this specific frequency across all longitudes
        peak_complex = fft_complex[peak_idx, :]
        
        # Calculate the phase (angle) in degrees
        # Note: Advanced analysis often uses phase unwrapping algorithms here 
        # to smooth out sudden jumps between 180 and -180 degrees.
        fluctuation_phase = angle.(peak_complex) .* (180 / pi)
        # -----------------------------------------------------
        
        
        # --- PLOTTING WITH MAKIE (4-panel layout) ---
        
        # Size modeled after the typical 8cm x 13cm layout used in publications
        size_inches = (8 / 2.54, 13 / 2.54) 
        dpi = 150
        size_pt = dpi .* size_inches
        
        fig = Figure(size=size_pt, fontsize=8)
        
        # Left Panel: Integrated fluctuation spectrum
        # Placed on the grid (rows 2:5, column 1)
        ax_left = Axis(fig[2:5, 1], ylabel=L"Frequency ($c/P$)", xlabel="Power",
                       xminorticksvisible=true, yminorticksvisible=true)
        lines!(ax_left, int_spectrum, freqs, color=:gray50, linewidth=1.5)
        ylims!(ax_left, 0.0, 0.5)
        ax_left.xreversed = true # Reverse X-axis to "stick" it to the heatmap
        
        # Main Panel: LRFS Heatmap
        # Placed on the grid (rows 2:5, columns 2:4)
        vmax = darkness * maximum(lrfs_power)
        ax_main = Axis(fig[2:5, 2:4], yticklabelsvisible=false,
                       xminorticksvisible=true, yminorticksvisible=true)
        heatmap!(ax_main, psr.longitudes, freqs, transpose(lrfs_power), colormap=:viridis, colorrange=(0, vmax))
        xlims!(ax_main, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_main, 0.0, 0.5)
        
        # Bottom Panel: Average pulsar profile
        # Placed on the grid (row 6, columns 2:4)
        norm_avg_profile = (avg_profile .- minimum(avg_profile)) ./ maximum(avg_profile .- minimum(avg_profile))
        ax_bottom = Axis(fig[6, 2:4], xlabel=L"Longitude ($^\circ$)", ylabel="Intensity",
                         xminorticksvisible=true, yminorticksvisible=true)
        lines!(ax_bottom, psr.longitudes, norm_avg_profile, color=:gray50, linewidth=1.5)
        xlims!(ax_bottom, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_bottom, 0.0, 1.1)

        # Top Panel: Fluctuation phase
        # Placed on the grid (row 1, columns 2:4)
        ax_top = Axis(fig[1, 2:4], ylabel="Phase\n(deg)", xaxisposition=:top,
                      xminorticksvisible=true, yminorticksvisible=true)
        scatter!(ax_top, psr.longitudes, fluctuation_phase, color=:red, markersize=4)
        xlims!(ax_top, psr.longitudes[1], psr.longitudes[end])
        
        # Synchronize axes so zooming/panning works across all relevant panels
        linkyaxes!(ax_main, ax_left)
        linkxaxes!(ax_main, ax_bottom, ax_top)
        
        hideydecorations!(ax_main, ticks=false, minorticks=false, grid=false)
        hidexdecorations!(ax_main, ticks=false, minorticks=false, grid=false)
        hidexdecorations!(ax_top, ticks=false, minorticks=false, grid=false)
   
        hidespines!(ax_main, :l, :b) 
        hidespines!(ax_left, :r)   
        hidespines!(ax_bottom, :t)  
        hidespines!(ax_top, :b)    
        
        # Adjust gaps between panels
        colgap!(fig.layout, 0)
        rowgap!(fig.layout, 0)
        
        display(fig)
        readline(stdin; keep=false)
        GLMakie.closeall()
    end


"""
    two_dfs(psr; mask_radius=2, darkness=0.8, feature_window=5)
    
    Generates a Two-Dimensional Fluctuation Spectrum (2DFS) plot with an isolated 
    feature integration for the bottom panel.
    
    Layout improvements based on TPA conventions:
    - Uses a color scale (:viridis) for better dynamic range visibility.
    - Bottom panel integrates only the dominant spectral feature in 1/P3, 
      rather than the entire noisy axis.
    - Horizontal lines highlight the integration window.
    """
    function two_dfs(psr; mask_radius=2, darkness=0.8, feature_window=5)
        # 1. Retrieve and center data
        data = psr.pulses
        n_pulses, n_bins = size(data)
        
        avg_profile = mean(data, dims=1)
        data_centered = data .- avg_profile
        
        # 2. 2D Discrete Fourier Transform
        fft_2d = fft(data_centered)
        power_2d = abs2.(fftshift(fft_2d))
        
        # 3. Axis calibration
        nu_3_full = fftshift(fftfreq(n_pulses))
        d_lon = psr.longitudes[2] - psr.longitudes[1]
        nu_2_full = fftshift(fftfreq(n_bins)) .* (360.0 / d_lon)
        
        # Keep only the positive half of 1/P3 (0.0 to 0.5 cpp)
        idx_nu3_pos = findall(x -> x >= 0.0, nu_3_full)
        nu_3 = nu_3_full[idx_nu3_pos]
        power_half = power_2d[idx_nu3_pos, :]
        
        # 4. Mask the DC component
        cx_2 = argmin(abs.(nu_2_full))
        
        rx = min(mask_radius, size(power_half, 1) - 1)
        ry = min(mask_radius, cx_2 - 1, length(nu_2_full) - cx_2)
        
        power_half[1:1+rx, cx_2-ry:cx_2+ry] .= 0.0
        
        # --- ADDED CODE: Determine P2 and P3 from the 2DFS maximum ---
        peak_cartesian = argmax(power_half)
        peak_idx_3 = peak_cartesian[1]
        peak_idx_2 = peak_cartesian[2]

        f_3_peak = nu_3[peak_idx_3]
        f_2_peak = nu_2_full[peak_idx_2]

        # Calculate physical values (P_3 in periods, P_2 in degrees)
        P3_val = f_3_peak != 0.0 ? 1.0 / f_3_peak : Inf
        P2_val = f_2_peak != 0.0 ? 360.0 / abs(f_2_peak) : Inf

        println("\n--------------------------------------------------")
        println("Detected main 2DFS peak:")
        println(" Temporal frequency (1/P3): $(round(f_3_peak, digits=4)) cpp => P3 = $(round(P3_val, digits=2)) periods (P_0)")
        println(" Spatial frequency (1/P2): $(round(f_2_peak, digits=4)) cpp => P2 = $(round(P2_val, digits=2)) degrees")
        println(" Drift direction: $(f_2_peak < 0 ? "negative (-)" : "positive (+)")")
        println("--------------------------------------------------\n")
        # -------------------------------------------------------------

        # 5. Horizontally integrated spectrum (left panel)
        int_power_3 = vec(sum(power_half, dims=2))
        
        # --- IMPROVEMENT: Feature detection for bottom plot ---
        # Find the dominant frequency peak in 1/P3 (ignoring the DC mask region)
        search_start = 1 + rx + 1
        if search_start <= length(int_power_3)
            peak_idx = argmax(int_power_3[search_start:end]) + search_start - 1
            
            # Define window around the peak
            win_start = max(search_start, peak_idx - feature_window)
            win_end = min(length(int_power_3), peak_idx + feature_window)
            
            # Vertically integrate ONLY the dominant feature slice
            int_power_2 = vec(sum(power_half[win_start:win_end, :], dims=1))
            
            peak_freq_start = nu_3[win_start]
            peak_freq_end = nu_3[win_end]
        else
            int_power_2 = vec(sum(power_half, dims=1))
            peak_freq_start = 0.0
            peak_freq_end = 0.5
        end
        
        # 6. Plotting
        size_inches = (10 / 2.54, 10 / 2.54)
        dpi = 150
        size_pt = dpi .* size_inches
        
        fig = Figure(size=size_pt, fontsize=10)
        
        # Left Panel
        ax_left = Axis(fig[1:4, 1], ylabel=L"Temporal Frequency $1/P_3$ [cpp]", xlabel="Power",
                       xminorticksvisible=true, yminorticksvisible=true)
        lines!(ax_left, int_power_3, nu_3, color=:gray30, linewidth=1.5)
        ax_left.xreversed = true 
        ylims!(ax_left, 0.0, 0.5)
        
        # Main Panel (Heatmap)
        vmax = darkness * maximum(power_half)
        ax_main = Axis(fig[1:4, 2:5], yticklabelsvisible=false, xticklabelsvisible=false,
                       xminorticksvisible=true, yminorticksvisible=true)
                       
        # Changed to Viridis colormap for better faint feature distinction
        heatmap!(ax_main, nu_2_full, nu_3, transpose(power_half), 
                 colormap=:viridis, colorrange=(0, vmax))
                 
        # --- ADDED CODE: Draw a marker at the maximum intensity location ---
        if f_2_peak != 0.0 && f_3_peak != 0.0
            scatter!(ax_main, [f_2_peak], [f_3_peak], color=:red, marker=:cross, markersize=15)
        end
        # -------------------------------------------------------------------
                 
        # Draw horizontal lines to indicate the integration window
        if peak_freq_start != 0.0
            hlines!(ax_main, [peak_freq_start, peak_freq_end], color=:white, linestyle=:dot, linewidth=1.2)
            hlines!(ax_left, [peak_freq_start, peak_freq_end], color=:red, linestyle=:dot, linewidth=1.2)
        end
        
        ylims!(ax_main, 0.0, 0.5)
        xlims!(ax_main, nu_2_full[1], nu_2_full[end])
        
        # Bottom Panel
        ax_bottom = Axis(fig[5, 2:5], xlabel=L"Spatial Frequency $1/P_2$ [cpp]", ylabel="Power",
                         xminorticksvisible=true, yminorticksvisible=true)
                         
        # Original integrated spectrum (from the feature window)
        lines!(ax_bottom, nu_2_full, int_power_2, color=:black, linewidth=1.5)
        
        # Mirrored spectrum for drift asymmetry
        lines!(ax_bottom, -nu_2_full, int_power_2, color=:red, linestyle=:dash, linewidth=1.2)
        xlims!(ax_bottom, nu_2_full[1], nu_2_full[end])
        
        # Linking and formatting
        linkyaxes!(ax_main, ax_left)
        linkxaxes!(ax_main, ax_bottom)
        
        hidespines!(ax_main, :l, :b)
        hidespines!(ax_left, :r)
        hidespines!(ax_bottom, :t)
        hidexdecorations!(ax_main, ticks=false, grid=false, minorticks=false)
        hideydecorations!(ax_main, ticks=false, grid=false, minorticks=false)
        
        colgap!(fig.layout, 0)
        rowgap!(fig.layout, 0)
        
        display(fig)
        return fig
    end
    function lrfs2(psr; start=1, number=nothing, darkness=1.0,
                  n_peaks=2, min_df=0.02, peak_freqs=nothing,
                  phase_thresh=0.05, flip_phase=false, show_mod_amp=true)
        data = psr.pulses

        num_total, bins = size(data)
        if isnothing(number)
            number = num_total - start + 1
        end

        # Extract the relevant block of single pulses
        da = data[start:start+number-1, :]
        n_pulses, n_bins = size(da)

        # Step 1: Average profile
        avg_profile = vec(mean(da, dims=1))

        # Step 2: Mean subtraction (removes DC at f = 0)
        da_mean_sub = da .- avg_profile'

        # Step 3: 1D real FFT along the pulse-number axis
        fft_complex = rfft(da_mean_sub, 1)
        lrfs_power = abs.(fft_complex).^2

        # Step 4: Integrated fluctuation spectrum and exact frequency axis
        int_spectrum = vec(sum(lrfs_power, dims=2))
        n_freq = size(lrfs_power, 1)
        freqs = collect(rfftfreq(n_pulses))

        # -----------------------------------------------------
        # FIND TOP n_peaks DISTINCT MAXIMA (OR USE peak_freqs)
        # -----------------------------------------------------
        peak_indices = Int[]
        if !isnothing(peak_freqs)
            for f_target in peak_freqs
                idx = argmin(abs.(freqs .- f_target))
                push!(peak_indices, max(2, idx)) # avoid DC bin
            end
        else
            df = n_freq > 1 ? (freqs[2] - freqs[1]) : 1.0
            min_sep_bins = max(2, round(Int, min_df / df))

            # Collect local maxima (skipping DC index 1)
            candidates = Int[]
            for i in 2:(n_freq - 1)
                if int_spectrum[i] > int_spectrum[i-1] && int_spectrum[i] >= int_spectrum[i+1]
                    push!(candidates, i)
                end
            end
            if n_freq >= 2 && int_spectrum[end] > int_spectrum[end-1]
                push!(candidates, n_freq)
            end

            # Sort strongest first and enforce minimum separation
            sort!(candidates, by = i -> int_spectrum[i], rev = true)
            for idx in candidates
                if all(abs(idx - p) >= min_sep_bins for p in peak_indices)
                    push!(peak_indices, idx)
                    length(peak_indices) == n_peaks && break
                end
            end

            # Fallback if no local max found
            if isempty(peak_indices) && n_freq >= 2
                push!(peak_indices, argmax(int_spectrum[2:end]) + 1)
            end
        end

        println("LRFS detected peaks:")
        for (k, idx) in enumerate(peak_indices)
            f_pk = freqs[idx]
            p3_pk = f_pk > 0 ? 1.0 / f_pk : Inf
            println("  F$k: 1/P3 = $(round(f_pk, digits=4)) cpp  (P3 = $(round(p3_pk, digits=2)) P0)")
        end

        # Color/marker styling (matches F1=blue, F2=green, F3=red from Hsu et al. Fig. 4)
        colors     = [:royalblue, :forestgreen, :crimson, :darkorange, :purple]
        markers    = [:circle, :utriangle, :rect, :diamond, :star5]
        linestyles = [:dash, :dashdot, :dot, :dashdotdot, :solid]

        # --- PLOTTING WITH MAKIE (4-panel layout) ---
        size_inches = (9 / 2.54, 14 / 2.54)
        dpi = 150
        size_pt = dpi .* size_inches

        fig = Figure(size=size_pt, fontsize=8)

        # Left Panel: Integrated fluctuation spectrum
        ax_left = Axis(fig[2:5, 1], ylabel=L"Frequency ($c/P$)", xlabel="Power",
                       xminorticksvisible=true, yminorticksvisible=true)
        lines!(ax_left, int_spectrum, freqs, color=:gray50, linewidth=1.5)
        ylims!(ax_left, 0.0, 0.5)
        ax_left.xreversed = true

        # Main Panel: LRFS Heatmap
        vmax = darkness * maximum(lrfs_power)
        ax_main = Axis(fig[2:5, 2:4], yticklabelsvisible=false,
                       xminorticksvisible=true, yminorticksvisible=true)
        heatmap!(ax_main, psr.longitudes, freqs, transpose(lrfs_power),
                 colormap=:viridis, colorrange=(0, vmax))
        xlims!(ax_main, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_main, 0.0, 0.5)

        # Bottom Panel: Average pulsar profile (+ optional modulation amplitudes)
        prof_min = minimum(avg_profile)
        prof_span = maximum(avg_profile) - prof_min
        norm_avg_profile = prof_span > 0 ? (avg_profile .- prof_min) ./ prof_span : zeros(n_bins)

        ax_bottom = Axis(fig[6, 2:4], xlabel=L"Longitude ($^\circ$)", ylabel="Intensity",
                         xminorticksvisible=true, yminorticksvisible=true)
        lines!(ax_bottom, psr.longitudes, norm_avg_profile, color=:black, linewidth=1.5)
        xlims!(ax_bottom, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_bottom, 0.0, 1.15)

        # Top Panel: Fluctuation phase for F1, F2, ...
        ax_top = Axis(fig[1, 2:4], ylabel="Phase\n(deg)", xaxisposition=:top,
                      yticks=[-180, -90, 0, 90, 180],
                      xminorticksvisible=true, yminorticksvisible=true)
        xlims!(ax_top, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_top, -195, 195)

        phase_sign = flip_phase ? -1.0 : 1.0

        for (k, idx) in enumerate(peak_indices)
            col = colors[mod1(k, length(colors))]
            mkr = markers[mod1(k, length(markers))]
            ls  = linestyles[mod1(k, length(linestyles))]
            f_pk = freqs[idx]

            # Fluctuation phase masked by local power
            peak_complex = fft_complex[idx, :]
            fluct_phase = phase_sign .* angle.(peak_complex) .* (180 / pi)

            p_row = lrfs_power[idx, :]
            p_max = maximum(p_row)
            mask = p_max > 0 ? (p_row .>= phase_thresh * p_max) : trues(n_bins)

            lbl = "F$k ($(round(f_pk, digits=3)))"
            scatter!(ax_top, psr.longitudes[mask], fluct_phase[mask],
                     color=col, marker=mkr, markersize=4.5, label=lbl)

            # Matching frequency marker in the left integrated spectrum
            hlines!(ax_left, [f_pk], color=(col, 0.8), linestyle=:dash, linewidth=1.0)

            # Optional longitude-resolved modulation amplitude (Fig. 4 style)
            if show_mod_amp && prof_span > 0
                scale_factor = (idx == n_freq && iseven(n_pulses)) ? (1.0 / n_pulses) : (2.0 / n_pulses)
                mod_amp = abs.(peak_complex) .* scale_factor
                lines!(ax_bottom, psr.longitudes, mod_amp ./ prof_span,
                       color=col, linestyle=ls, linewidth=1.2)
            end
        end

        # Place legend in the empty top-left corner above the integrated spectrum
        Legend(fig[1, 1], ax_top, framevisible=false, labelsize=6.5,
               patchsize=(8, 8), rowgap=1, padding=(2, 2, 2, 2),
               tellwidth=false, tellheight=false)

        # Synchronize axes
        linkyaxes!(ax_main, ax_left)
        linkxaxes!(ax_main, ax_bottom, ax_top)

        hideydecorations!(ax_main, ticks=false, minorticks=false, grid=false)
        hidexdecorations!(ax_main, ticks=false, minorticks=false, grid=false)
        hidexdecorations!(ax_top, ticks=false, minorticks=false, grid=false)

        hidespines!(ax_main, :l, :b)
        hidespines!(ax_left, :r)
        hidespines!(ax_bottom, :t)
        hidespines!(ax_top, :b)

        colgap!(fig.layout, 0)
        rowgap!(fig.layout, 0)

        display(fig)
        readline(stdin; keep=false)
        GLMakie.closeall()
    end

 """
    slrfs(psr; start=1, number=nothing, win=64, step=1, pad_factor=2,
          window=:hann, lon_range=nothing, normalize_windows=false,
          darkness=0.8, cmap=:viridis, show_true_p3=true, yreversed=false)

    Static 4-panel Sliding LRFS plot:
      - Top-left:     Single-pulse stack (longitude vs pulse number)
      - Top-right:    SLRFS spectrogram (1/P3 vs pulse number) + optional true 1/P3 overlay
      - Bottom-left:  Normalized average pulse profile (highlighting lon_range if set)
      - Bottom-right: Time-averaged fluctuation spectrum
    """
    function slrfs(psr; start=1, number=nothing, win=64, step=1, pad_factor=2,
                   window=:hann, lon_range=nothing, normalize_windows=false,
                   darkness=0.8, cmap=:viridis, show_true_p3=true, yreversed=false)
        data = psr.pulses
        num_total, _ = size(data)
        if isnothing(number)
            number = num_total - start + 1
        end

        da = data[start:(start + number - 1), :]
        n_pulses = size(da, 1)
        pulse_indices = collect(start:(start + n_pulses - 1))

        res = compute_slrfs(da, psr.longitudes;
                            win=win, step=step, pad_factor=pad_factor,
                            window=window, lon_range=lon_range,
                            normalize_windows=normalize_windows)

        # Shift window centers into global pulse numbers
        p_centers = res.pulse_centers .+ (start - 1)

        # Bottom-panel reference curves
        avg_profile = vec(mean(da, dims=1))
        pmin, pmax = extrema(avg_profile)
        norm_prof = (pmax > pmin) ? (avg_profile .- pmin) ./ (pmax - pmin) : zeros(length(avg_profile))

        mean_spec = vec(mean(res.power, dims=2))
        smax = maximum(mean_spec)
        norm_spec = smax > 0 ? mean_spec ./ smax : mean_spec

        # --- Figure Layout (2×2 linked grid) ---
        size_inches = (14 / 2.54, 13 / 2.54)
        fig = Figure(size = 150 .* size_inches, fontsize = 9)

        # Top-Left: Pulse stack
        ax_stack = Axis(fig[1:4, 1:2], ylabel = "Pulse number",
                        xticklabelsvisible = false,
                        xminorticksvisible = true, yminorticksvisible = true)
        heatmap!(ax_stack, psr.longitudes, pulse_indices, permutedims(da), colormap = cmap)
        xlims!(ax_stack, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_stack, pulse_indices[1], pulse_indices[end])

        # Top-Right: SLRFS Spectrogram
        ax_slrfs = Axis(fig[1:4, 3:5],
                        title = "Sliding LRFS (win=$(res.win), step=$step, $(window))",
                        xticklabelsvisible = false, yticklabelsvisible = false,
                        xminorticksvisible = true, yminorticksvisible = true)
        vmax = max(darkness * maximum(res.power), 1e-12)
        heatmap!(ax_slrfs, res.freqs, p_centers, res.power,
                 colormap = cmap, colorrange = (0.0, vmax))
        xlims!(ax_slrfs, 0.0, 0.5)
        ylims!(ax_slrfs, pulse_indices[1], pulse_indices[end])

        # Optional overlay of ground-truth simulated 1/P3(t) from psr.p3
        if show_true_p3 && !isnothing(psr.p3) && isa(psr.p3, AbstractVector) && length(psr.p3) >= start + n_pulses - 1
            p3_slice = psr.p3[start:(start + n_pulses - 1)]
            f3_true = [abs(v) > 0 ? abs(1.0 / v) : NaN for v in p3_slice]
            lines!(ax_slrfs, f3_true, Float64.(pulse_indices),
                   color = (:white, 0.9), linestyle = :dash, linewidth = 1.5,
                   label = L"True $1/P_3(t)$")
            axislegend(ax_slrfs, position = :rt, framevisible = false,
                       labelcolor = :white, labelsize = 8)
        end

        # Bottom-Left: Average pulse profile
        ax_prof = Axis(fig[5, 1:2], xlabel = L"Longitude ($^\circ$)", ylabel = "Intensity",
                       xminorticksvisible = true, yminorticksvisible = true)
        lines!(ax_prof, psr.longitudes, norm_prof, color = :black, linewidth = 1.5)
        xlims!(ax_prof, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_prof, 0.0, 1.1)

        # Highlight longitude slice if lon_range was specified
        if !isnothing(lon_range)
            vlines!(ax_stack, [lon_range[1], lon_range[2]], color = (:white, 0.8), linestyle = :dash, linewidth = 1.0)
            vspan!(ax_prof, lon_range[1], lon_range[2], color = (:orange, 0.25))
        end

        # Bottom-Right: Time-averaged fluctuation spectrum
        ax_spec = Axis(fig[5, 3:5], xlabel = L"Frequency $1/P_3$ ($c/P$)", ylabel = "Avg Power",
                       yticklabelsvisible = false,
                       xminorticksvisible = true, yminorticksvisible = true)
        lines!(ax_spec, res.freqs, norm_spec, color = :black, linewidth = 1.5)
        xlims!(ax_spec, 0.0, 0.5)
        ylims!(ax_spec, 0.0, 1.1)

        if yreversed
            ax_stack.yreversed = true
            ax_slrfs.yreversed = true
        end

        linkyaxes!(ax_stack, ax_slrfs)
        linkxaxes!(ax_stack, ax_prof)
        linkxaxes!(ax_slrfs, ax_spec)

        colgap!(fig.layout, 4)
        rowgap!(fig.layout, 2)

        screen = display(fig)
        wait(screen)
        return fig
    end

    """
    slrfs_interactive(psr; start=1, number=nothing, win=64, step=2, pad_factor=2,
                      window=:hann, lon_range=nothing, darkness=0.5,
                      n_peaks=3, min_df=0.02, peak_freqs=nothing,
                      phase_thresh=0.05, flip_phase=false,
                      show_true_p3=true, delay=0.08)

    Interactive GLMakie SLRFS explorer:
      - Left:  SLRFS spectrogram with a live shaded indicator of the active pulse window
      - Right: Live 4-panel LRFS (multi-peak fluctuation phases, 2D power heatmap,
               integrated spectrum, and local vs. global average profile)
      - Bottom: Window slider + Play/Pause animation button
    """
    function slrfs_interactive(psr; start=1, number=nothing, win=64, step=2, pad_factor=2,
                               window=:hann, lon_range=nothing, darkness=0.5,
                               n_peaks=3, min_df=0.02, peak_freqs=nothing,
                               phase_thresh=0.05, flip_phase=false,
                               show_true_p3=true, delay=0.08)
        data = psr.pulses
        num_total, n_bins = size(data)
        if isnothing(number)
            number = num_total - start + 1
        end

        da = data[start:(start + number - 1), :]
        n_pulses = size(da, 1)
        pulse_indices = collect(start:(start + n_pulses - 1))

        res = compute_slrfs(da, psr.longitudes;
                            win=win, step=step, pad_factor=pad_factor,
                            window=window, lon_range=lon_range)

        p_centers = res.pulse_centers .+ (start - 1)
        n_win     = length(res.starts)
        n_freq    = length(res.freqs)
        phase_sign = flip_phase ? -1.0 : 1.0
        n_tracks  = isnothing(peak_freqs) ? n_peaks : length(peak_freqs)

        # Global normalized average profile (static grey reference in bottom-right panel)
        global_prof = vec(mean(da, dims=1))
        gmin, gmax = extrema(global_prof)
        norm_global = (gmax > gmin) ? (global_prof .- gmin) ./ (gmax - gmin) : zeros(n_bins)

        # Helper to compute one window's full 2D LRFS + multi-peak phases
        buf = zeros(Float64, res.n_fft, n_bins)
        function eval_window(k::Int)
            s = res.starts[k]
            block = @view da[s:(s + res.win - 1), :]
            loc_prof = vec(mean(block, dims=1))

            fill!(buf, 0.0)
            @views buf[1:res.win, :] .= (block .- loc_prof') .* res.weights

            ft = rfft(buf, 1)
            pwr = abs2.(ft)
            pwr[1, :] .= 0.0

            # Sum over selected longitude bins for peak finding
            int_spec = vec(sum(@view(pwr[:, res.bin_mask]), dims=2))
            pks = _find_lrfs_peaks(int_spec, res.freqs;
                                   n_peaks=n_peaks, min_df=min_df, peak_freqs=peak_freqs)

            # Normalize local heatmap & integrated spectrum for steady color/axis limits
            p_max = maximum(pwr)
            pwr_norm = p_max > 0 ? permutedims(pwr) ./ p_max : zeros(Float64, n_bins, n_freq)

            s_max = maximum(int_spec)
            spec_norm = s_max > 0 ? int_spec ./ s_max : zeros(Float64, n_freq)

            lmin, lmax = extrema(loc_prof)
            prof_norm = (lmax > lmin) ? (loc_prof .- lmin) ./ (lmax - lmin) : zeros(n_bins)

            # Build per-peak phase Point2f vectors and frequency markers
            phase_pts = Vector{Vector{Point2f}}(undef, n_tracks)
            pk_freqs_out = fill(NaN, n_tracks)
            for j in 1:n_tracks
                if j <= length(pks)
                    idx = pks[j]
                    pk_freqs_out[j] = res.freqs[idx]
                    ph_deg = phase_sign .* angle.(ft[idx, :]) .* (180 / pi)
                    row_p  = pwr[idx, :]
                    rmax   = maximum(row_p)
                    mask   = rmax > 0 ? (row_p .>= phase_thresh * rmax) : falses(n_bins)
                    if any(mask)
                        phase_pts[j] = Point2f.(psr.longitudes[mask], ph_deg[mask])
                    else
                        phase_pts[j] = Point2f[Point2f(NaN, NaN)]
                    end
                else
                    phase_pts[j] = Point2f[Point2f(NaN, NaN)]
                end
            end

            p_lo = Float64(start + s - 1)
            p_hi = Float64(start + s + res.win - 2)
            return pwr_norm, spec_norm, prof_norm, phase_pts, pk_freqs_out, p_lo, p_hi
        end

        # Evaluate initial window (k = 1)
        pwr0, spec0, prof0, pts0, f0, lo0, hi0 = eval_window(1)

        # --- Observables ---
        heatmap_obs = Observable(pwr0)
        spec_obs    = Observable(spec0)
        prof_obs    = Observable(prof0)
        lo_obs      = Observable(lo0)
        hi_obs      = Observable(hi0)
        phase_obs   = [Observable(pts0[j]) for j in 1:n_tracks]
        freq_obs    = [Observable([f0[j]]) for j in 1:n_tracks]

        function format_status(lo, hi, fvals)
            parts = String[]
            for (j, fv) in enumerate(fvals)
                if isfinite(fv) && fv > 0
                    push!(parts, "F$j=$(round(fv, digits=3)) cpp (P3=$(round(1/fv, digits=1)))")
                end
            end
            return "Window: pulses $(Int(lo))–$(Int(hi))   |   " * join(parts, "   ")
        end
        status_obs = Observable(format_status(lo0, hi0, f0))

        # --- Figure Layout ---
        fig = Figure(size = (1200, 720), fontsize = 10)
        gl_left  = fig[1, 1] = GridLayout()
        gl_right = fig[1, 2] = GridLayout()
        gl_ctrl  = fig[2, 1:2] = GridLayout()

        # Left Panel: SLRFS Overview
        ax_slrfs = Axis(gl_left[1, 1],
                        title = "SLRFS Overview (win=$(res.win))",
                        xlabel = L"Frequency $1/P_3$ ($c/P$)",
                        ylabel = "Pulse number",
                        xminorticksvisible = true, yminorticksvisible = true)
        vmax_s = max(darkness * maximum(res.power), 1e-12)
        heatmap!(ax_slrfs, res.freqs, p_centers, res.power,
                 colormap = :viridis, colorrange = (0.0, vmax_s))
        xlims!(ax_slrfs, 0.0, 0.5)
        ylims!(ax_slrfs, pulse_indices[1], pulse_indices[end])

        if show_true_p3 && !isnothing(psr.p3) && isa(psr.p3, AbstractVector) && length(psr.p3) >= start + n_pulses - 1
            p3_slice = psr.p3[start:(start + n_pulses - 1)]
            f3_true = [abs(v) > 0 ? abs(1.0 / v) : NaN for v in p3_slice]
            lines!(ax_slrfs, f3_true, Float64.(pulse_indices),
                   color = (:white, 0.85), linestyle = :dash, linewidth = 1.5)
        end

        # Live window highlight box on SLRFS
        hspan!(ax_slrfs, lo_obs, hi_obs, color = (:white, 0.22))
        hlines!(ax_slrfs, lo_obs, color = :white, linewidth = 1.2)
        hlines!(ax_slrfs, hi_obs, color = :white, linewidth = 1.2)

        # Right Panels: Live 4-Panel LRFS
        ax_top = Axis(gl_right[1, 2:4], ylabel = "Phase (deg)", xaxisposition = :top,
                      yticks = [-180, -90, 0, 90, 180],
                      xminorticksvisible = true, yminorticksvisible = true)
        xlims!(ax_top, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_top, -195, 195)

        ax_int = Axis(gl_right[2:5, 1], ylabel = L"Frequency ($c/P$)", xlabel = "Power",
                      xminorticksvisible = true, yminorticksvisible = true)
        lines!(ax_int, spec_obs, res.freqs, color = :gray40, linewidth = 1.5)
        ax_int.xreversed = true
        xlims!(ax_int, 0.0, 1.08)
        ylims!(ax_int, 0.0, 0.5)

        ax_main = Axis(gl_right[2:5, 2:4], yticklabelsvisible = false,
                       xminorticksvisible = true, yminorticksvisible = true)
        heatmap!(ax_main, psr.longitudes, res.freqs, heatmap_obs,
                 colormap = :viridis, colorrange = (0.0, darkness))
        xlims!(ax_main, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_main, 0.0, 0.5)

        ax_bot = Axis(gl_right[6, 2:4], xlabel = L"Longitude ($^\circ$)", ylabel = "Intensity",
                      xminorticksvisible = true, yminorticksvisible = true)
        lines!(ax_bot, psr.longitudes, norm_global, color = (:gray60, 0.7), linestyle = :dash, linewidth = 1.2)
        lines!(ax_bot, psr.longitudes, prof_obs, color = :black, linewidth = 1.5)
        xlims!(ax_bot, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_bot, 0.0, 1.15)

        colors  = [:royalblue, :forestgreen, :crimson, :darkorange, :purple]
        markers = [:circle, :utriangle, :rect, :diamond, :star5]
        for j in 1:n_tracks
            col = colors[mod1(j, length(colors))]
            mkr = markers[mod1(j, length(markers))]
            scatter!(ax_top, phase_obs[j], color = col, marker = mkr, markersize = 6, label = "F$j")
            hlines!(ax_int, freq_obs[j], color = (col, 0.85), linestyle = :dash, linewidth = 1.2)
        end

        Legend(gl_right[1, 1], ax_top, framevisible = false, labelsize = 9,
               tellwidth = false, tellheight = false)

        linkyaxes!(ax_main, ax_int)
        linkxaxes!(ax_main, ax_bot, ax_top)
        hideydecorations!(ax_main, ticks = false, minorticks = false, grid = false)
        hidexdecorations!(ax_main, ticks = false, minorticks = false, grid = false)
        hidexdecorations!(ax_top, ticks = false, minorticks = false, grid = false)
        hidespines!(ax_main, :l, :b)
        hidespines!(ax_int, :r)
        hidespines!(ax_bot, :t)
        hidespines!(ax_top, :b)
        colgap!(gl_right, 0)
        rowgap!(gl_right, 0)

        colsize!(fig.layout, 1, Relative(0.38))
        rowsize!(fig.layout, 2, Auto())

        # --- Bottom Controls: Play/Pause + Slider + Status ---
        playing = Ref(false)
        btn_play = Button(gl_ctrl[1, 1], label = "Play", width = 75)
        sl = Slider(gl_ctrl[1, 2], range = 1:n_win, startvalue = 1)
        Label(gl_ctrl[2, 1:2], status_obs, fontsize = 10, halign = :left)

        on(sl.value) do k
            pwr_k, spec_k, prof_k, pts_k, f_k, lo_k, hi_k = eval_window(k)
            heatmap_obs[] = pwr_k
            spec_obs[]    = spec_k
            prof_obs[]    = prof_k
            lo_obs[]      = lo_k
            hi_obs[]      = hi_k
            for j in 1:n_tracks
                phase_obs[j][] = pts_k[j]
                freq_obs[j][]  = [f_k[j]]
            end
            status_obs[] = format_status(lo_k, hi_k, f_k)
        end

        on(btn_play.clicks) do _
            playing[] = !playing[]
            btn_play.label = playing[] ? "Pause" : "Play"
        end

        screen = display(fig)

        @async begin
            while events(fig).window_open[]
                if playing[] && n_win > 1
                    next_k = sl.value[] >= n_win ? 1 : sl.value[] + 1
                    set_close_to!(sl, next_k)
                end
                sleep(delay)
            end
        end

        wait(screen)
        return fig
    end
function _align_and_unwrap_phase(lons, raw_deg, mask, ref_bin)
        n = length(lons)
        out = fill(NaN, n)

        # Find contiguous valid longitude runs (e.g. leading & trailing components)
        runs = UnitRange{Int}[]
        i = 1
        while i <= n
            if mask[i]
                j = i
                while j < n && mask[j + 1]
                    j += 1
                end
                push!(runs, i:j)
                i = j + 1
            else
                i += 1
            end
        end
        isempty(runs) && return out

        # Pick reference bin inside mask (prefer profile peak, fallback to nearest valid bin)
        ref_run_idx = findfirst(r -> ref_bin in r, runs)
        anchor_ref = ref_bin
        if isnothing(ref_run_idx)
            valid_bins = findall(mask)
            anchor_ref = valid_bins[argmin(abs.(valid_bins .- ref_bin))]
            ref_run_idx = findfirst(r -> anchor_ref in r, runs)
        end

        # Shift so phase = 0° at the reference peak
        rel = raw_deg .- raw_deg[anchor_ref]

        # Unwrap inside each contiguous component
        for (r_idx, r) in enumerate(runs)
            anchor = (r_idx == ref_run_idx) ? anchor_ref : first(r)
            out[anchor] = (r_idx == ref_run_idx) ? 0.0 : (mod(rel[anchor] + 180.0, 360.0) - 180.0)
            for b in (anchor + 1):last(r)
                d = mod(rel[b] - rel[b - 1] + 180.0, 360.0) - 180.0
                out[b] = out[b - 1] + d
            end
            for b in (anchor - 1):-1:first(r)
                d = mod(rel[b] - rel[b + 1] + 180.0, 360.0) - 180.0
                out[b] = out[b + 1] + d
            end
        end

        # Extrapolate slope from reference component to place other components on the right 360° branch
        rr = runs[ref_run_idx]
        slope = 0.0
        if length(rr) >= 3
            xs = lons[rr]
            ys = out[rr]
            xm, ym = mean(xs), mean(ys)
            denom = sum((xs .- xm).^2)
            slope = denom > 0 ? sum((xs .- xm) .* (ys .- ym)) / denom : 0.0
        end

        for (r_idx, r) in enumerate(runs)
            r_idx == ref_run_idx && continue
            expected = slope * (mean(lons[r]) - lons[anchor_ref])
            shift_k = round((expected - mean(out[r])) / 360.0)
            out[r] .+= 360.0 * shift_k
        end

        return out
    end

    function lrfs3(psr; start=1, number=nothing, darkness=1.0,
               n_peaks=2, min_df=0.02, peak_freqs=nothing,
               phase_thresh=0.05, flip_phase=false, show_mod_amp=true,
               align_unwrap=true, sort_by_freq=true)
        data = psr.pulses

        num_total, bins = size(data)
        if isnothing(number)
            number = num_total - start + 1
        end

        # Extract the relevant block of single pulses
        da = data[start:start+number-1, :]
        n_pulses, n_bins = size(da)

        # Step 1: Average profile
        avg_profile = vec(mean(da, dims=1))

        # Step 2: Mean subtraction (removes DC at f = 0)
        da_mean_sub = da .- avg_profile'

        # Step 3: 1D real FFT along the pulse-number axis
        fft_complex = rfft(da_mean_sub, 1)
        lrfs_power = abs.(fft_complex).^2

        # Step 4: Integrated fluctuation spectrum and exact frequency axis
        int_spectrum = vec(sum(lrfs_power, dims=2))
        n_freq = size(lrfs_power, 1)
        freqs = collect(rfftfreq(n_pulses))

        # -----------------------------------------------------
        # FIND TOP n_peaks DISTINCT MAXIMA (OR USE peak_freqs)
        # -----------------------------------------------------
        peak_indices = Int[]
        if !isnothing(peak_freqs)
            for f_target in peak_freqs
                idx = argmin(abs.(freqs .- f_target))
                push!(peak_indices, max(2, idx)) # avoid DC bin
            end
        else
            df = n_freq > 1 ? (freqs[2] - freqs[1]) : 1.0
            min_sep_bins = max(2, round(Int, min_df / df))

            # Collect local maxima (skipping DC index 1)
            candidates = Int[]
            for i in 2:(n_freq - 1)
                if int_spectrum[i] > int_spectrum[i-1] && int_spectrum[i] >= int_spectrum[i+1]
                    push!(candidates, i)
                end
            end
            if n_freq >= 2 && int_spectrum[end] > int_spectrum[end-1]
                push!(candidates, n_freq)
            end

            # Sort strongest first and enforce minimum separation
            sort!(candidates, by = i -> int_spectrum[i], rev = true)
            for idx in candidates
                if all(abs(idx - p) >= min_sep_bins for p in peak_indices)
                    push!(peak_indices, idx)
                    length(peak_indices) == n_peaks && break
                end
            end

            # Fallback if no local max found
            if isempty(peak_indices) && n_freq >= 2
                push!(peak_indices, argmax(int_spectrum[2:end]) + 1)
            end
        end

        sort_by_freq && sort!(peak_indices)

        println("LRFS detected peaks:")
        for (k, idx) in enumerate(peak_indices)
            f_pk = freqs[idx]
            p3_pk = f_pk > 0 ? 1.0 / f_pk : Inf
            println("  F$k: 1/P3 = $(round(f_pk, digits=4)) cpp  (P3 = $(round(p3_pk, digits=2)) P0)")
        end

        # Color/marker styling (matches F1=blue, F2=green, F3=red from Hsu et al. Fig. 4)
        colors     = [:royalblue, :forestgreen, :crimson, :darkorange, :purple]
        markers    = [:circle, :utriangle, :rect, :diamond, :star5]
        linestyles = [:dash, :dashdot, :dot, :dashdotdot, :solid]

        # --- PLOTTING WITH MAKIE (4-panel layout) ---
        size_inches = (9 / 2.54, 14 / 2.54)
        dpi = 150
        size_pt = dpi .* size_inches

        fig = Figure(size=size_pt, fontsize=8)

        # Left Panel: Integrated fluctuation spectrum
        ax_left = Axis(fig[2:5, 1], ylabel=L"Frequency ($c/P$)", xlabel="Power",
                       xminorticksvisible=true, yminorticksvisible=true)
        int_norm = maximum(int_spectrum) > 0 ? int_spectrum ./ maximum(int_spectrum) : int_spectrum  # <-- REPLACE
        lines!(ax_left, int_norm, freqs, color=:gray50, linewidth=1.5)                             # <-- WITH THESE 2 LINES
        ylims!(ax_left, 0.0, 0.5)
        ax_left.xreversed = true

        # Main Panel: LRFS Heatmap
        vmax = darkness * maximum(lrfs_power)
        ax_main = Axis(fig[2:5, 2:4], yticklabelsvisible=false,
                       xminorticksvisible=true, yminorticksvisible=true)
        heatmap!(ax_main, psr.longitudes, freqs, transpose(lrfs_power),
                 colormap=:viridis, colorrange=(0, vmax))
        xlims!(ax_main, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_main, 0.0, 0.5)

        # Bottom Panel: Average pulsar profile (+ optional modulation amplitudes)
        prof_min = minimum(avg_profile)
        prof_span = maximum(avg_profile) - prof_min
        norm_avg_profile = prof_span > 0 ? (avg_profile .- prof_min) ./ prof_span : zeros(n_bins)

        ax_bottom = Axis(fig[6, 2:4], xlabel=L"Longitude ($^\circ$)", ylabel="Intensity",
                         xminorticksvisible=true, yminorticksvisible=true)
        lines!(ax_bottom, psr.longitudes, norm_avg_profile, color=:black, linewidth=1.5)
        xlims!(ax_bottom, psr.longitudes[1], psr.longitudes[end])
        ylims!(ax_bottom, 0.0, 1.15)

        # Top Panel: Fluctuation phase
        ax_top = Axis(fig[1, 2:4], ylabel="Phase\n(deg)", xaxisposition=:top,
                    xminorticksvisible=true, yminorticksvisible=true)
        xlims!(ax_top, psr.longitudes[1], psr.longitudes[end])
        if !align_unwrap
            ax_top.yticks = [-180, -90, 0, 90, 180]
            ylims!(ax_top, -195, 195)
        end

        # Anchor phase to the trailing peak (matching Basu Fig. 4)
        mid_bin = div(n_bins, 2)
        ref_bin = argmax(avg_profile[mid_bin:end]) + mid_bin - 1
        phase_sign = flip_phase ? -1.0 : 1.0

        for (k, idx) in enumerate(peak_indices)
            col = colors[mod1(k, length(colors))]
            mkr = markers[mod1(k, length(markers))]
            ls  = linestyles[mod1(k, length(linestyles))]
            f_pk = freqs[idx]

            peak_complex = fft_complex[idx, :]
            raw_phase = phase_sign .* angle.(peak_complex) .* (180 / pi)

            p_row = lrfs_power[idx, :]
            p_max = maximum(p_row)
            mask = p_max > 0 ? (p_row .>= phase_thresh * p_max) : trues(n_bins)

            fluct_phase = align_unwrap ?
                _align_and_unwrap_phase(psr.longitudes, raw_phase, mask, ref_bin) :
                raw_phase

            valid = mask .& isfinite.(fluct_phase)
            lbl = "F$k ($(round(f_pk, digits=3)))"
            scatter!(ax_top, psr.longitudes[valid], fluct_phase[valid],
                    color=col, marker=mkr, markersize=4.5, label=lbl)

            hlines!(ax_left, [f_pk], color=(col, 0.8), linestyle=:dash, linewidth=1.0)

            if show_mod_amp && prof_span > 0
                scale_factor = (idx == n_freq && iseven(n_pulses)) ? (1.0 / n_pulses) : (2.0 / n_pulses)
                mod_amp = abs.(peak_complex) .* scale_factor
                lines!(ax_bottom, psr.longitudes, mod_amp ./ prof_span,
                    color=col, linestyle=ls, linewidth=1.2)
            end
        end

        # Place legend in the empty top-left corner above the integrated spectrum
        Legend(fig[1, 1], ax_top, framevisible=false, labelsize=6.5,
               patchsize=(8, 8), rowgap=1, padding=(2, 2, 2, 2),
               tellwidth=false, tellheight=false)

        # Synchronize axes
        linkyaxes!(ax_main, ax_left)
        linkxaxes!(ax_main, ax_bottom, ax_top)

        hideydecorations!(ax_main, ticks=false, minorticks=false, grid=false)
        hidexdecorations!(ax_main, ticks=false, minorticks=false, grid=false)
        hidexdecorations!(ax_top, ticks=false, minorticks=false, grid=false)

        hidespines!(ax_main, :l, :b)
        hidespines!(ax_left, :r)
        hidespines!(ax_bottom, :t)
        hidespines!(ax_top, :b)

        colgap!(fig.layout, 0)
        rowgap!(fig.layout, 0)

        display(fig)
        readline(stdin; keep=false)
        GLMakie.closeall()
    end
 """
    compute_slrfs(data, longitudes; win=64, step=1, pad_factor=2,
                  window=:hann, lon_range=nothing, normalize_windows=false)

    Computes a Sliding Longitude-Resolved Fluctuation Spectrum (SLRFS) by running
    a windowed 1D real FFT along the pulse-number axis and summing power across
    the selected longitude bins.
    """
    function compute_slrfs(data, longitudes; win=64, step=1, pad_factor=2,
                           window=:hann, lon_range=nothing, normalize_windows=false)
        n_pulses, n_bins = size(data)
        win = min(win, n_pulses)
        @assert win >= 4 "Window length `win` must be at least 4 pulses."

        # Optional longitude filter (e.g. isolate component C1 or C2)
        if isnothing(lon_range)
            bin_mask = trues(n_bins)
        else
            bin_mask = (longitudes .>= lon_range[1]) .& (longitudes .<= lon_range[2])
            @assert any(bin_mask) "No longitude bins fall inside lon_range=$(lon_range)."
        end
        bin_idx = findall(bin_mask)

        # Taper weights along pulse axis
        w = if window == :hann
            0.5 .* (1.0 .- cos.(2π .* (0:win-1) ./ (win - 1)))
        elseif window == :hamming
            0.54 .- 0.46 .* cos.(2π .* (0:win-1) ./ (win - 1))
        elseif window == :rect
            ones(Float64, win)
        else
            error("Unknown window `:$(window)`. Use :hann, :hamming, or :rect.")
        end

        n_fft  = win * max(1, Int(pad_factor))
        freqs  = collect(rfftfreq(n_fft))
        n_freq = length(freqs)

        starts = collect(1:step:(n_pulses - win + 1))
        n_win  = length(starts)
        pulse_centers = starts .+ (win - 1) / 2.0

        slrfs_map = zeros(Float64, n_freq, n_win)
        buf = zeros(Float64, n_fft, length(bin_idx))

        for (k, s) in enumerate(starts)
            block = @view data[s:(s + win - 1), bin_idx]
            block_mean = mean(block, dims=1)

            fill!(buf, 0.0)
            @views buf[1:win, :] .= (block .- block_mean) .* w

            ft = rfft(buf, 1)
            pwr = vec(sum(abs2.(ft), dims=2))
            pwr[1] = 0.0 # explicitly zero DC bin

            if normalize_windows
                pmax = maximum(pwr)
                pmax > 0 && (pwr ./= pmax)
            end

            slrfs_map[:, k] = pwr
        end

        return (
            pulse_centers = pulse_centers,
            freqs         = freqs,
            power         = slrfs_map,
            starts        = starts,
            win           = win,
            n_fft         = n_fft,
            weights       = w,
            bin_mask      = bin_mask
        )
    end

end # module end