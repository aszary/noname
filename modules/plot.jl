module Plot
    #using CairoMakie
    using LsqFit
    using GLMakie
    using GeometryBasics
    using Glob
    using Statistics
    include("functions.jl")
    include("sparks.jl")



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
        for i in eachindex(psr.open_lines)
            for j in eachindex(psr.open_lines[i])
                lines!(ax, psr.open_lines[i][j][1], psr.open_lines[i][j][2], psr.open_lines[i][j][3], color=:green)         
            end
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
                scatter!(ax, gr[1][i], gr[2][j], gr[3][i, j], marker=:xcross, color=:red)
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
        for i in eachindex(psr.open_lines)
            for j in eachindex(psr.open_lines[i])
                lines!(ax, psr.open_lines[i][j][1], psr.open_lines[i][j][2], psr.open_lines[i][j][3], color=:green)         
            end
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
        arrows!(x, y, vx, vy, color=:white)

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

    function signal(psr; delay=0.1)

        # better accuracy for the sphere 
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true) # better camera control (Scene), but zlims does not work

        # plot polar cap
        lines!(ax, psr.pc[1], psr.pc[2], psr.pc[3])

        # line of sight end points
        for line in psr.los_lines
            scatter!(ax, line[1][end], line[2][end], line[3][end], color=:blue, marker=:xcross)
        end

        # spark positions Observable
        spark_positions_obs = Observable(Point3f.(psr.sparks_locations[1]))
        # plotting sparks as spheres
        #meshscatter!(ax, spark_positions_obs, markersize=psr.spark_radius, color=:red)
        scatter!(ax, spark_positions_obs, markersize=psr.spark_radius, color=:red)



        # PULSAR SIGNAL BELOW
        ax_bottom = Axis(fig[2, 1], xlabel="longitude [deg.]", ylabel="Intensity")

        # setting panel sizes
        #colsize!(fig.layout, 2, Relative(0.25))  # lewy panel 25% szerokości
        rowsize!(fig.layout, 2, Relative(0.25))  # dolny panel 25% wysokości
        
        # set y range based on maximum signal
        signal_max = maximum(psr.signal)
        ylims!(ax_bottom, -0.1*signal_max, signal_max * 1.1) 

        signal_obs = Observable(psr.signal[1, :])
        lines!(ax_bottom, signal_obs , color=:black, linewidth=2)




 
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

            signal_obs[] = psr.signal[i, :]

            sleep(delay)  # Opóźnienie między krokami (w sekundach)
            
            i = i+1
            if i == n_steps -1 # infinite loop
                i = 1
            end
        end




    end


end # module end