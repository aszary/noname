module Plot
    using GLMakie
    using GeometryBasics
    using LinearAlgebra
    using Glob
    include("functions.jl")
    include("field.jl")
    include("sparks.jl")
    include("transformations.jl")


    function pulsar(psr)
        fig, ax, p = mesh(Sphere(Point3f(0,0,0), psr.r), color = (:teal, 0.7), transparency = true, shading = false)
        #f = Figure()
        #ax = Axis3(f[1, 1], aspect = :equal)


        # Draw a sphere centered at (0,0,0) with radius r
        #mesh!(ax, Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        #println(psr.rotation_axis)
        #println(psr.magnetic_axis)

        #println(rot_vec)
        #println(mag_vec)
        #drawing arrows
        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)
        #arrows3d!(ax, [10000], [10], [10], [1, 0.5], color = [:red, :blue]) #,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)


        r = 1.2e4
        theta = 30 # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)

        #draw points
        #scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)
        #magnetic_field!(ax, psr)
        magnetic_lines!(ax, psr)
        #plot_grids(psr, ax)
        #plot_sparks(psr, ax)
        plot_sparks2(psr, ax)
        polarcap!(ax, psr)
        bx, by, bz = Field.bc(psr, 500000.0) 
        lines!(ax, bx, by, bz, color=:cyan, linewidth=3, label="Beam")

        # 2. Plot LOS (centered on Rotation Axis)
        # This green circle should now CROSS the cyan circle
        plot_los!(ax, psr, height=500000.0)
        # Calculate the beam circle points at height = 500,000 meters
        # This uses the new Field.bc function which handles rotation automatically
        #bx, by, bz = Field.bc(psr, 500000.0) 
        
        # Plot the ring
        
        # Optional: Draw lines connecting the surface polar cap to the 500km ring
        # to visualize the cone shape. (Assumes psr.pc is defined)
        if psr.pc !== nothing
             # Connect every 5th point to keep the plot clean
             for i in 1:5:length(bx)
                 # Note: This assumes psr.pc points and bx points correspond to the same phi indices.
                 # If psr.pc was generated with a different phi_num, this might look twisted.
                 # Ideally, ensure phi_num in Field.bc matches the one used for psr.pc.
                 
                 # Accessing psr.pc points (handling the Vector of Vectors structure)
                 p_surf_x = psr.pc[1][i]
                 p_surf_y = psr.pc[2][i]
                 p_surf_z = psr.pc[3][i]
                 
                 lines!(ax, [p_surf_x, bx[i]], [p_surf_y, by[i]], [p_surf_z, bz[i]], 
                        color=(:cyan, 0.3), linewidth=1)
             end
        end
        mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        # Draw light cylinder
        #light_cylinder(psr, ax)
        cam = cam3d!(ax.scene, eyeposition=[10000, 10000, 20000], lookat =[0, 0, 10000], upvector=[0,0,1], center = false)
    
        display(fig)
    end
    function small_grid(psr) #pulsar function for smallgrids
        sphere_mesh= GeometryBasics.mesh(Tessellation(Sphere(Point3f(0,0,0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true, shading = false)
        #f = Figure()
        #ax = Axis3(f[1, 1], aspect = :equal)


        # Draw a sphere centered at (0,0,0) with radius r
        #mesh!(ax, Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)

        #println(psr.rotation_axis)
        #println(psr.magnetic_axis)

        #println(rot_vec)
        #println(mag_vec)
        #drawing arrows
        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)
        #arrows3d!(ax, [10000], [10], [10], [1, 0.5], color = [:red, :blue]) #,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)


        r = 1.2e4
        theta = 30 # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)

        #draw points
        #scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)
        #magnetic_field!(ax, psr)
        #magnetic_lines!(ax, psr)
        #plot_grids(psr, ax)
        #plot_sparks(psr, ax)
        plot_sparks2(psr, ax)
        polarcap!(ax, psr)
        mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)
        if psr.grid !== nothing
            for gr in psr.grid
                scatter!(ax, gr[1], gr[2], gr[3], marker=:circle, color=:blue)
            end
        end
        # Draw light cylinder
        #light_cylinder(psr, ax)
        cam = cam3d!(ax.scene, eyeposition=[10000, 10000, 20000], lookat =[0, 0, 10000], upvector=[0,0,1], center = false)
    
        display(fig)
    end


    """
        Plots light cylinder
    
    function light_cylinder(psr, ax)
        # Create wireframe cylinder with lines
        θ = range(0, 2π, length=200)
        height = 3 * psr.r
        radius = 2* psr.r
        #radius = 2* psr.r_lc
        
        # Vertical lines
        for t in θ[1:end-1]  # Skip last point to avoid duplication
            x_line = [radius * cos(t), radius * cos(t)]
            y_line = [radius * sin(t), radius * sin(t)]
            z_line = [-2 * psr.r, height]
            lines!(ax, Point3f.(x_line, y_line, z_line), color = :red, alpha = 0.5, linewidth = 2)
        end
        
        # Horizontal circles at different heights
        z_levels = range(0, height, length=5)
        for z in z_levels
            x_circle = [radius * cos(t) for t in θ]
            y_circle = [radius * sin(t) for t in θ]
            z_circle = fill(z, length(θ))
            lines!(ax, Point3f.(x_circle, y_circle, z_circle), color = :red, alpha = 0.5, linewidth = 2)
        end
        # Set equal limits for all axes
        max_extent = max(height, radius)  # Assuming cylinder height is 2*psr.r
        limits!(ax, -max_extent, max_extent, -max_extent, max_extent, -max_extent, max_extent)

    end
    """
    #calculates magnetic field at positions and colors them by magnitude
    function magnetic_field!(ax, psr)
    fv = psr.fields

    positions = fv.locations
    magnitudes = [norm(B) for B in fv.magnetic_fields]

    xs = [p[1] for p in positions]
    ys = [p[2] for p in positions]
    zs = [p[3] for p in positions]

    scatter!(ax, xs, ys, zs, color=magnitudes, colormap=:viridis, markersize=5)
end
#calulates and plots magnetic field lines
function magnetic_lines!(ax, psr)
    fv = psr.fields
    for line in fv.magnetic_lines
        xs, ys, zs = line[1], line[2], line[3]
        lines!(ax, xs, ys, zs, color=:blue, linewidth=1)
    end
end
function polarcap!(ax, psr; color=:orange, linewidth=2)
    if psr.pc === nothing
        @warn "psr.pc is nothing — polar cap not drawn"
        return
    end

    # Rozpoznaj format: albo [pts...] gdzie pts = [x,y,z], albo [xs, ys, zs]
    pc = psr.pc
    if length(pc) == 3 && isa(pc[1], AbstractVector) && isa(pc[2], AbstractVector) && isa(pc[3], AbstractVector)
        # format Field.pc = [xs, ys, zs]
        xs = pc[1]
        ys = pc[2]
        zs = pc[3]
    else
        # oczekujemy listy punktów [[x,y,z], ...]
        xs = [p[1] for p in pc]
        ys = [p[2] for p in pc]
        zs = [p[3] for p in pc]
    end

    # ensure closed loop
    xs_closed = [xs; xs[1]]
    ys_closed = [ys; ys[1]]
    zs_closed = [zs; zs[1]]

    lines!(ax, xs_closed, ys_closed, zs_closed, color=color, linewidth=linewidth)

    # draw southern cap by z -> -z
    zs_s = -zs
    zs_s_closed = [zs_s; zs_s[1]]
    lines!(ax, xs_closed, ys_closed, zs_s_closed, color=color, linewidth=linewidth)

    # optionally expand limits so cap is visible
    try
        autolimits!(ax)
    catch
        # autolimits! may not be available in some contexts; ignore if fails
    end
end
function plot_grids(psr, ax)
        if psr.grid !== nothing
            scatter!(ax, psr.grid[1], psr.grid[2], psr.grid[3], marker=:circle, color=:blue)
        end
    end


    function plot_sparks(psr, ax)  gr = psr.grid
        if psr.sparks !== nothing
            for (i, j) in psr.sparks
                scatter!(ax, gr[1][i], gr[2][j], gr[3][i, j], marker=:xcross, color=:red)
            end
        end
    end
    function plot_sparks2(psr, ax) 
        if psr.sparks !== nothing
            gr = psr.grid # Retrieve grid in case we are using indices
            
            for s in psr.sparks
                # Case 1: Spark is stored as Grid Indices [i, j]
                if length(s) == 2
                    i = Int(s[1])
                    j = Int(s[2])
                    # Look up coordinates from the grid
                    x = gr[1][i]
                    y = gr[2][j]
                    z = gr[3][i, j]
                    scatter!(ax, x, y, z, marker=:xcross, color=:red, markersize = 20)
                
                # Case 2: Spark is stored as Cartesian Coordinates [x, y, z]
                elseif length(s) == 3
                    scatter!(ax, s[1], s[2], s[3], marker=:xcross, color=:red, markersize = 20)
                end
            end
        end
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
        #filename = get_newfilename("output", "potential2D_", "png")
        #println("Filename: $filename")

        #save(filename, fig)
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
    function steps(psr; n_steps=500, skip_steps=10, speedup=10, delay=0.01)
        sphere_mesh = GeometryBasics.mesh(Tesselation(Sphere(Point3f(0, 0, 0), psr.r), 128))
        fig, ax, p = mesh(sphere_mesh, color = (:teal, 0.7), transparency = true)
        
        # rotation axis
        rot_vec = Functions.spherical2cartesian(psr.rotation_axis)
        # magnetic axis
        mag_vec = Functions.spherical2cartesian(psr.magnetic_axis)
        arrows3d!(ax, [0, 0], [0, 0], [0, 0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])
        
        # draw polar caps
        """
        for pc in psr.pc
            lines!(ax, pc[1], pc[2], pc[3], color=:red, linewidth=1)
        end
        """
        polarcap!(ax, psr)
        spark_plots = []
        # plot sparks
        for sp in psr.sparks
            sp_plot = scatter!(ax, sp[1], sp[2], sp[3], marker=:xcross, color=:red)
            push!(spark_plots, sp_plot)
        end
        
        cam = cam3d!(ax.scene, eyeposition=[629.1865281011553, 799.8239786011065, 10488.971033158241], lookat=[16.260756148965456, 120.13202591109358, 9893.546008743007], upvector=[0.02591359646710347, 0.028736371988157157, 0.9992510727755556], center = false)
        
        display(fig)

        for iteration in 1:n_steps
            println("\n--- Step $iteration/$n_steps ---")
            
            for i in 1:skip_steps
                Sparks.step(psr; speedup=speedup)
                Sparks.create_grids!(psr)
                Sparks.calculate_potentials!(psr)
            end
            
            # update spark positions
            for (i, sp) in enumerate(psr.sparks)
                spark_plots[i][1] = sp[1]
                spark_plots[i][2] = sp[2]
                spark_plots[i][3] = sp[3]
            end
            
            sleep(delay)  # delay between steps
        end
        
        println("\n--- Animation complete ---")
        return fig
    end    


    function steps2D(psr; delay=0.1)

        
        gr = psr.grid
        
        x_range = gr[1]
        y_range = gr[2]
        grid_size_x = length(x_range)
        grid_size_y = length(y_range)
        mid_idx_x = grid_size_x ÷ 2
        mid_idx_y = grid_size_y ÷ 2


        # PLOTTING
        GLMakie.activate!()
        #CairoMakie.activate!()

        # Figure size
        size_inches = (17/2.54, 11/2.54) # 17cm x 11cm
        size_pt = 72 .* size_inches
     
        fig = Figure(size=size_pt, fontsize=8, figure_padding=(1, 2, 3, 3)) 
        
        ax_main = Axis(fig[1, 1]; 
            aspect=DataAspect(), 
            ylabel="y (m)", 
            xminorticksvisible=true, 
            yminorticksvisible=true, 
            xaxisposition=:bottom)

        ax_bottom = Axis(fig[2, 1];)
        ax_right = Axis(fig[1, 2];
            yaxisposition=:right)

        linkxaxes!(ax_main, ax_bottom)
        linkyaxes!(ax_main, ax_right)

        hidexdecorations!(ax_main, grid=false)
        hideydecorations!(ax_right, grid=false)
        
        # Adjust layout sizes
        colsize!(fig.layout, 2, Relative(0.2)) 
        rowsize!(fig.layout, 2, Relative(0.2)) 
        # plot polar cap
        lines!(ax_main, psr.pc[1], psr.pc[2], psr.pc[3])

        v_observable = Observable(psr.potential_simulation[1])
        
        v_slice_x = lift(v_obs -> v_obs[mid_idx_x, :], v_observable)
        v_slice_y = lift(v_obs -> v_obs[:, mid_idx_y], v_observable)

        spark_positions = Observable([Point2f(sp[1], sp[2]) for sp in psr.sparks_locations[1]])

        # === PLOTOWANIE ===
       
        heatmap!(ax_main, x_range, y_range, v_observable, interpolate=false)

        
        scatter!(ax_main, spark_positions, marker=:circle, color=(:red, 0), markersize=15, strokewidth=2, strokecolor=:red)

        
        lines!(ax_bottom, x_range, v_slice_x, color=:blue)
        lines!(ax_right, v_slice_y, y_range, color=:blue)
        
        # autolimits!(ax_bottom, 0.1, 0.1)
        # autolimits!(ax_right, 0.1, 0.1)


        
        button1 = Button(fig[1, 1], label = "Speed up", 
               width = 60, height = 25,
               halign = :right, valign = :top,
               tellwidth = false, tellheight = false)
        on(button1.clicks) do n
            delay = delay / 2
            println("delay: $delay")
        end


        button2 = Button(fig[1, 1], label = "Speed down", 
               width = 60, height = 25,
               halign = :right, valign = :bottom,
               tellwidth = false, tellheight = false)
        on(button2.clicks) do n
            delay = delay * 2
            println("delay: $delay")
        end


        display(fig)

        n_steps = length(psr.sparks_locations)

        i = 1
        while (i < n_steps)
            println("\n--- Animation step $i/$n_steps ---")
            
            # potential update
            v_observable[] = psr.potential_simulation[i] 
            
            # spark positions update
            spark_positions[] = [Point2f(sp[1], sp[2]) for sp in psr.sparks_locations[i]]
           
            sleep(delay)
            
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
"""
    Plots the Line of Sight (LOS) trajectory around the Rotation Axis.
    
    To see the LOS 'cross' the beam, we must plot it correctly:
    - Beam: Circle around Magnetic Axis (Cyan)
    - LOS: Circle around Rotation Axis (Green)
    The intersection of these two non-parallel circles is the pulse.
    """
    function plot_los!(ax, psr; height=500000.0)
        # 1. Calculate Viewing Angle (Zeta)
        # Convert to radians if they are in degrees!
        # Assuming psr.alpha and psr.beta are already in radians based on previous context.
        # If they are in degrees, use: deg2rad(psr.alpha + psr.beta)
        zeta = deg2rad(psr.alpha) + psr.beta
        
        # 2. Define the Rotation Axis (The center of the LOS cone)
        k = Functions.spherical2cartesian(psr.rotation_axis)
        k = k / norm(k) # Normalize rotation axis

        # 3. Construct a basis (u, v) perpendicular to Rotation Axis
        # Arbitrary helper vector (Z-axis)
        aux = [0.0, 0.0, 1.0]
        if abs(dot(k, aux)) > 0.99; aux = [1.0, 0.0, 0.0]; end
        
        u = cross(k, aux)
        u = u / norm(u)
        v = cross(k, u)
        
        # 4. Generate the Circle
        # The LOS vector rotates around k at angle zeta.
        # L(phi) = k*cos(zeta) + (u*cos(phi) + v*sin(phi))*sin(zeta)
        
        phis = range(0, 2pi, length=200)
        xs = Float64[]
        ys = Float64[]
        zs = Float64[]
        
        # Total radius: Star Radius + Altitude
        # MAKE SURE height IS IN METERS (e.g., 500000.0 for 500km)
        total_r = psr.r + height
        
        # Calculate sine and cosine once
        cz = cos(zeta)
        sz = sin(zeta)
        
        for phi in phis
            # Construct unit vector
            dir = k * cz + (u * cos(phi) + v * sin(phi)) * sz
            
            # Scale by radius
            p = dir * total_r
            
            push!(xs, p[1])
            push!(ys, p[2])
            push!(zs, p[3])
        end
        
        # 5. Plot
        lines!(ax, xs, ys, zs, color=:green, linestyle=:dash, linewidth=3, label="LOS Path")
    end
end # module end