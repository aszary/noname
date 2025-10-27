module Plot
    using GLMakie
    using GeometryBasics
    using LinearAlgebra
    using Glob

    include("functions.jl")
    include("field.jl")
    include("sparks.jl")


    function pulsar(psr)
        fi = psr.fields

        fig, ax, p = mesh(Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
         ax.scene.plots[1].plots[1].attributes[:aspect] = :data
        
        #f = Figure()
        
        #ax = Axis3(f[1, 1], aspect = :equal, perspectiveness = 0)
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

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)
        #arrows3d!(ax, [10000], [10], [10], [1, 0.5], color = [:red, :blue]) #,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        r = 1.2e4
        theta = 30 # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)

        #draw points
        scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)

        mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        
        #plot_magnetic_field!(ax, psr.fields)
        #plot_magnetic_lines!(ax, psr.fields)
        plot_polar_cap!(ax, psr; color=:red, linewidth=2.0)
        plot_magnetic_lines_from_polar_cap!(ax, psr.fields)
        #plot_sparks2(psr, ax)
        plot_grids(psr, ax)
        plot_sparks(psr, ax)

        cam = cam3d!(ax.scene, eyeposition=[10000, 10000, 20000], lookat =[0, 0, 10000], upvector=[0,0,1], center = false)

        
        # Draw light cylinder
        #light_cylinder(psr, ax)
        display(fig)
    end

    function pulsar2(psr)
        fi = psr.fields


        fig, ax, p = mesh(Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
        
        #f = Figure()
        
        #ax = Axis3(f[1, 1], aspect = :equal, perspectiveness = 0)
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

        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)
        #arrows3d!(ax, [10000], [10], [10], [1, 0.5], color = [:red, :blue]) #,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)

        r = 1.2e4
        theta = 30 # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)

        #draw points
        scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)

        mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        
        plot_magnetic_field!(ax, psr.fields)
        #plot_magnetic_lines!(ax, psr.fields)
        plot_polar_cap!(ax, psr; color=:red, linewidth=2.0)
        plot_magnetic_lines_from_polar_cap!(ax, psr.fields)
        plot_sparks2(psr, ax)
        plot_grids2(psr, ax)
        #plot_sparks(psr, ax)

        cam = cam3d!(ax.scene, eyeposition=[10000, 10000, 20000], lookat =[0, 0, 10000], upvector=[0,0,1], center = false)

        
        # Draw light cylinder
        #light_cylinder(psr, ax)
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

    function steps2DD(psr; delay=0.1)
        # --- Data preparation ---
        gr = psr.grid
        grid_size = size(gr[1])[1]
        mid = Int(grid_size ÷ 2)

        x = Array{Float64}(undef, grid_size * grid_size)
        y = Array{Float64}(undef, grid_size * grid_size)
        v = Array{Float64}(undef, grid_size * grid_size)
        vs = []

        for ii in 1:length(psr.potential_simulation)
            ind = 0
            for i in 1:grid_size
                for j in 1:grid_size
                    ind += 1
                    x[ind] = gr[1][i]
                    y[ind] = gr[2][j]
                    v[ind] = psr.potential_simulation[ii][i, j]
                end
            end
            push!(vs, deepcopy(v))
        end

        # --- PLOTTING ---
        GLMakie.activate!()

        size_inches = (17/2.54, 11/2.54)
        size_pt = 72 .* size_inches

        fig = Figure(size=size_pt, fontsize=8, figure_padding=(5, 5, 5, 5))

        # Create sub-layouts
        grid      = fig[1, 1] = GridLayout()
        grid2     = fig[2, 1] = GridLayout()
        rightgrid = fig[1, 2] = GridLayout()

        # Adjust layout ratios (main square + half-size slices)
        colsize!(fig.layout, 1, Auto(1))
        colsize!(fig.layout, 2, Auto(0.5))
        rowsize!(fig.layout, 1, Auto(1))
        rowsize!(fig.layout, 2, Auto(0.5))

        # --- AXES ---
        ax_main = Axis(grid[1, 1];
            aspect=DataAspect(),
            xlabel="x (m)", ylabel="y (m)",
            xminorticksvisible=true, yminorticksvisible=true
        )

        ax_vert  = Axis(rightgrid[1, 1]; ylabel="Potential (V)")
        ax_horiz = Axis(grid2[1, 1]; xlabel="x (m)")

        # --- PLOTS ---
        lines!(ax_main, psr.pc[1], psr.pc[2], psr.pc[3])

        v_obs = Observable(vs[1])
        heatmap!(ax_main, x, y, v_obs, interpolate=false)

        spark_positions = Observable([Point2f(sp[1], sp[2]) for sp in psr.sparks_locations[1]])
        scatter!(ax_main, spark_positions, marker=:circle, color=(:red, 0),
                markersize=15, strokewidth=2, strokecolor=:red)

        V0 = reshape(vs[1], (grid_size, grid_size))
        v_vert_obs  = Observable(V0[:, mid])
        v_horiz_obs = Observable(V0[mid, :])

        #lines!(ax_vert, gr[2], v_vert_obs, color=:blue)
        # Vertical slice (right) – flipped axes for readability
        #lines!(ax_vert, v_vert_obs, gr[2], color=:blue)
        #ax_vert.xlabel = "Potential (V)"
        #ax_vert.ylabel = "y (m)"

        # Vertical slice (right) – flipped and mirrored
        lines!(ax_vert, v_vert_obs, gr[2], color=:blue)
        ax_vert.xlabel = "Potential (V)"
        ax_vert.ylabel = "y (m)"
        ax_vert.xreversed = true  # mirror horizontally

        
        lines!(ax_horiz, gr[1], v_horiz_obs, color=:blue)

        display(fig)

        # --- ANIMATION ---
        n_steps = min(length(psr.sparks_locations), length(vs))
        i = 1

        

        while true
            println("\n--- Animation step $i/$n_steps ---")

            v_obs[] = vs[i]
            spark_positions[] = [Point2f(sp[1], sp[2]) for sp in psr.sparks_locations[i]]

            V = reshape(vs[i], (grid_size, grid_size))
            v_vert_obs[] = V[:, mid]   # still the same data
            v_horiz_obs[] = V[mid, :]

            sleep(delay)
            i = i == n_steps ? 1 : i + 1
        end
    end




    function plot_magnetic_field!(ax, fi; marker_size=10, colormap=:plasma)
        # Convert positions to flat Float32 arrays
        xs = Float32[p[1] for p in fi.locations]
        ys = Float32[p[2] for p in fi.locations]
        zs = Float32[p[3] for p in fi.locations]

        # Magnitude of B for coloring
        Bmag = Float32[norm(v) for v in fi.magnetic_fields]

        # Plot as scatter (dots)
        scatter!(ax, xs, ys, zs;
                markersize=marker_size,
                color=Bmag,
                colormap=colormap,
                transparency=true,
                alpha=0.5)
    end



    
    function plot_magnetic_lines!(ax, fv)
        # fv = psr.fields
        for ml in fv.magnetic_lines
            lines!(ax, ml[1], ml[2], ml[3], color=:blue, linewidth=1.0, transparency=true)
        end
    end

    function plot_polar_cap!(ax, psr; color=:red, linewidth=2.0)
        # group points by cap (north/south)
        caps = Dict(:north => (Float64[], Float64[], Float64[]),
                    :south => (Float64[], Float64[], Float64[]))

        for (r_pc, θ, φ) in psr.polar_cap
            r = psr.r
            x = r * sin(θ) * cos(φ)
            y = r * sin(θ) * sin(φ)
            z = r * cos(θ)

            if z ≥ 0
                push!(caps[:north][1], x)
                push!(caps[:north][2], y)
                push!(caps[:north][3], z)
            else
                push!(caps[:south][1], x)
                push!(caps[:south][2], y)
                push!(caps[:south][3], z)
            end
        end

        # plot separately
        lines!(ax, caps[:north]..., color=color, linewidth=linewidth)
        lines!(ax, caps[:south]..., color=color, linewidth=linewidth)
    end


   function plot_magnetic_lines_from_polar_cap!(ax, f; color=:blue, linewidth=1.0)
        for ml in f.magnetic_lines_from_cap
            lines!(ax, ml[1], ml[2], ml[3], color=color, linewidth=linewidth)
        end
    end

    function plot_sparks2(psr, ax)
        if psr.sparks != nothing
            for s in psr.sparks
                scatter!(ax, s[1], s[2], s[3], marker=:xcross, color=:red)
            end
        end
    end
    
    function plot_grids2(psr, ax)
        if psr.grid !== nothing
            for gr in psr.grid
                 scatter!(ax, gr[1], gr[2], gr[3], marker=:diamond, color=:black)
            end
        end
    end

    function plot_grids(psr, ax)
        if psr.grid !== nothing
            scatter!(ax, psr.grid[1], psr.grid[2], psr.grid[3], marker=:circle, color=:blue)
        end
    end


    function plot_sparks(psr, ax)  gr = psr.grid
        if psr.sparks !== nothing
            for (i, j) in psr.sparks
                scatter!(ax, gr[1][i], gr[2][j], gr[3][i, j], marker=:xcross, color=:red)
            end
        end
    end
     
    function light_cylinder(psr, ax)
        # Create wireframe cylinder with lines
        θ = range(0, 2π, length=200)
        height = 3 * psr.r
        radius = 2* psr.r
        #radius = 2* psr.r_lc
        
        # Vertical lines
        for t in θ[1:end-1]  # Skip last point to avoid duplication
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
        max_extent = max(height, radius)  # Assuming cylinder height is 2*psr.r
        limits!(ax, -max_extent, max_extent, -max_extent, max_extent, -max_extent, max_extent)

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
        #fig = Figure(resolution = (700, 700)) # changes resolution
        #resize!(fig, (700, 700)) # changes resolution
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

end # module end