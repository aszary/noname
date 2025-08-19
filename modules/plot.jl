module Plot
    using GLMakie
    using GeometryBasics
    using LinearAlgebra

    include("functions.jl")
    include("field.jl")


    function pulsar(psr)
        fi = psr.fields

        f = Figure()
        ax = Axis3(f[1, 1], aspect = :equal)
        # Draw a sphere centered at (0,0,0) with radius r
        mesh!(ax, Sphere(Point3f(0, 0, 0), psr.r), color = (:teal, 0.7), transparency = true)
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
        limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        plot_magnetic_field!(ax, psr.fields)
        #plot_magnetic_lines!(ax, psr.fields)
        
        # Draw light cylinder
        #light_cylinder(psr, ax)
        display(f)
    end

    function plot_magnetic_field!(ax, fi; scale=2e3, scatter=true)
        # Extract positions and vectors
        xs = [p[1] for p in fi.locations]
        ys = [p[2] for p in fi.locations]
        zs = [p[3] for p in fi.locations]

        us = [v[1] for v in fi.magnetic_fields]
        vs = [v[2] for v in fi.magnetic_fields]
        ws = [v[3] for v in fi.magnetic_fields]

        # Draw arrows
        for (p, v) in zip(fi.locations, fi.magnetic_fields)
            p1 = Point3f(p...)
            p2 = Point3f((p .+ scale * (v ./ norm(v)))...)
            lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]], [p1[3], p2[3]], color = (:blue, 0.4))
        end

        # Optional scatter colored by |B|
        if scatter
            Bmag = [norm(v) for v in fi.magnetic_fields]
            meshscatter!(ax, xs, ys, zs; markersize=0.5, color=Bmag, colormap=:seismic)
        end

        #for i in 1:length(fi.locations)
        #    x0, y0, z0 = fi.locations[i]
        #    bx, by, bz = fi.magnetic_fields[i]
        #    # scale the vector for visibility
        #    scale = 0.5 * psr.r / maximum(norm.(fi.magnetic_fields))
        #    x1, y1, z1 = x0 + bx*scale, y0 + by*scale, z0 + bz*scale
        #    lines!(ax, [x0, x1], [y0, y1], [z0, z1], color = (:blue, 0.5))
        #end
    end

    
    function plot_magnetic_lines!(ax, fv)
        # fv = psr.fields
        for ml in fv.magnetic_lines
            lines!(ax, ml[1], ml[2], ml[3], color=:blue, linewidth=1.0, transparency=true)
        end
    end




     

    """
        Plots light cylinder
    """
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

end # module end