module Plot
    using GLMakie
    using GeometryBasics
    using LinearAlgebra
    include("functions.jl")



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
        arrows3d!(ax,[0,], [0,0], [0,0], [rot_vec[1], mag_vec[1]], [rot_vec[2], mag_vec[2]], [rot_vec[3], mag_vec[3]], color = [:red, :blue])#,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)
        #arrows3d!(ax, [10000], [10], [10], [1, 0.5], color = [:red, :blue]) #,  shaftradius = 0.01, tipradius = 0.01, tiplength=0.01)


        r = 1.2e4
        theta = 30 # w stopniach
        phi = 0
        p_sp = [r, deg2rad(theta), phi]
        p_car = Functions.spherical2cartesian(p_sp)

        #draw points
        #scatter!(ax, [0, p_car[1]], [0, p_car[2]], [1.2e4, p_car[3]], markersize=10, color=:red)
        #magnetic_field!(ax, psr)
        #magnetic_lines!(ax, psr)
        polarcap!(ax, psr)
        mx = 2e4
        #limits!(ax, -mx, mx, -mx, mx, -mx, mx)

        # Draw light cylinder
        #light_cylinder(psr, ax)

        # ustaw kamerę
        cam3d!(ax.scene,
            eyeposition =  Vec3f(9.78646, 184.23686, 140),
            #lookat      = Vec3f(0, 0, 0),
            upvector    = Vec3f(0, 1, 0)
        )
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

    #calculates northern polar cap (circle)
    cart_points = [Functions.spherical2cartesian(sp) for sp in psr.pc]

    xs = [p[1] for p in cart_points]
    ys = [p[2] for p in cart_points]
    zs = [p[3] for p in cart_points]

    lines!(ax, [xs; xs[1]], [ys; ys[1]], [zs; zs[1]], color=color, linewidth=linewidth)
    #calculates southern polar cap (circle)
    pc_south = [[psr.r, pi - sp[2], sp[3]] for sp in psr.pc]
    cart_points_south = [Functions.spherical2cartesian(sp) for sp in pc_south]

    xs = [p[1] for p in cart_points_south]
    ys = [p[2] for p in cart_points_south]
    zs = [p[3] for p in cart_points_south]

    lines!(ax, [xs; xs[1]], [ys; ys[1]], [zs; zs[1]], color=color, linewidth=2)
end


end # module end