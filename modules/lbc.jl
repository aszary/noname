# Based on Rahul Basu original C code.
module LBC2

using GLMakie


# -----------------------------------------------------------------------------
# coortrans — 2-D coordinate rotation by angle θ (inverse / transpose matrix)
#
#   [ cos θ   sin θ ] [ x ]
#   [-sin θ   cos θ ] [ y ]
#
# This rotates the coordinate axes by +θ, which is equivalent to rotating
# the point itself by -θ.  Used to switch between the polar-cap frame
# (tilted by th_cap) and the display frame.
# -----------------------------------------------------------------------------
function coortrans(x_in::Float64, y_in::Float64, theta::Float64)
    x_out =  x_in * cos(theta) + y_in * sin(theta)
    y_out = -x_in * sin(theta) + y_in * cos(theta)
    return x_out, y_out
end


# -----------------------------------------------------------------------------
# sparkconfig — compute spark positions and rasterize them onto a grid
#
# Arguments
# ---------
#   th_sprk_u  : flat array of angular positions for upper-track sparks
#   th_sprk_d  : flat array of angular positions for lower-track sparks
#                (ring ii occupies indices (ii-1)*trk_max+1 .. (ii-1)*trk_max+trk_max)
#   N_up, N_dn : number of active sparks on the upper / lower track per ring
#   theta_sp   : angular spacing between sparks for each ring [rad]
#   h_sprk     : spark semi-axis (radius) [m]
#   h_drft     : raster grid step size [m]
#   a_cap      : polar cap major semi-axis [m]
#   b_cap      : polar cap minor semi-axis [m]
#   th_cap     : tilt angle of the cap ellipse [rad]
#   co_angl    : co-rotation phase offset [rad]
#   x_cent     : x-coordinate of the cap centre in display frame [m]
#   y_cent     : y-coordinate of the cap centre in display frame [m]
#   N_trk      : number of concentric spark rings
#   trk_max    : stride (maximum sparks per half-ring) in the flat arrays
#
# Returns
# -------
#   spark_x, spark_y : x/y coordinates of each spark centre [m]
#   spark_size        : effective semi-axis of each spark [m]
# -----------------------------------------------------------------------------
function sparkconfig(th_sprk_u, th_sprk_d, N_up, N_dn, theta_sp,
                     h_sprk, h_drft, a_cap, b_cap, th_cap,
                     co_angl, x_cent, y_cent, N_trk, trk_max)

    # Accumulate spark centre positions and their semi-axes
    spark_x    = Float64[]   # x-coordinate of each spark centre
    spark_y    = Float64[]   # y-coordinate of each spark centre
    spark_size = Float64[]   # effective semi-axis of each spark (may shrink near cap edge)

    # Start from the outermost ring and work inward
    a_out = a_cap               # outer major semi-axis of current ring [m]
    a_in  = a_out - 2.0 * h_sprk  # inner major semi-axis of current ring [m]
    b_out = b_cap
    b_in  = b_out - b_cap / a_cap * 2.0 * h_sprk

    for ring in 1:N_trk

        # Mid-line of the current annular ring
        a_trk = 0.5 * (a_out + a_in)   # major semi-axis of track mid-line
        b_trk = 0.5 * (b_out + b_in)   # minor semi-axis of track mid-line

        # Index offset into the flat angle arrays for this ring
        u_off = (ring - 1) * trk_max
        d_off = (ring - 1) * trk_max

        # ------------------------------------------------------------------
        # Upper half-ring: angles decrease from π toward 0 (clockwise drift)
        # ------------------------------------------------------------------
        for jj in 1:N_up[ring]

            ang = th_sprk_u[u_off + jj]   # angular position of this spark [rad]

            # Default spark centre on the track mid-line
            xi = a_trk * cos(ang - th_cap)
            yi = b_trk * sin(ang - th_cap)
            xs, ys = coortrans(xi, yi, -th_cap)
            xs += x_cent;  ys += y_cent
            a_s = h_sprk   # default spark semi-axis

            # --- Spark overlaps the leading edge of the cap (ang ≈ π) ---
            # When the spark is too close to the boundary it is clipped:
            # the semi-axis shrinks and the centre shifts outward.
            if π - co_angl - ang <= theta_sp[ring] / 2
                half_gap  = 0.5 * (π - co_angl - ang) + theta_sp[ring] / 4
                a_s       = a_trk * sin(half_gap)
                trk_a     = a_trk + h_sprk - a_s          # shifted track radius
                xi = trk_a * cos(π - co_angl - half_gap - th_cap)
                yi = trk_a * b_trk / a_trk * sin(π - co_angl - half_gap - th_cap)
                xs, ys = coortrans(xi, yi, -th_cap)
                xs += x_cent;  ys += y_cent
            end

            # --- Spark overlaps the trailing edge of the cap (ang ≈ 0) ---
            if ang <= theta_sp[ring] / 2 - co_angl
                half_gap  = 0.5 * (ang + theta_sp[ring] / 2 + co_angl)
                a_s       = a_trk * sin(half_gap)
                trk_a     = a_trk + h_sprk - a_s
                xi = trk_a * cos(half_gap - th_cap - co_angl)
                yi = trk_a * b_trk / a_trk * sin(half_gap - th_cap - co_angl)
                xs, ys = coortrans(xi, yi, -th_cap)
                xs += x_cent;  ys += y_cent
            end

            push!(spark_x, xs);  push!(spark_y, ys);  push!(spark_size, a_s)
        end

        # ------------------------------------------------------------------
        # Lower half-ring: angles increase from π toward 2π (anti-clockwise)
        # ------------------------------------------------------------------
        for jj in 1:N_dn[ring]

            ang = th_sprk_d[d_off + jj]

            # Default spark centre on the track mid-line
            xi = a_trk * cos(ang - th_cap)
            yi = b_trk * sin(ang - th_cap)
            xs, ys = coortrans(xi, yi, -th_cap)
            xs += x_cent;  ys += y_cent
            a_s = h_sprk

            # --- Spark overlaps the leading edge (ang ≈ π) ---
            if ang - π + co_angl <= theta_sp[ring] / 2
                half_gap = 0.5 * (ang - π + co_angl) + theta_sp[ring] / 4
                a_s      = a_trk * sin(half_gap)
                trk_a    = a_trk + h_sprk - a_s
                xi = trk_a * cos(π - co_angl + half_gap - th_cap)
                yi = trk_a * b_trk / a_trk * sin(π - co_angl + half_gap - th_cap)
                xs, ys = coortrans(xi, yi, -th_cap)
                xs += x_cent;  ys += y_cent
            end

            # --- Spark overlaps the trailing edge (ang ≈ 2π) ---
            if 2π - ang - co_angl <= theta_sp[ring] / 2
                half_gap = 0.5 * (2π - ang - co_angl) + theta_sp[ring] / 4
                a_s      = a_trk * sin(half_gap)
                trk_a    = a_trk + h_sprk - a_s
                xi = trk_a * cos(2π - half_gap - th_cap - co_angl)
                yi = trk_a * b_trk / a_trk * sin(2π - half_gap - th_cap - co_angl)
                xs, ys = coortrans(xi, yi, -th_cap)
                xs += x_cent;  ys += y_cent
            end

            push!(spark_x, xs);  push!(spark_y, ys);  push!(spark_size, a_s)
        end

        # ------------------------------------------------------------------
        # Fill the gap at the θ = π junction between upper and lower tracks.
        # A bridging spark is inserted if the angular gap there is wide enough.
        # ------------------------------------------------------------------
        if th_sprk_d[d_off + 1] - th_sprk_u[u_off + 1] >= theta_sp[ring]
            gap_half = 0.5 * (th_sprk_d[d_off + 1] - th_sprk_u[u_off + 1])
            a_s      = 0.5 * (2 * a_trk * sin(gap_half) - 2.0 * h_sprk)
            e_trk    = a_trk + h_sprk - a_s
            xi = e_trk * cos(π - th_cap - co_angl)
            yi = e_trk * b_trk / a_trk * sin(π - th_cap - co_angl)
            xs, ys = coortrans(xi, yi, -th_cap)
            push!(spark_x, xs + x_cent);  push!(spark_y, ys + y_cent);  push!(spark_size, a_s)
        end

        # ------------------------------------------------------------------
        # Fill the gap at the θ = 0 / 2π junction (end of the ring).
        # ------------------------------------------------------------------
        last_ang_upper = th_sprk_u[u_off + N_up[ring]]          # last upper-track angle
        last_ang_lower = 2π - th_sprk_d[d_off + N_dn[ring]]     # mirrored last lower angle
        gap_sum        = last_ang_lower + last_ang_upper

        if theta_sp[ring] <= gap_sum <= 2 * theta_sp[ring]
            a_s   = 0.5 * (2 * a_trk * sin(0.5 * gap_sum) - 2.0 * h_sprk)
            e_trk = a_trk + h_sprk - a_s
            xi = e_trk * cos(2π - th_cap - co_angl)
            yi = e_trk * b_trk / a_trk * sin(2π - th_cap - co_angl)
            xs, ys = coortrans(xi, yi, -th_cap)
            push!(spark_x, xs + x_cent);  push!(spark_y, ys + y_cent);  push!(spark_size, a_s)
        end

        # Step inward to the next ring
        a_out -= 2.0 * h_sprk;  a_in -= 2.0 * h_sprk
        b_out -= b_cap / a_cap * 2.0 * h_sprk
        b_in  -= b_cap / a_cap * 2.0 * h_sprk
    end # ring end

    # ------------------------------------------------------------------
    # Central (core) spark that fills the innermost region of the cap.
    # Its size is whatever remains after all rings have been laid down.
    # ------------------------------------------------------------------
    push!(spark_x, x_cent);  push!(spark_y, y_cent)
    push!(spark_size, a_cap - N_trk * 2.0 * h_sprk - 1.5 * h_drft)

    return spark_x, spark_y, spark_size
end




function animate(;ntime=200, th_cap=30.0, a_cap=15.0, b_cap=5.0, co_angl=45.0, h_sprk=2.6, h_drft=0.1)
    
    th_cap = deg2rad(th_cap)
    co_angl = deg2rad(co_angl)

    #=
    # Hard-coded defaults 
    ntime   = 200      # number of animation frames
    th_cap  = deg2rad(30.0)      # polar cap tilt angle [rad]  (0 = circular cap in display frame)
    a_cap   = 15.0     # polar cap major semi-axis [m]
    b_cap   = 5.0     # polar cap minor semi-axis [m]  (= a_cap → circular)
    co_angl = deg2rad(45.0)      # co-rotation phase offset [rad]

    # -------------------------------------------------------------------------
    # Derived / physical parameters
    # -------------------------------------------------------------------------
    h_sprk = 2.6   # spark semi-axis (half-size) [m]
    h_drft = 0.1   # raster grid step and drift increment per frame [m]
    =#

    # Spark ellipticity follows the cap ellipticity
    a_sprk = h_sprk
    b_sprk = a_sprk * b_cap / a_cap

    # Centre of the polar cap ellipse in the display (rotated) coordinate frame.
    # These are the projection lengths of the cap semi-axes onto each display axis.
    x_cent = sqrt((a_cap * cos(th_cap))^2 + (b_cap * sin(th_cap))^2)
    y_cent = sqrt((a_cap * sin(th_cap))^2 + (b_cap * cos(th_cap))^2)

    # Number of concentric spark rings that fit inside the cap
    N_trk = floor(Int, b_cap / (2 * b_sprk))   # ≡ floor(a_cap / (2*h_sprk))

    println(stderr, "N_trk = $N_trk   a_cap = $a_cap   b_cap = $b_cap")

    # -------------------------------------------------------------------------
    # Angular spacing between neighbouring sparks on each ring.
    # The number of sparks N_s on a ring is estimated from ring area / spark area
    # (with a 0.75 packing factor), giving spacing 2π / N_s.
    # -------------------------------------------------------------------------
    theta_sp = zeros(Float64, N_trk)
    let a_o = a_cap, a_i = a_o - 2*a_sprk,
        b_o = b_cap, b_i = b_o - 2*b_sprk
        for ring in 1:N_trk
            N_s = floor(Int, 0.75 * (a_o*b_o - a_i*b_i) / (a_sprk*b_sprk))
            theta_sp[ring] = 2π / N_s
            a_o -= 2*a_sprk;  a_i -= 2*a_sprk
            b_o -= 2*b_sprk;  b_i -= 2*b_sprk
        end
    end

    # -------------------------------------------------------------------------
    # Flat (1-D) storage for angular positions of sparks on each ring.
    # Ring `ring` (1-based) occupies slice (ring-1)*trk_max+1 .. (ring-1)*trk_max+trk_max.
    # trk_max is an upper bound on how many sparks can fit on a single half-ring.
    # -------------------------------------------------------------------------
    trk_max = floor(Int,
        0.75 * (a_cap*b_cap - (a_cap - 2*a_sprk)*(b_cap - 2*b_sprk)) /
        (a_sprk * b_sprk) / 2) + 1

    sz        = 2 * trk_max * N_trk + 2
    th_sprk_u = zeros(Float64, sz)   # angular positions of sparks on upper half-rings
    th_sprk_d = zeros(Float64, sz)   # angular positions of sparks on lower half-rings
    N_up      = zeros(Int, N_trk)    # active spark count on upper track, per ring
    N_dn      = zeros(Int, N_trk)    # active spark count on lower track, per ring

    # -------------------------------------------------------------------------
    # Initialise spark angular positions.
    # Upper track starts at π and steps down by theta_sp toward 0.
    # Lower track starts at π and steps up by theta_sp toward 2π.
    # -------------------------------------------------------------------------
    for ring in 1:N_trk
        u_off = (ring - 1) * trk_max
        d_off = (ring - 1) * trk_max

        # Upper track: first spark just inside the leading edge (θ < π)
        th_sprk_u[u_off + 1] = π - theta_sp[ring] / 2 - co_angl
        N_up[ring] = 1
        while th_sprk_u[u_off + N_up[ring]] >= theta_sp[ring] - co_angl
            th_sprk_u[u_off + N_up[ring] + 1] = th_sprk_u[u_off + N_up[ring]] - theta_sp[ring]
            N_up[ring] += 1
        end

        # Lower track: first spark just past the leading edge (θ > π)
        th_sprk_d[d_off + 1] = π + theta_sp[ring] / 2 - co_angl
        N_dn[ring] = 1
        while th_sprk_d[d_off + N_dn[ring]] <= 2π - theta_sp[ring] - co_angl
            th_sprk_d[d_off + N_dn[ring] + 1] = th_sprk_d[d_off + N_dn[ring]] + theta_sp[ring]
            N_dn[ring] += 1
        end
    end

    # -------------------------------------------------------------------------
    # Polar cap boundary ellipse (for the plot outline)
    # -------------------------------------------------------------------------
    ts      = range(0.0, 2π; length = 101)
    x_elips = [coortrans(a_cap*cos(t), b_cap*sin(t), -th_cap)[1] + x_cent for t in ts]
    y_elips = [coortrans(a_cap*cos(t), b_cap*sin(t), -th_cap)[2] + y_cent for t in ts]

    # -------------------------------------------------------------------------
    # GLMakie figure and axes
    # -------------------------------------------------------------------------
    fig       = Figure(size = (600, 600))
    title_obs = Observable("Iteration # 1")     # Observable so the title updates live
    ax = Axis(fig[1, 1];
              xlabel  = "X (m)",
              ylabel  = "Y (m)",
              title   = title_obs,
              aspect  = DataAspect())
    xlims!(ax, 0.0, 2 * x_cent)
    ylims!(ax, 0.0, 2 * y_cent)

    lines!(ax, x_elips, y_elips; color = :black, linewidth = 2)

    spark_pts   = Observable(Point2f[])   # spark centre positions
    spark_sizes = Observable(Vec2f[])     # spark diameters (major, minor) in data units
    scatter!(ax, spark_pts; markersize = spark_sizes, markerspace = :data,
             color = :steelblue, marker = :circle, rotation = -th_cap)

    display(fig)

    # -------------------------------------------------------------------------
    # Angular drift step per frame.
    # Uses the average track radius of the outermost ring (faithful to the C
    # original, where a_out/a_in are not recomputed inside the per-spark loop).
    # -------------------------------------------------------------------------
    mean_outer_radius = 0.5 * (a_cap + (a_cap - 2 * a_sprk))
    del_theta_drift   = h_drft / mean_outer_radius   # arc length → angle [rad/frame]

    # =========================================================================
    # Time-evolution loop
    # =========================================================================
    for step in 1:ntime

        # Compute spark positions and sizes for the current angles
        sx, sy, ss = sparkconfig(th_sprk_u, th_sprk_d, N_up, N_dn, theta_sp,
                                  h_sprk, h_drft, a_cap, b_cap, th_cap,
                                  co_angl, x_cent, y_cent, N_trk, trk_max)

        # Push new spark data to the plot and update the frame counter in the title
        spark_pts[]   = Point2f.(sx, sy)
        spark_sizes[] = Vec2f.(2.0 .* ss, 2.0 .* ss .* b_cap ./ a_cap)   # elliptical: (major, minor) diameters
        title_obs[]   = "Iteration # $step"
        sleep(0.05)   # ~50 ms per frame — matches the delay(50) call in the C original
        println("Number of sparks: ", length(sx))

        # ---------------------------------------------------------------------
        # Advance spark angles by one drift step and rebuild the spark lists.
        # Upper track drifts clockwise (angle decreases).
        # Lower track drifts anti-clockwise (angle increases).
        # When the lead spark passes the boundary it wraps by +/- theta_sp.
        # ---------------------------------------------------------------------
        for ring in 1:N_trk
            u_off = (ring - 1) * trk_max
            d_off = (ring - 1) * trk_max

            # -- Upper track --
            th_sprk_u[u_off + 1] -= del_theta_drift
            # Wrap if the first spark has drifted past the allowed range
            if th_sprk_u[u_off + 1] < π - theta_sp[ring] - co_angl
                th_sprk_u[u_off + 1] += theta_sp[ring]
            end
            # Recompute the full list of upper-track angles from the first spark
            N_up[ring] = 1
            while th_sprk_u[u_off + N_up[ring]] >= theta_sp[ring] - co_angl
                th_sprk_u[u_off + N_up[ring] + 1] = th_sprk_u[u_off + N_up[ring]] - theta_sp[ring]
                N_up[ring] += 1
            end

            # -- Lower track --
            th_sprk_d[d_off + 1] += del_theta_drift
            # Wrap if the first spark has drifted past the allowed range
            if th_sprk_d[d_off + 1] > π + theta_sp[ring] - co_angl
                th_sprk_d[d_off + 1] -= theta_sp[ring]
            end
            # Recompute the full list of lower-track angles from the first spark
            N_dn[ring] = 1
            while th_sprk_d[d_off + N_dn[ring]] <= 2π - theta_sp[ring] - co_angl
                th_sprk_d[d_off + N_dn[ring] + 1] = th_sprk_d[d_off + N_dn[ring]] + theta_sp[ring]
                N_dn[ring] += 1
            end
        end
    end

    println("Animation complete. Close the window to exit.")
    # Keep the process alive until the user closes the Makie window
    while events(fig).window_open[]
        sleep(0.1)
    end
end




end # module LBC
