#=
    drift_check.jl -- diagnostyka bi-driftingu dla zadanej konfiguracji

    Uruchomienie:
        julia --project=. drift_check.jl input/B1839-04.json
        julia --project=. drift_check.jl input/B1839-04.json --scan   # + tabela rf

    Co liczy:
      1. Elipsa czapy polarnej + kontrola, czy nie jest wieksza od dipolowego r_pc.
      2. Przeciecia torow iskier (rfs z JSON-a) z plasma line: dlugosc komponentu
         oraz tempo dryfu wazone intensywnoscia po oknie widocznosci iskry
         (a nie chwilowa pochodna w punkcie najmniejszej odleglosci -- ta potrafi
         miec inny znak niz to, co widac w danych).
      3. LRFS na modelowanych pojedynczych impulsach: gradient fazy przy P3
         w kazdym komponencie. To jest ta sama wielkosc, ktora czyta sie z gornego
         panelu wykresu LRFS, i to ona rozstrzyga o kierunku dryfu.
      4. Porownanie z Tab. 1 z Szary et al. 2020 (ApJ 896, 168).

    Opcja --scan dodaje tabele przeciec dla rf = 0.20..0.98, zeby dobrac "rfs"
    tak, by komponenty wypadly na obserwowanych dlugosciach.
=#

module DriftCheck

    using JSON3
    using LinearAlgebra
    using Statistics
    using FFTW
    using DelimitedFiles

    include("modules/functions.jl")
    include("modules/plot.jl")
    include("modules/field.jl")
    include("modules/nsfield.jl")
    include("modules/lines.jl")
    include("modules/sparks.jl")
    include("modules/signal.jl")
    include("modules/lbc.jl")
    include("modules/pulsar.jl")

    # PSR B1839-04, Tab. 1
    const L_OBS = [-26.0, -16.0, 17.0, 25.0]
    const R_OBS = [-0.69, -0.16, 0.13, 0.13]
    const C_OBS = [18.0, 8.0, 4.0, 2.0]

    ellipse_3d(ef, R, t, rf) = begin
        u = ef.center_local[1] + rf*ef.a*cos(t)*cos(ef.θ) - rf*ef.b*sin(t)*sin(ef.θ)
        v = ef.center_local[2] + rf*ef.a*cos(t)*sin(ef.θ) + rf*ef.b*sin(t)*cos(ef.θ)
        p = ef.centroid + u*ef.x_hat + v*ef.y_hat
        p / norm(p) * R
    end

    """Gesto probkowana plasma line: (punkty 3D, odpowiadajace im dlugosci)."""
    function plasma_line(psr; nref=40)
        los = [[l[1][end], l[2][end], l[3][end]] for l in psr.los_lines]
        lon = psr.longitudes
        pts = Vector{Vector{Float64}}(); lngs = Float64[]
        for i in 1:length(los)-1, k in 0:nref-1
            w = k/nref
            push!(pts, los[i] .+ w .* (los[i+1] .- los[i]))
            push!(lngs, lon[i] + w*(lon[i+1]-lon[i]))
        end
        push!(pts, los[end]); push!(lngs, lon[end])
        return pts, lngs
    end

    """Przeciecia toru rf z plasma line.
    Zwraca [(dlugosc, dphi/dpsi [deg/rad] wazone, szerokosc pasma [deg], min. odleglosc)]."""
    function crossings(psr, ef, rf, pts, lngs; nt=1440)
        R = psr.r
        r_avg = sqrt(ef.a*ef.b); cs = cos(ef.θ); sn = sin(ef.θ)
        base_sigma = psr.spark_radius / 3.72
        phis = zeros(nt); ds = zeros(nt); sig = zeros(nt)
        for k in 1:nt
            s = ellipse_3d(ef, R, 2pi*(k-1)/nt, rf)
            best = Inf; bi = 1
            for (i,p) in enumerate(pts)
                d = 0.0; @inbounds for c in 1:3; d += (p[c]-s[c])^2; end
                if d < best; best = d; bi = i; end
            end
            phis[k] = lngs[bi]; ds[k] = sqrt(best)
            vec = pts[bi] .- s
            u = dot(vec, ef.x_hat); v = dot(vec, ef.y_hat)
            ue = u*cs + v*sn; ve = -u*sn + v*cs
            r2d = hypot(ue, ve)
            sig[k] = r2d < 1e-10 ? base_sigma :
                     base_sigma * (ef.a*ef.b/hypot(ef.b*ue/r2d, ef.a*ve/r2d)) / r_avg
        end
        w = [exp(-ds[k]^2/(2*sig[k]^2)) for k in 1:nt]
        mins = Int[]
        for k in 1:nt
            km = mod1(k-1,nt); kp = mod1(k+1,nt)
            if w[k] >= w[km] && w[k] > w[kp] && w[k] > 0.2; push!(mins, k); end
        end
        sort!(mins, by = k -> -w[k])
        keep = Int[]
        for k in mins
            if all(min(abs(k-j), nt-abs(k-j)) > nt÷12 for j in keep); push!(keep, k); end
            length(keep) == 2 && break
        end
        out = []
        for k in keep
            idx = Int[k]
            for dir in (-1, 1), step in 1:(nt÷6)
                j = mod1(k + dir*step, nt)
                w[j] < 0.05*w[k] && break
                push!(idx, j)
            end
            length(idx) < 5 && continue
            ψ = [2pi*(mod1(j-k+nt÷2, nt) - nt÷2)/nt for j in idx]
            φ = phis[idx]; ww = w[idx]
            Σw = sum(ww); mψ = sum(ww.*ψ)/Σw; mφ = sum(ww.*φ)/Σw
            den = sum(ww.*(ψ .- mψ).^2)
            den < 1e-12 && continue
            slope = sum(ww.*(ψ .- mψ).*(φ .- mφ))/den
            push!(out, (phis[k], slope, maximum(φ)-minimum(φ), ds[k]))
        end
        sort!(out, by = c -> c[1])
        return out
    end

    """LRFS na modelowanych impulsach: zwraca (czestosc P3, moc(bin), faza(bin) w stopniach)."""
    function lrfs(psr)
        d = psr.pulses
        n, nb = size(d)
        dm = d .- mean(d, dims=1)
        F  = rfft(dm, 1)
        pw = dropdims(sum(abs2, F, dims=2), dims=2)
        pw[1] = 0.0                      # bez skladowej stalej
        kpk = argmax(pw)
        f   = (kpk-1)/n                  # w jednostkach 1/P
        return f, vec(abs.(F[kpk, :])), vec(rad2deg.(angle.(F[kpk, :])))
    end

    """Gradient fazy LRFS wokol danej dlugosci, wazony moca.

    Okno dobierane adaptacyjnie: spojny obszar wokol lokalnego maksimum mocy,
    w ktorym moc > frac * moc szczytowa. Staly szeroki przedzial nie dziala --
    pasma dryfu bywaja wezsze niz 4 deg i fit lapie wtedy sasiednie struktury."""
    function phase_gradient(lon, power, phase, centre; frac=0.5, maxhalf=8.0)
        near = findall(i -> abs(lon[i]-centre) <= maxhalf, eachindex(lon))
        isempty(near) && return NaN, 0.0
        k = near[argmax(power[near])]
        pk = power[k]
        lo = k; while lo > 1        && power[lo-1] > frac*pk; lo -= 1; end
        hi = k; while hi < length(lon) && power[hi+1] > frac*pk; hi += 1; end
        idx = collect(lo:hi)
        length(idx) < 4 && return NaN, 0.0
        # rozwiniecie fazy wzdluz okna
        ph = copy(phase[idx])
        for i in 2:length(ph)
            while ph[i] - ph[i-1] >  180; ph[i] -= 360; end
            while ph[i] - ph[i-1] < -180; ph[i] += 360; end
        end
        x = lon[idx]; ww = power[idx]
        Σw = sum(ww); mx = sum(ww.*x)/Σw; mp = sum(ww.*ph)/Σw
        den = sum(ww.*(x .- mx).^2)
        den < 1e-12 && return NaN, 0.0
        return sum(ww.*(x .- mx).*(ph .- mp))/den, maximum(ww)
    end

    function run(jsonfile; do_scan=false)
        psr = Pulsar(jsonfile)
        sc  = psr.sparks_config; si = sc.init
        rfs = collect(si.rfs); N = si.num
        P   = psr.p; P3 = psr.p3[1]

        println("="^78)
        println("  $jsonfile")
        println("  alpha=$(psr.alpha)  beta=$(psr.beta)  D=$(round(Int, psr.r_em/1000)) km  ",
                "P3=$(P3) P  spark_radius=$(psr.spark_radius) m")
        println("  rfs=$(rfs)  num=$N  num_mode=$(get(si, :num_mode, "arc"))")
        println("="^78)

        Lines.init_line_of_sight(psr, num=psr.nsfield.nlos)
        Lines.calculate_line_of_sight(psr)
        Lines.generate_open!(psr, num=psr.nsfield.nopen)
        ef = psr.ellipse_fit

        println("\n[1] CZAPA POLARNA")
        println("    a = $(round(ef.a,digits=2)) m   b = $(round(ef.b,digits=2)) m   ",
                "b/a = $(round(ef.b/ef.a,digits=3))   theta = $(round(rad2deg(ef.θ),digits=1)) deg")
        println("    dipolowe r_pc = $(round(psr.r_pc,digits=1)) m",
                ef.a > psr.r_pc ? "   <<< UWAGA: a > r_pc, dopasowanie elipsy podejrzane" : "")

        pts, lngs = plasma_line(psr)
        dpsi = 2pi/(N*P3)                       # przyrost azymutu iskry na okres

        println("\n[2] PRZECIECIA TOROW Z PLASMA LINE  (tempo wazone po oknie widocznosci)")
        println("    tor      dlugosc     tempo dryfu     szer. pasma   min. odl.")
        model_L = Float64[]; model_R = Float64[]
        for rf in rfs
            for c in crossings(psr, ef, rf, pts, lngs)
                rate = c[2]*dpsi/P
                push!(model_L, c[1]); push!(model_R, rate)
                println("    rf=$(rpad(rf,5))  $(lpad(round(c[1],digits=1),7)) deg  ",
                        "$(lpad(round(rate,digits=3),8)) deg/s  ",
                        "$(lpad(round(c[3],digits=1),8)) deg  ",
                        "$(lpad(round(c[4],digits=2),7)) m")
            end
        end
        println("    OBSERWOWANE:  $(L_OBS) deg,  $(R_OBS) deg/s,  szer. $(C_OBS) deg")
        if length(model_R) == 4
            p = sortperm(model_L)
            println("    R = $(round(sum((R_OBS .- model_R[p]).^2), digits=3))   (prog z publikacji: R < 0.1)")
        end

        if do_scan
            println("\n[2b] SKAN rf  (do doboru \"rfs\": ktory tor tnie plasma line gdzie)")
            for rf in 0.20:0.02:0.98
                cr = crossings(psr, ef, rf, pts, lngs)
                isempty(cr) && continue
                println("    rf=$(rpad(round(rf,digits=2),5))  ",
                        join(["$(lpad(round(c[1],digits=1),7))deg tempo=$(lpad(round(c[2]*dpsi/P,digits=3),7))" for c in cr], " | "))
            end
        end

        println("\n[3] LRFS NA MODELOWANYCH IMPULSACH  (to samo, co gorny panel wykresu)")
        if si.method == "ellipse"
            Sparks.init_sparks1_ellipse!(psr; rfs=rfs, num=N, spacing=get(si,:spacing,"t"),
                phase=get(si,:phase,0.0), num_mode=get(si,:num_mode,"arc"))
        else
            error("drift_check obsluguje na razie tylko init.method = \"ellipse\"")
        end
        Sparks.simulate_sparks_solidbody(psr)
        Signal.generate_signal_new(psr; noise_level=psr.noise_level, v_scale=0.3)
        Signal.generate_pulses(psr)

        # zrzut stosu impulsow -- do bezposredniego sledzenia subpulsow
        dump = "output/drift_check_$(basename(jsonfile)[1:end-5])"
        try
            writedlm(dump*"_pulses.csv", psr.pulses, ',')
            writedlm(dump*"_long.csv", psr.longitudes, ',')
            println("    zrzut: $(dump)_pulses.csv")
        catch e
            println("    (zrzut pominiety: $(sprint(showerror, e)))")
        end

        f, power, phase = lrfs(psr)
        println("    P3 z LRFS = $(round(1/f, digits=2)) P   (zadane: $(P3) P)")
        println("    komponent   moc      grad. fazy [deg/deg]   kierunek")
        pmax = maximum(power)
        senses = String[]
        for L in (isempty(model_L) ? L_OBS : sort(model_L))
            g, w = phase_gradient(psr.longitudes, power, phase, L)
            rel = w/pmax
            # UWAGA: "?" = za malo mocy. Nie uzywac tu "-", bo koliduje ze znakiem ujemnym.
            s = (isnan(g) || rel < 0.05) ? "?" : (g > 0 ? "+" : "-")
            push!(senses, s)
            println("    $(lpad(round(L,digits=1),8)) deg  $(lpad(round(rel,digits=2),5))  ",
                    "$(lpad(isnan(g) ? NaN : round(g,digits=2),12))            $s")
        end
        uniq = unique(filter(!=("-"), senses))
        println("\n    WERDYKT: ", length(uniq) > 1 ?
                ">>> BI-DRIFT (gradient fazy zmienia znak miedzy komponentami) <<<" :
                "brak bi-driftu (jeden kierunek w calym profilu)")
        println("    (znaki wzgledne; sama konwencja znaku FFT nie ma znaczenia,")
        println("     liczy sie roznica miedzy komponentami wiodacymi a tylnymi)")
        return psr
    end
end

let
    file = length(ARGS) >= 1 ? ARGS[1] : "input/B1839-04.json"
    DriftCheck.run(file; do_scan = "--scan" in ARGS)
end
