# Plan: dodanie polarization position angle (PA) do symulacji

Dokument opisuje plan implementacji generowania PA w symulacji. Dwa warianty:
**A** — geometryczny RVM (jedna krzywa PA na profil), **B** — pełna polaryzacja
per-pulse z sumowaniem Stokesów (Q, U, depolaryzacja, opcjonalnie OPM).

## Status wyjściowy

- `Signal.generate_signal_radii` produkuje wyłącznie intensywność `psr.signal`.
- W `Lines.calculate_line_of_sight` w pierwszej iteracji pętli liczone jest **B**
  w punkcie emisji `psr.line_of_sight[i]`, ale wartość jest tracona po jednym kroku.
- Oś rotacji w układzie osi magnetycznej jest znana: Ω̂ = (sin α, 0, cos α).
- LOS w tym układzie jest radialny: n̂_j = `psr.line_of_sight[j] / r_em`.

## Wariant A — geometryczny RVM

Jedna krzywa PA(φ) wyliczona z **B** w punkcie LOS. Nie zależy od sparków.

### Zmiany

1. **Struktura `Pulsar`** ([main.jl:14](../main.jl#L14)):
   - `b_em::Vector{Vector{Float64}}` — **B** kartezjańsko w każdym punkcie emisji LOS.
   - `pa::Vector{Float64}` — PA na bin longitude (długość = `psr.nsfield.nlos`).
   - Inicjalizacja w obu konstruktorach.

2. **`Lines.calculate_line_of_sight`** ([lines.jl:218](../modules/lines.jl#L218))
   i wersja dipolowa ([lines.jl:163](../modules/lines.jl#L163)):
   - Po pierwszym wyliczeniu `b` w punkcie emisji: `push!(psr.b_em, b)`.

3. **Nowa funkcja `Signal.calculate_pa!(psr)`** w [signal.jl](../modules/signal.jl):

   ```julia
   function calculate_pa!(psr)
       α = deg2rad(psr.alpha)
       Ω = [sin(α), 0.0, cos(α)]
       n = length(psr.line_of_sight)
       psr.pa = zeros(n)
       for i in 1:n
           nhat = psr.line_of_sight[i] / norm(psr.line_of_sight[i])
           B    = psr.b_em[i]
           Bp   = B - dot(B, nhat) * nhat
           Ωp   = Ω - dot(Ω, nhat) * nhat
           e0   = Ωp / norm(Ωp)
           e1   = cross(nhat, e0)
           psr.pa[i] = atan(dot(Bp, e1), dot(Bp, e0))
           # Konwencja: PA wektora E = PA(B) + π/2; wybrać raz i udokumentować.
       end
   end
   ```

4. **Wywołanie** w `generate_signal()` ([main.jl:226](../main.jl#L226))
   i `generate_signal_dipole()` po `Signal.generate_signal_radii(...)`:
   ```julia
   Signal.calculate_pa!(psr)
   ```

5. **Plot** w `plot.jl`: `psr.longitudes` vs `rad2deg.(psr.pa)`.

### Skala roboty: ~30-50 linii.

---

## Wariant B — pełna polaryzacja per-pulse (Stokes Q, U)

Każdy spark emituje z własnego punktu na linii pola; **B** w tym punkcie wyznacza
PA fotonu. Wkłady sumowane niekoherentnie jako Q, U → naturalna depolaryzacja
i odchyłki od czystego RVM, gdy LOS próbkuje sparki o różnych χ.

### Model

Dla sparku *k* w pozycji **s**_k:

1. e_k — punkt przebicia linii pola przez sferę r=r_em.
2. **B**_k = B(e_k).
3. n̂_k = e_k / |e_k|.
4. χ_k z **B**_k, n̂_k, Ω̂ (jak w wariancie A, ale per-spark).
5. Wkład do biny LOS *j* w pulsie *i*:
   ```
   I += w_kj
   Q += w_kj · cos(2 χ_k)
   U += w_kj · sin(2 χ_k)
   ```
   gdzie `w_kj = exp(-d²/2σ²)` (waga Gaussa, jak obecnie).

6. Z tego:
   ```
   PA(i,j) = 0.5 * atan2(U, Q)
   L(i,j)  = sqrt(Q² + U²)
   p_lin   = L / I
   ```

### Mapa „spark → punkt emisji + B"

Dwa warianty implementacji:

- **Dipol (zamknięta forma)**: linia dipolowa r = r_max·sin²θ; znając (θ_s, ψ_s)
  na powierzchni dostajesz θ_em z `sin²θ_em = (r_em/r_*)·sin²θ_s`, ψ_em = ψ_s.
- **Nie-dipol (trasowanie)**: pętla od **s**_k w górę (`st = +b/|b|·step`) do
  r=r_em — odwrotny krok do tego co w `calculate_line_of_sight`. Drogie per-spark
  per-pulse → preferowany **lookup**: raz na siatce powierzchniowej (np. 50×50)
  policzyć `(θ_s, ψ_s) → (θ_em, ψ_em, B)`, w pętli interpolować biliniowo.

### Zmiany

1. **Struktura `Pulsar`**:
   - `spark_emission_points::Vector{Vector{Vector{Float64}}}` indeksowane `[step][k]`.
   - `spark_b_em::Vector{Vector{Vector{Float64}}}` — analogicznie.
   - `stokes_q`, `stokes_u` — macierze równoległe do `psr.signal`.
   - `pa_pulses`, `lin_pol` — wyniki polaryzacji per (puls, longitude).

2. **`Sparks.compute_emission_map!(psr)`** — populuje `spark_emission_points`
   i `spark_b_em` (lookup lub closed-form dipolowo).

3. **Rozszerzenie `Signal.generate_signal_radii`** ([signal.jl:33](../modules/signal.jl#L33)):
   ```julia
   χ = pa_from_B(psr.spark_b_em[j][k], psr.spark_emission_points[j][k], α)
   contrib = exp(-dist^2 / (2σ^2))
   psr.signal[j,i]   += contrib
   psr.stokes_q[j,i] += contrib * cos(2χ)
   psr.stokes_u[j,i] += contrib * sin(2χ)
   ```

4. **`Signal.compute_polarization!(psr)`**: z Q, U liczy `pa_pulses`, `lin_pol`,
   integrowany profil.

5. **Ploty**: PA per-puls (chmura punktów PA(φ) ze wszystkich pulsów),
   profil L/I uśredniony po pulsach.

### Niuanse

1. **Geometria wagi**: obecny Gauss mierzy odległość *na powierzchni*. Spójniej
   fizycznie byłoby liczyć w *punkcie emisji* (`|psr.line_of_sight[j] - e_k|`).
   Decyzja merytoryczna do podjęcia przed implementacją.

2. **Skoki ortogonalne (OPM)**: czysty model nie wygeneruje skoków o 90°. Aby je
   uzyskać, dodać per-sparkowi flagę `mode ∈ {0, 1}` z ustalonym prawdopodobieństwem
   i zamienić `χ_k → χ_k + mode·π/2`. ~3 linie.

3. **Stokes V** (polaryzacja kołowa): wymaga osobnego mechanizmu (gradient
   częstotliwościowy, propagacja). Poza zakresem.

4. **Konwencja PA**: zdecydować raz, czy `pa` odnosi się do wektora **B** czy
   wektora **E** (różnica π/2). Standard pulsarowy → wektor **E**.

### Skala roboty

- Wariant dipolowy + OPM: ~popołudnie.
- Wariant nie-dipolowy z lookupem: ~dzień.
- Najwięcej myślenia idzie w decyzję z punktu 1 (powierzchnia vs r_em w wadze).

---

## Kolejność implementacji

1. Wariant A (RVM) — szybko, bezpiecznie, daje sanity check do walidacji wariantu B.
2. Wariant B dipolowy — sprawdzenie modelu sumowania Stokesów.
3. Wariant B nie-dipolowy z lookupem — pełna funkcjonalność.
4. OPM jako opcjonalna nakładka.
