# =============================================================================
# Statistiques et visualisations du pendule double
# =============================================================================

include("pendule.jl")

# =============================================================================
# Fonctions de métriques pour comparer simulation et mesures
# =============================================================================

function calculate_R2(y_measured, y_simulated)
    """Calcule le coefficient de détermination R²"""
    ss_res = sum((y_measured .- y_simulated).^2)
    ss_tot = sum((y_measured .- mean(y_measured)).^2)
    return 1 - ss_res / ss_tot
end

function calculate_RMSE(y_measured, y_simulated)
    """Calcule l'erreur quadratique moyenne (RMSE)"""
    return sqrt(mean((y_measured .- y_simulated).^2))
end

function calculate_MAE(y_measured, y_simulated)
    """Calcule l'erreur absolue moyenne (MAE)"""
    return mean(abs.(y_measured .- y_simulated))
end

function calculate_max_error(y_measured, y_simulated)
    """Calcule l'erreur maximale"""
    return maximum(abs.(y_measured .- y_simulated))
end

# =============================================================================
# Exposant de Lyapunov (chaos théorique)
# =============================================================================

function lyapunov_exponent(m1, m2, L1, L2, g, θ1_0, θ2_0, T_max; δ0=1e-8, dt=0.01)
    """
    Calcule l'exposant de Lyapunov maximal du pendule double.
    
    L'exposant de Lyapunov mesure la sensibilité aux conditions initiales.
    Un exposant positif indique un comportement chaotique.
    
    Méthode: On simule deux trajectoires avec une petite perturbation δ0,
    puis on mesure comment cette séparation évolue au cours du temps.
    
    λ = lim(t→∞) (1/t) * ln(|δ(t)| / |δ0|)
    """
    
    # Conditions initiales pour la trajectoire de référence
    u0_ref = [θ1_0, 0.0, θ2_0, 0.0]
    
    # Conditions initiales perturbées (petite perturbation sur θ1)
    u0_pert = [θ1_0 + δ0, 0.0, θ2_0, 0.0]
    
    p = [m1, m2, L1, L2, g]
    tspan = (0.0, T_max)
    t_eval = collect(0:dt:T_max)
    
    # Simuler les deux trajectoires (modèle idéal sans frottements)
    prob_ref = ODEProblem(equations_double_pendulum_ideal!, u0_ref, tspan, p)
    prob_pert = ODEProblem(equations_double_pendulum_ideal!, u0_pert, tspan, p)
    
    sol_ref = solve(prob_ref, Tsit5(), saveat=t_eval)
    sol_pert = solve(prob_pert, Tsit5(), saveat=t_eval)
    
    # Calculer la séparation au cours du temps
    n_points = length(t_eval)
    separations = zeros(n_points)
    
    for i in 1:n_points
        # Distance dans l'espace des phases (θ1, ω1, θ2, ω2)
        δθ1 = sol_pert[1, i] - sol_ref[1, i]
        δω1 = sol_pert[2, i] - sol_ref[2, i]
        δθ2 = sol_pert[3, i] - sol_ref[3, i]
        δω2 = sol_pert[4, i] - sol_ref[4, i]
        
        separations[i] = sqrt(δθ1^2 + δω1^2 + δθ2^2 + δω2^2)
    end
    
    # Éviter les valeurs nulles ou négatives pour le log
    separations = max.(separations, 1e-15)
    
    # Calculer l'exposant de Lyapunov par régression linéaire de ln(δ) vs t
    log_sep = log.(separations)
    
    # Régression linéaire: ln(δ) = λ*t + ln(δ0)
    # On utilise seulement les points où la séparation n'a pas saturé
    valid_idx = findall(separations .< 10.0)  # Éviter la saturation
    
    if length(valid_idx) > 10
        t_valid = t_eval[valid_idx]
        log_sep_valid = log_sep[valid_idx]
        
        # Régression linéaire simple
        n = length(t_valid)
        sum_t = sum(t_valid)
        sum_log = sum(log_sep_valid)
        sum_t2 = sum(t_valid.^2)
        sum_t_log = sum(t_valid .* log_sep_valid)
        
        λ = (n * sum_t_log - sum_t * sum_log) / (n * sum_t2 - sum_t^2)
    else
        λ = NaN
    end
    
    return λ, t_eval, separations
end

function lyapunov_spectrum(m1, m2, L1, L2, g, θ1_0, θ2_0, T_max; n_perturbations=10, δ0=1e-8, dt=0.01)
    """
    Estime le spectre de Lyapunov en moyennant sur plusieurs directions de perturbation.
    """
    λ_values = Float64[]
    
    # Perturbations dans différentes directions de l'espace des phases
    perturbation_directions = [
        [1.0, 0.0, 0.0, 0.0],  # θ1
        [0.0, 1.0, 0.0, 0.0],  # ω1
        [0.0, 0.0, 1.0, 0.0],  # θ2
        [0.0, 0.0, 0.0, 1.0],  # ω2
    ]
    
    for dir in perturbation_directions
        u0_ref = [θ1_0, 0.0, θ2_0, 0.0]
        u0_pert = u0_ref .+ δ0 .* dir
        
        p = [m1, m2, L1, L2, g]
        tspan = (0.0, T_max)
        t_eval = collect(0:dt:T_max)
        
        prob_ref = ODEProblem(equations_double_pendulum_ideal!, u0_ref, tspan, p)
        prob_pert = ODEProblem(equations_double_pendulum_ideal!, u0_pert, tspan, p)
        
        sol_ref = solve(prob_ref, Tsit5(), saveat=t_eval)
        sol_pert = solve(prob_pert, Tsit5(), saveat=t_eval)
        
        # Séparation finale
        δ_final = sqrt(sum((sol_pert[:, end] .- sol_ref[:, end]).^2))
        
        if δ_final > 1e-14
            λ = log(δ_final / δ0) / T_max
            push!(λ_values, λ)
        end
    end
    
    return λ_values
end

# =============================================================================
# Divergence mesures vs simulation
# =============================================================================

function analyze_divergence(t_data, θ1_data, θ2_data, θ1_sim, θ2_sim; threshold_deg=10.0)
    """
    Analyse la divergence entre les mesures expérimentales et la simulation.
    
    Retourne:
    - t_divergence: temps où l'erreur dépasse le seuil
    - errors_θ1, errors_θ2: erreurs au cours du temps
    - error_total: erreur combinée (norme)
    - divergence_rate: taux de croissance de l'erreur (pseudo-Lyapunov)
    """
    
    n = length(t_data)
    errors_θ1 = abs.(θ1_data .- θ1_sim)
    errors_θ2 = abs.(θ2_data .- θ2_sim)
    error_total = sqrt.(errors_θ1.^2 .+ errors_θ2.^2)
    
    # Seuil en radians
    threshold_rad = deg2rad(threshold_deg)
    
    # Trouver le moment de divergence (quand l'erreur dépasse le seuil)
    t_divergence = NaN
    idx_divergence = findfirst(error_total .> threshold_rad)
    if idx_divergence !== nothing
        t_divergence = t_data[idx_divergence]
    end
    
    # Calculer le taux de divergence (fit exponentiel sur la partie croissante)
    # On cherche: error(t) ≈ error(0) * exp(λ_div * t)
    # Donc: ln(error) ≈ ln(error(0)) + λ_div * t
    
    # Utiliser seulement les points où l'erreur est significative et croissante
    valid_idx = findall((error_total .> 1e-6) .& (error_total .< 1.0))
    
    divergence_rate = NaN
    if length(valid_idx) > 10
        t_valid = t_data[valid_idx]
        log_err_valid = log.(error_total[valid_idx])
        
        # Régression linéaire
        n_pts = length(t_valid)
        sum_t = sum(t_valid)
        sum_log = sum(log_err_valid)
        sum_t2 = sum(t_valid.^2)
        sum_t_log = sum(t_valid .* log_err_valid)
        
        divergence_rate = (n_pts * sum_t_log - sum_t * sum_log) / (n_pts * sum_t2 - sum_t^2)
    end
    
    # Temps de prédictibilité (temps pour que l'erreur soit multipliée par e)
    predictability_time = divergence_rate > 0 ? 1.0 / divergence_rate : Inf
    
    return (
        t_divergence = t_divergence,
        errors_θ1 = errors_θ1,
        errors_θ2 = errors_θ2,
        error_total = error_total,
        divergence_rate = divergence_rate,
        predictability_time = predictability_time,
        threshold_deg = threshold_deg
    )
end

function find_best_fit_window(t_data, θ1_data, θ2_data, θ1_sim, θ2_sim; max_error_deg=5.0)
    """
    Trouve la fenêtre temporelle où la simulation reste proche des mesures.
    Utile pour déterminer jusqu'à quand le modèle est fiable.
    """
    max_error_rad = deg2rad(max_error_deg)
    errors = sqrt.((θ1_data .- θ1_sim).^2 .+ (θ2_data .- θ2_sim).^2)
    
    # Trouver le dernier index où l'erreur est acceptable
    last_good_idx = findlast(errors .<= max_error_rad)
    
    if last_good_idx === nothing
        return 0.0, 0
    else
        return t_data[last_good_idx], last_good_idx
    end
end

# =============================================================================
# Énergie du système
# =============================================================================

function calculate_energy(θ1, ω1, θ2, ω2, m1, m2, L1, L2, g)
    """Calcule l'énergie totale du pendule double"""
    
    # Énergie cinétique
    v1_sq = (L1 * ω1)^2
    v2_sq = (L1 * ω1)^2 + (L2 * ω2)^2 + 2 * L1 * L2 * ω1 * ω2 * cos(θ1 - θ2)
    
    T = 0.5 * m1 * v1_sq + 0.5 * m2 * v2_sq
    
    # Énergie potentielle (référence: pivot)
    y1 = -L1 * cos(θ1)
    y2 = y1 - L2 * cos(θ2)
    
    V = m1 * g * y1 + m2 * g * y2
    
    return T + V, T, V
end

function energy_over_time(t_data, sol, m1, m2, L1, L2, g)
    """Calcule l'évolution de l'énergie au cours du temps"""
    E_total = Float64[]
    E_kinetic = Float64[]
    E_potential = Float64[]
    
    for i in 1:length(t_data)
        E, T, V = calculate_energy(sol[1,i], sol[2,i], sol[3,i], sol[4,i], m1, m2, L1, L2, g)
        push!(E_total, E)
        push!(E_kinetic, T)
        push!(E_potential, V)
    end
    
    return E_total, E_kinetic, E_potential
end

# =============================================================================
# Animation du pendule double
# =============================================================================

function animate_pendulum(t_data, θ1_data, θ2_data, θ1_sim, θ2_sim, L1, L2; fps=30, skip=2)
    """Crée une animation GIF du pendule double (mesuré vs simulé)"""
    
    function pendulum_positions(θ1, θ2, L1, L2)
        x1 = -L1 * sin(θ1)
        y1 = -L1 * cos(θ1)
        x2 = x1 - L2 * sin(θ2)
        y2 = y1 - L2 * cos(θ2)
        return x1, y1, x2, y2
    end
    
    L_total = (L1 + L2) * 1.3
    
    # Pré-calculer toutes les positions de la masse 2 pour le tracé
    trace_x2_m = Float64[]
    trace_y2_m = Float64[]
    trace_x2_s = Float64[]
    trace_y2_s = Float64[]
    
    for j in 1:length(t_data)
        _, _, x2_m, y2_m = pendulum_positions(θ1_data[j], θ2_data[j], L1, L2)
        _, _, x2_s, y2_s = pendulum_positions(θ1_sim[j], θ2_sim[j], L1, L2)
        push!(trace_x2_m, x2_m)
        push!(trace_y2_m, y2_m)
        push!(trace_x2_s, x2_s)
        push!(trace_y2_s, y2_s)
    end
    
    println("Création de l'animation...")
    
    anim = @animate for i in 1:skip:length(t_data)
        x1_m, y1_m, x2_m, y2_m = pendulum_positions(θ1_data[i], θ2_data[i], L1, L2)
        x1_s, y1_s, x2_s, y2_s = pendulum_positions(θ1_sim[i], θ2_sim[i], L1, L2)
        
        plt = plot(
            xlim=(-L_total, L_total),
            ylim=(-L_total, L_total),
            aspect_ratio=:equal,
            legend=:topright,
            title="Pendule Double - t = $(round(t_data[i], digits=2)) s",
            xlabel="x (m)",
            ylabel="y (m)",
            size=(600, 600)
        )
        
        # Tracé de la masse 2
        plot!(plt, trace_x2_m[1:i], trace_y2_m[1:i], lw=1, color=:blue, alpha=0.5, label="Tracé mesuré")
        plot!(plt, trace_x2_s[1:i], trace_y2_s[1:i], lw=1, color=:red, alpha=0.5, label="Tracé simulé")
        
        # Pendules
        plot!(plt, [0, x1_m, x2_m], [0, y1_m, y2_m], lw=3, color=:blue, label="Mesuré", marker=:circle, markersize=8)
        plot!(plt, [0, x1_s, x2_s], [0, y1_s, y2_s], lw=3, color=:red, ls=:dash, label="Simulé", marker=:circle, markersize=8)
        
        scatter!(plt, [0], [0], color=:black, markersize=10, label="Pivot")
    end
    
    gif(anim, "assets/pendule_animation.gif", fps=fps)
    println("Animation sauvegardée: assets/pendule_animation.gif")
    
    # Image du tracé final
    plt_final = plot(
        xlim=(-L_total, L_total),
        ylim=(-L_total, L_total),
        aspect_ratio=:equal,
        legend=:topright,
        title="Tracé complet de la masse 2",
        xlabel="x (m)",
        ylabel="y (m)",
        size=(600, 600)
    )
    
    plot!(plt_final, trace_x2_m, trace_y2_m, lw=2, color=:blue, label="Mesuré")
    plot!(plt_final, trace_x2_s, trace_y2_s, lw=2, color=:red, ls=:dash, label="Simulé")
    scatter!(plt_final, [0], [0], color=:black, markersize=10, label="Pivot")
    
    savefig(plt_final, "assets/pendule_trace.png")
    println("Tracé final sauvegardé: assets/pendule_trace.png")
end

# =============================================================================
# Calcul et affichage des statistiques
# =============================================================================

println("\n" * "="^60)
println("STATISTIQUES DU PENDULE DOUBLE")
println("="^60)

# Métriques de comparaison
R2_θ1 = calculate_R2(θ1_data, θ1_sim_corrected)
R2_θ2 = calculate_R2(θ2_data, θ2_sim_corrected)
RMSE_θ1 = calculate_RMSE(θ1_data, θ1_sim_corrected)
RMSE_θ2 = calculate_RMSE(θ2_data, θ2_sim_corrected)
MAE_θ1 = calculate_MAE(θ1_data, θ1_sim_corrected)
MAE_θ2 = calculate_MAE(θ2_data, θ2_sim_corrected)
max_err_θ1 = calculate_max_error(θ1_data, θ1_sim_corrected)
max_err_θ2 = calculate_max_error(θ2_data, θ2_sim_corrected)

println("\n--- Métriques de comparaison simulation/mesures ---")
println("Coefficient de détermination R²:")
println("  θ1: $R2_θ1")
println("  θ2: $R2_θ2")
println("\nErreur quadratique moyenne (RMSE):")
println("  θ1: $(rad2deg(RMSE_θ1))°  ($(RMSE_θ1) rad)")
println("  θ2: $(rad2deg(RMSE_θ2))°  ($(RMSE_θ2) rad)")
println("\nErreur absolue moyenne (MAE):")
println("  θ1: $(rad2deg(MAE_θ1))°")
println("  θ2: $(rad2deg(MAE_θ2))°")
println("\nErreur maximale:")
println("  θ1: $(rad2deg(max_err_θ1))°")
println("  θ2: $(rad2deg(max_err_θ2))°")

# Exposant de Lyapunov
println("\n--- Exposant de Lyapunov ---")
λ, t_lyap, separations = lyapunov_exponent(m_opt[1], m_opt[2], L1, L2, g, -θ1_0, θ2_0, 2.0)
println("Exposant de Lyapunov maximal: λ = $λ s⁻¹")

if λ > 0
    println("  → Système CHAOTIQUE (λ > 0)")
    println("  → Temps de prédictibilité (Lyapunov time): $(1/λ) s")
else
    println("  → Système régulier (λ ≤ 0)")
end

# Spectre de Lyapunov
λ_spectrum = lyapunov_spectrum(m_opt[1], m_opt[2], L1, L2, g, -θ1_0, θ2_0, 2.0)
println("\nSpectre de Lyapunov (4 directions):")
for (i, λi) in enumerate(λ_spectrum)
    dir_names = ["θ1", "ω1", "θ2", "ω2"]
    println("  λ_$(dir_names[i]) = $λi s⁻¹")
end

# Énergie
println("\n--- Conservation de l'énergie (simulation) ---")
E_total, E_kinetic, E_potential = energy_over_time(t_data, sol_opt, m_opt[1], m_opt[2], L1, L2, g)
E_variation = (maximum(E_total) - minimum(E_total)) / abs(mean(E_total)) * 100
println("Énergie initiale: $(E_total[1]) J")
println("Énergie finale: $(E_total[end]) J")
println("Variation relative: $(E_variation)%")

# Divergence mesures vs simulation
println("\n--- Divergence mesures vs simulation ---")
div_analysis = analyze_divergence(t_data, θ1_data, θ2_data, θ1_sim_corrected, θ2_sim_corrected, threshold_deg=10.0)

println("Seuil de divergence: $(div_analysis.threshold_deg)°")
if !isnan(div_analysis.t_divergence)
    println("⚠ Moment de divergence: t = $(round(div_analysis.t_divergence, digits=3)) s")
    println("  (l'erreur dépasse $(div_analysis.threshold_deg)° à partir de ce moment)")
else
    println("✓ Pas de divergence > $(div_analysis.threshold_deg)° sur la durée mesurée")
end

println("\nTaux de divergence (pseudo-Lyapunov): $(round(div_analysis.divergence_rate, digits=3)) s⁻¹")
if div_analysis.divergence_rate > 0
    println("  → Temps de prédictibilité du modèle: $(round(div_analysis.predictability_time, digits=3)) s")
    println("  → L'erreur double tous les $(round(log(2)/div_analysis.divergence_rate, digits=3)) s")
end

# Fenêtre de validité du modèle
t_valid, idx_valid = find_best_fit_window(t_data, θ1_data, θ2_data, θ1_sim_corrected, θ2_sim_corrected, max_error_deg=5.0)
println("\nFenêtre de validité (erreur < 5°): 0 à $(round(t_valid, digits=3)) s")
if idx_valid > 0
    println("  → Le modèle est fiable pour $(round(100*idx_valid/length(t_data), digits=1))% de la durée")
end

# =============================================================================
# Graphiques
# =============================================================================

# Graphique de comparaison des angles
p_angles = plot(t_data, θ1_data, label="θ1 mesuré", lw=2, title="Comparaison simulation vs données")
plot!(p_angles, t_data, θ1_sim_corrected, label="θ1 simulé", ls=:dash, lw=2)
plot!(p_angles, t_data, θ2_data, label="θ2 mesuré", lw=2)
plot!(p_angles, t_data, θ2_sim_corrected, label="θ2 simulé", ls=:dash, lw=2)
xlabel!(p_angles, "Temps (s)")
ylabel!(p_angles, "Angle (rad)")
savefig(p_angles, "assets/pendule_comparison.png")
println("\nGraphique sauvegardé: assets/pendule_comparison.png")

# Graphique de l'exposant de Lyapunov (divergence)
p_lyap = plot(t_lyap, log.(separations), label="ln(δ)", lw=2, title="Divergence des trajectoires")
plot!(p_lyap, t_lyap, λ .* t_lyap .+ log(1e-8), label="Fit: λt + ln(δ₀)", ls=:dash, lw=2)
xlabel!(p_lyap, "Temps (s)")
ylabel!(p_lyap, "ln(séparation)")
savefig(p_lyap, "assets/lyapunov_divergence.png")
println("Graphique sauvegardé: assets/lyapunov_divergence.png")

# Graphique de l'énergie
p_energy = plot(t_data, E_total, label="E totale", lw=2, title="Énergie du système (simulation)")
plot!(p_energy, t_data, E_kinetic, label="E cinétique", lw=2)
plot!(p_energy, t_data, E_potential, label="E potentielle", lw=2)
xlabel!(p_energy, "Temps (s)")
ylabel!(p_energy, "Énergie (J)")
savefig(p_energy, "assets/energie.png")
println("Graphique sauvegardé: assets/energie.png")

# Graphique de la divergence mesures vs simulation
p_div = plot(t_data, rad2deg.(div_analysis.errors_θ1), label="Erreur θ1", lw=2, 
             title="Divergence mesures vs simulation", color=:blue)
plot!(p_div, t_data, rad2deg.(div_analysis.errors_θ2), label="Erreur θ2", lw=2, color=:orange)
plot!(p_div, t_data, rad2deg.(div_analysis.error_total), label="Erreur totale", lw=2, color=:red)
hline!(p_div, [div_analysis.threshold_deg], label="Seuil $(div_analysis.threshold_deg)°", ls=:dash, color=:black)
if !isnan(div_analysis.t_divergence)
    vline!(p_div, [div_analysis.t_divergence], label="Divergence", ls=:dot, color=:gray, lw=2)
end
xlabel!(p_div, "Temps (s)")
ylabel!(p_div, "Erreur (°)")
savefig(p_div, "assets/divergence_mesures_simulation.png")
println("Graphique sauvegardé: assets/divergence_mesures_simulation.png")

# Graphique log de la divergence (pour voir le taux exponentiel)
# Éviter les valeurs nulles pour l'échelle log
error_for_log = max.(div_analysis.error_total, 1e-10)
p_div_log = plot(t_data, error_for_log, label="Erreur totale", lw=2, 
                 title="Divergence mesures vs simulation (échelle log)", yscale=:log10)
xlabel!(p_div_log, "Temps (s)")
ylabel!(p_div_log, "Erreur (rad)")
savefig(p_div_log, "assets/divergence_log.png")
println("Graphique sauvegardé: assets/divergence_log.png")

# Animation
animate_pendulum(t_data, θ1_data, θ2_data, θ1_sim_corrected, θ2_sim_corrected, L1, L2, fps=10, skip=2)

println("\n" * "="^60)
println("Statistiques terminées!")
println("="^60)
