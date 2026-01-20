include("pendule.jl")

# =============================================================================
# Fonctions de métriques
# =============================================================================

function calculate_R2(y_measured, y_simulated)
    ss_res = sum((y_measured .- y_simulated).^2)
    ss_tot = sum((y_measured .- mean(y_measured)).^2)
    return 1 - ss_res / ss_tot
end

function calculate_RMSE(y_measured, y_simulated)
    return sqrt(mean((y_measured .- y_simulated).^2))
end

function calculate_MAE(y_measured, y_simulated)
    return mean(abs.(y_measured .- y_simulated))
end

function calculate_max_error(y_measured, y_simulated)
    return maximum(abs.(y_measured .- y_simulated))
end

# =============================================================================
# Exposant de Lyapunov
# =============================================================================

function lyapunov_exponent(m1, m2, L1, L2, g, θ1_0, θ2_0, T_max; δ0=1e-8, dt=0.01)
    u0_ref = [θ1_0, 0.0, θ2_0, 0.0]
    u0_pert = [θ1_0 + δ0, 0.0, θ2_0, 0.0]
    
    p = [m1, m2, L1, L2, g]
    tspan = (0.0, T_max)
    t_eval = collect(0:dt:T_max)
    
    prob_ref = ODEProblem(equations_double_pendulum, u0_ref, tspan, p)
    prob_pert = ODEProblem(equations_double_pendulum, u0_pert, tspan, p)
    
    sol_ref = solve(prob_ref, Tsit5(), saveat=t_eval)
    sol_pert = solve(prob_pert, Tsit5(), saveat=t_eval)
    
    n_points = length(t_eval)
    separations = zeros(n_points)
    
    for i in 1:n_points
        δθ1 = sol_pert[1, i] - sol_ref[1, i]
        δω1 = sol_pert[2, i] - sol_ref[2, i]
        δθ2 = sol_pert[3, i] - sol_ref[3, i]
        δω2 = sol_pert[4, i] - sol_ref[4, i]
        
        separations[i] = sqrt(δθ1^2 + δω1^2 + δθ2^2 + δω2^2)
    end
    
    separations = max.(separations, 1e-15)
    log_sep = log.(separations)
    
    valid_idx = findall(separations .< 10.0)
    
    if length(valid_idx) > 10
        t_valid = t_eval[valid_idx]
        log_sep_valid = log_sep[valid_idx]
        
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
    λ_values = Float64[]
    
    perturbation_directions = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    
    for dir in perturbation_directions
        u0_ref = [θ1_0, 0.0, θ2_0, 0.0]
        u0_pert = u0_ref .+ δ0 .* dir
        
        p = [m1, m2, L1, L2, g]
        tspan = (0.0, T_max)
        t_eval = collect(0:dt:T_max)
        
        prob_ref = ODEProblem(equations_double_pendulum, u0_ref, tspan, p)
        prob_pert = ODEProblem(equations_double_pendulum, u0_pert, tspan, p)
        
        sol_ref = solve(prob_ref, Tsit5(), saveat=t_eval)
        sol_pert = solve(prob_pert, Tsit5(), saveat=t_eval)
        
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
    n = length(t_data)
    errors_θ1 = abs.(θ1_data .- θ1_sim)
    errors_θ2 = abs.(θ2_data .- θ2_sim)
    error_total = sqrt.(errors_θ1.^2 .+ errors_θ2.^2)
    
    threshold_rad = deg2rad(threshold_deg)
    
    t_divergence = NaN
    idx_divergence = findfirst(error_total .> threshold_rad)
    if idx_divergence !== nothing
        t_divergence = t_data[idx_divergence]
    end
    
    valid_idx = findall((error_total .> 1e-6) .& (error_total .< 1.0))
    
    divergence_rate = NaN
    if length(valid_idx) > 10
        t_valid = t_data[valid_idx]
        log_err_valid = log.(error_total[valid_idx])
        
        n_pts = length(t_valid)
        sum_t = sum(t_valid)
        sum_log = sum(log_err_valid)
        sum_t2 = sum(t_valid.^2)
        sum_t_log = sum(t_valid .* log_err_valid)
        
        divergence_rate = (n_pts * sum_t_log - sum_t * sum_log) / (n_pts * sum_t2 - sum_t^2)
    end
    
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
    max_error_rad = deg2rad(max_error_deg)
    errors = sqrt.((θ1_data .- θ1_sim).^2 .+ (θ2_data .- θ2_sim).^2)
    
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
    v1_sq = (L1 * ω1)^2
    v2_sq = (L1 * ω1)^2 + (L2 * ω2)^2 + 2 * L1 * L2 * ω1 * ω2 * cos(θ1 - θ2)
    
    T = 0.5 * m1 * v1_sq + 0.5 * m2 * v2_sq
    
    y1 = -L1 * cos(θ1)
    y2 = y1 - L2 * cos(θ2)
    
    V = m1 * g * y1 + m2 * g * y2
    
    return T + V, T, V
end

function energy_over_time(t_data, sol, m1, m2, L1, L2, g)
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
# Animations
# =============================================================================

function animate_pendulum(t_data, θ1_data, θ2_data, θ1_sim, θ2_sim, L1, L2; fps=30, skip=2)
    function pendulum_positions(θ1, θ2, L1, L2)
        x1 = L1 * sin(θ1)
        y1 = -L1 * cos(θ1)
        x2 = x1 + L2 * sin(θ2)
        y2 = y1 - L2 * cos(θ2)
        return x1, y1, x2, y2
    end
    
    L_total = (L1 + L2) * 1.3
    
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
            title="t = $(round(t_data[i], digits=2)) s",
            xlabel="x (m)",
            ylabel="y (m)",
            size=(600, 600)
        )
        
        plot!(plt, trace_x2_m[1:i], trace_y2_m[1:i], lw=1, color=:blue, alpha=0.5, label="Tracé mesuré")
        plot!(plt, trace_x2_s[1:i], trace_y2_s[1:i], lw=1, color=:red, alpha=0.5, label="Tracé simulé")
        
        plot!(plt, [0, x1_m, x2_m], [0, y1_m, y2_m], lw=3, color=:blue, label="Mesuré", marker=:circle, markersize=8)
        plot!(plt, [0, x1_s, x2_s], [0, y1_s, y2_s], lw=3, color=:red, ls=:dash, label="Simulé", marker=:circle, markersize=8)
        
        scatter!(plt, [0], [0], color=:black, markersize=10, label="Pivot")
    end
    
    gif(anim, "assets/pendule_animation.gif", fps=fps)
    println("✓ Animation sauvegardée: pendule_animation.gif")
    
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
    println("✓ Tracé final sauvegardé: pendule_trace.png")
end

function animate_extrapolation(t_sim, θ1_sim, θ2_sim, L1, L2, t_data_end; fps=30, skip=2)
    function pendulum_positions(θ1, θ2, L1, L2)
        x1 = L1 * sin(θ1)
        y1 = -L1 * cos(θ1)
        x2 = x1 + L2 * sin(θ2)
        y2 = y1 - L2 * cos(θ2)
        return x1, y1, x2, y2
    end
    
    L_total = (L1 + L2) * 1.3
    
    trace_x2 = Float64[]
    trace_y2 = Float64[]
    
    for j in 1:length(t_sim)
        _, _, x2, y2 = pendulum_positions(θ1_sim[j], θ2_sim[j], L1, L2)
        push!(trace_x2, x2)
        push!(trace_y2, y2)
    end
    
    idx_data_end = findfirst(t -> t >= t_data_end, t_sim)
    if idx_data_end === nothing
        idx_data_end = length(t_sim)
    end
    
    println("Création de l'animation d'extrapolation...")
    
    anim = @animate for i in 1:skip:length(t_sim)
        x1, y1, x2, y2 = pendulum_positions(θ1_sim[i], θ2_sim[i], L1, L2)
        
        in_extrapolation = i > idx_data_end
        zone_text = in_extrapolation ? " [EXTRAPOLATION]" : ""
        title_color = in_extrapolation ? :red : :blue
        
        plt = plot(
            xlim=(-L_total, L_total),
            ylim=(-L_total, L_total),
            aspect_ratio=:equal,
            legend=:topright,
            title="t = $(round(t_sim[i], digits=2)) s" * zone_text,
            titlefontcolor=title_color,
            xlabel="x (m)",
            ylabel="y (m)",
            size=(600, 600)
        )
        
        if idx_data_end > 1
            plot!(plt, trace_x2[1:min(i, idx_data_end)], trace_y2[1:min(i, idx_data_end)], 
                  lw=2, color=:blue, alpha=0.7, label="Données réelles")
        end
        
        if i > idx_data_end
            plot!(plt, trace_x2[idx_data_end:i], trace_y2[idx_data_end:i], 
                  lw=2, color=:red, alpha=0.7, label="Extrapolation")
        end
        
        pendulum_color = in_extrapolation ? :red : :blue
        plot!(plt, [0, x1, x2], [0, y1, y2], lw=3, color=pendulum_color, 
              label="Simulé", marker=:circle, markersize=8)
        
        scatter!(plt, [0], [0], color=:black, markersize=10, label="Pivot")
        
        if i > idx_data_end
            scatter!(plt, [trace_x2[idx_data_end]], [trace_y2[idx_data_end]], 
                    color=:green, markersize=12, marker=:star, 
                    label="Fin données")
        end
        
        if in_extrapolation
            annotate!(plt, 0, L_total*0.9, 
                     text("⚠ Prédiction du modèle", 10, :center, :red))
        end
    end
    
    gif(anim, "assets/pendule_extrapolation.gif", fps=fps)
    println("✓ Animation extrapolation sauvegardée: pendule_extrapolation.gif")
end

# =============================================================================
# Statistiques
# =============================================================================

println("\n" * "="^60)
println("STATISTIQUES")
println("="^60)

R2_θ1 = calculate_R2(θ1_data, θ1_sim)
R2_θ2 = calculate_R2(θ2_data, θ2_sim)
RMSE_θ1 = calculate_RMSE(θ1_data, θ1_sim)
RMSE_θ2 = calculate_RMSE(θ2_data, θ2_sim)
MAE_θ1 = calculate_MAE(θ1_data, θ1_sim)
MAE_θ2 = calculate_MAE(θ2_data, θ2_sim)

println("\nMétriques simulation/mesures:")
println("  R²: θ1=$(round(R2_θ1, digits=4)), θ2=$(round(R2_θ2, digits=4))")
println("  RMSE: θ1=$(round(rad2deg(RMSE_θ1), digits=2))°, θ2=$(round(rad2deg(RMSE_θ2), digits=2))°")
println("  MAE: θ1=$(round(rad2deg(MAE_θ1), digits=2))°, θ2=$(round(rad2deg(MAE_θ2), digits=2))°")

println("\nExposant de Lyapunov:")
λ, t_lyap, separations = lyapunov_exponent(m_opt[1], m_opt[2], L1, L2, g, -θ1_0, θ2_0, 2.0)
println("  λ = $(round(λ, digits=3)) s⁻¹")
if λ > 0
    println("  → Système CHAOTIQUE (temps de prédictibilité: $(round(1/λ, digits=2))s)")
else
    println("  → Système régulier")
end

println("\nÉnergie (simulation):")
E_total, E_kinetic, E_potential = energy_over_time(t_data, sol_opt, m_opt[1], m_opt[2], L1, L2, g)
E_std = std(E_total)
E_std_relative = E_std / abs(mean(E_total)) * 100
println("  E initiale: $(round(E_total[1], digits=6)) J")
println("  σ(E): $(round(E_std, sigdigits=3)) J ($(round(E_std_relative, digits=3))%)")

println("\nDivergence mesures/simulation:")
div_analysis = analyze_divergence(t_data, θ1_data, θ2_data, θ1_sim, θ2_sim, threshold_deg=10.0)
if !isnan(div_analysis.t_divergence)
    println("  ⚠ Divergence à t = $(round(div_analysis.t_divergence, digits=2))s")
else
    println("  ✓ Pas de divergence >10°")
end

t_valid, idx_valid = find_best_fit_window(t_data, θ1_data, θ2_data, θ1_sim, θ2_sim, max_error_deg=5.0)
println("  Fenêtre valide (erreur <5°): 0-$(round(t_valid, digits=2))s ($(round(100*idx_valid/length(t_data), digits=1))%)")

# =============================================================================
# Export CSV
# =============================================================================

df_comparison = DataFrame(
    temps_s = t_data,
    theta1_mesure_deg = rad2deg.(θ1_data),
    theta1_simule_deg = rad2deg.(θ1_sim),
    theta2_mesure_deg = rad2deg.(θ2_data),
    theta2_simule_deg = rad2deg.(θ2_sim),
    erreur_theta1_deg = rad2deg.(abs.(θ1_data .- θ1_sim)),
    erreur_theta2_deg = rad2deg.(abs.(θ2_data .- θ2_sim))
)

CSV.write("assets/comparison_data.csv", df_comparison)
println("\n✓ Données exportées: comparison_data.csv")

# =============================================================================
# Graphiques
# =============================================================================

p_angles = plot(t_data, rad2deg.(θ1_data), label="θ1 mesuré", lw=2, title="Comparaison")
plot!(p_angles, t_data, rad2deg.(θ1_sim), label="θ1 simulé", ls=:dash, lw=2)
plot!(p_angles, t_data, rad2deg.(θ2_data), label="θ2 mesuré", lw=2)
plot!(p_angles, t_data, rad2deg.(θ2_sim), label="θ2 simulé", ls=:dash, lw=2)
xlabel!(p_angles, "Temps (s)")
ylabel!(p_angles, "Angle (°)")
savefig(p_angles, "assets/pendule_comparison.png")

p_energy = plot(t_data, E_total, label="E totale", lw=2, title="Énergie")
plot!(p_energy, t_data, E_kinetic, label="E cinétique", lw=2)
plot!(p_energy, t_data, E_potential, label="E potentielle", lw=2)
E_mean = mean(E_total)
annotate!(p_energy, t_data[end]*0.7, maximum(E_total)*0.95,
         text("σ(E) = $(round(E_std, sigdigits=3)) J ($(round(E_std_relative, digits=3))%)", 
              9, :left, :black))
xlabel!(p_energy, "Temps (s)")
ylabel!(p_energy, "Énergie (J)")
savefig(p_energy, "assets/energie.png")

p_div = plot(t_data, rad2deg.(div_analysis.errors_θ1), label="Erreur θ1", lw=2, title="Divergence")
plot!(p_div, t_data, rad2deg.(div_analysis.errors_θ2), label="Erreur θ2", lw=2)
plot!(p_div, t_data, rad2deg.(div_analysis.error_total), label="Erreur totale", lw=2, color=:red)
hline!(p_div, [10.0], label="Seuil 10°", ls=:dash, color=:black)
xlabel!(p_div, "Temps (s)")
ylabel!(p_div, "Erreur (°)")
savefig(p_div, "assets/divergence.png")

println("\n✓ Graphiques sauvegardés")

animate_pendulum(t_data, θ1_data, θ2_data, θ1_sim, θ2_sim, L1, L2, fps=10, skip=2)

# Frame 1 avec masses et longueurs
function save_frame_detailed(θ1_data, θ2_data, θ1_sim, θ2_sim, L1, L2, m1, m2, t_data; frame=1)
    function pendulum_positions(θ1, θ2, L1, L2)
        x1 = L1 * sin(θ1)
        y1 = -L1 * cos(θ1)
        x2 = x1 + L2 * sin(θ2)
        y2 = y1 - L2 * cos(θ2)
        return x1, y1, x2, y2
    end
    
    L_total = (L1 + L2) * 1.3
    
    x1_m, y1_m, x2_m, y2_m = pendulum_positions(θ1_data[frame], θ2_data[frame], L1, L2)
    x1_s, y1_s, x2_s, y2_s = pendulum_positions(θ1_sim[frame], θ2_sim[frame], L1, L2)
    
    plt = plot(
        xlim=(-L_total, L_total),
        ylim=(-L_total, L_total),
        aspect_ratio=:equal,
        legend=:topright,
        title="Pendule Double - t = $(round(t_data[frame], digits=3)) s",
        xlabel="x (m)",
        ylabel="y (m)",
        size=(700, 700)
    )
    
    plot!(plt, [0, x1_m, x2_m], [0, y1_m, y2_m], lw=3, color=:blue, label="Mesuré", marker=:circle, markersize=10)
    plot!(plt, [0, x1_s, x2_s], [0, y1_s, y2_s], lw=3, color=:red, ls=:dash, label="Simulé", marker=:circle, markersize=10)
    
    scatter!(plt, [0], [0], color=:black, markersize=12, label="Pivot")
    
    # Annotations des masses et longueurs
    #= annotate!(plt, x1_m*0.5, y1_m*0.5, text("L₁ = $(round(L1*1000, digits=1)) mm", 9, :blue))
    annotate!(plt, (x1_m+x2_m)*0.5, (y1_m+y2_m)*0.5, text("L₂ = $(round(L2*1000, digits=1)) mm", 9, :blue))
    
    annotate!(plt, x1_m, y1_m - 0.02, text("m₁ = $(round(m1*1000, digits=1)) g", 9, :blue, :bottom))
    annotate!(plt, x2_m, y2_m - 0.02, text("m₂ = $(round(m2*1000, digits=1)) g", 9, :blue, :bottom)) =#
    
    # Paramètres système
    annotate!(plt, -L_total*0.95, L_total*0.95, 
             text("Paramètres système:", 8, :left, :black, :bold))
    annotate!(plt, -L_total*0.95, L_total*0.88, 
             text("m₁ = $(round(m1*1000, digits=1)) g", 8, :left))
    annotate!(plt, -L_total*0.95, L_total*0.82, 
             text("m₂ = $(round(m2*1000, digits=1)) g", 8, :left))
    annotate!(plt, -L_total*0.95, L_total*0.76, 
             text("L₁ = $(round(L1*1000, digits=1)) mm", 8, :left))
    annotate!(plt, -L_total*0.95, L_total*0.70, 
             text("L₂ = $(round(L2*1000, digits=1)) mm", 8, :left))
    annotate!(plt, -L_total*0.95, L_total*0.64, 
             text("(m₁+m₂)/m₁ = $(round((m1+m2)/m1, digits=2))", 8, :left))

    δ_deg = rad2deg(θ2_data[frame] - θ1_data[frame])
    annotate!(plt, -L_total*0.95, L_total*0.58, 
             text("δ = θ₂ - θ₁ = $(round(δ_deg, digits=1))°", 8, :left))
    
    savefig(plt, "assets/pendule_frame1.png")
    println("✓ Frame détaillée sauvegardée: pendule_frame1.png")
end

save_frame_detailed(θ1_data, θ2_data, θ1_sim, θ2_sim, L1, L2, m1_opt, m2_opt, t_data, frame=1)

if @isdefined(θ1_sim_full) && @isdefined(θ2_sim_full) && @isdefined(t_sim_full)
    println("\nVisualisation extrapolation...")
    
    p_extrap = plot(title="Simulation avec extrapolation", size=(900, 500))
    plot!(p_extrap, t_data, rad2deg.(θ1_data), label="θ1 mesuré", lw=2, color=:blue)
    plot!(p_extrap, t_data, rad2deg.(θ2_data), label="θ2 mesuré", lw=2, color=:orange)
    plot!(p_extrap, t_sim_full, rad2deg.(θ1_sim_full), label="θ1 simulé", ls=:dash, lw=2, color=:blue, alpha=0.7)
    plot!(p_extrap, t_sim_full, rad2deg.(θ2_sim_full), label="θ2 simulé", ls=:dash, lw=2, color=:orange, alpha=0.7)
    vline!(p_extrap, [t_data[end]], label="Fin données", ls=:dot, color=:red, lw=2)
    vspan!(p_extrap, [t_data[end], t_sim_full[end]], alpha=0.1, color=:gray, label="Extrapolation")
    xlabel!(p_extrap, "Temps (s)")
    ylabel!(p_extrap, "Angle (°)")
    savefig(p_extrap, "assets/extrapolation.png")
    
    animate_extrapolation(t_sim_full, θ1_sim_full, θ2_sim_full, L1, L2, t_data[end], fps=20, skip=2)
    
    println("✓ Graphiques extrapolation sauvegardés")
end

println("\n" * "="^60)
println("✓ Analyse terminée!")
println("="^60)