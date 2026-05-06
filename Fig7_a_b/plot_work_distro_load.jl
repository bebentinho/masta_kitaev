# using Markdown
# using InteractiveUtils
using LinearAlgebra
using QuantumOptics
using CairoMakie
using ProgressMeter
using GeometryBasics
using JLD2
using DelimitedFiles
using Colors
using Plots

function criar_colormap_customizado(nome_cor_base)
    # 1. Transformar o nome ou símbolo na cor em si
    cor_base = parse(Colorant, nome_cor_base)
    
    # 2. Criar o tom 50% mais claro (misturando com branco)
    # O peso 0.5 define que será exatamente no meio do caminho entre a cor e o branco
    cor_clara = weighted_color_mean(0.45, colorant"white", cor_base)
    
    # 3. Criar o tom 50% mais escuro (misturando com preto)
    cor_escura = weighted_color_mean(0.75, colorant"black", cor_base)
    
    # 4. Construir o ColorGradient passando pelos 3 estágios
    # O cgrad cria uma transição suave: Clara -> Intermediária (Base) -> Escura
    return cgrad([cor_clara, cor_base, cor_escura])
end

function plot_work_evolution_3d_volumetric(all_tau_data, filename, colormap = :viridis)
    set_theme!(fonts = (; regular = "Computer Modern Roman", bold = "Computer Modern Bold", italic = "Computer Modern Italic", bold_italic = "Computer Modern Bold Italic"), fontsize = 20)
    fig = Figure(size = (1000, 800))
    ax = Axis3(fig[1, 1],
        xlabel = L"$W/\omega$",
        ylabel = L"$\omega \tau$",
        zlabel = L"$P(W/\omega)$",
        yticks = (log10.([0.5, 1, 3, 10, 30]), ["0.5", "1", "3", "10", "30"]),
        zticks = ([0.0:0.25:1;], ["0.0","0.25", "0.5", "0.75", "1.0"]),
        azimuth = 0.25 * π,
        elevation = 0.15 * π,
        # limits = ((-10,10),(0.0, log10.(tau+5.0)), (0,1)),
        xgridvisible = true, ygridvisible = true, zgridvisible = true,
        protrusions = (100, 100, 100, 100),
        xticklabelsize = 22.5, yticklabelsize = 22.5, xlabelsize = 22.5, ylabelsize = 22.5
    )
    full_grad = cgrad(colormap)
    # colors = cgrad(full_grad[0.5:2.0], length(all_tau_data), categorical = true)
    colors = cgrad(full_grad, length(all_tau_data), categorical = true)
    for (i, data) in enumerate(all_tau_data)
        τ_val, W, P = data
        
        # scatter!(ax, W, fill(τ_val, length(W)), P, 
                # color = colors[i], markersize = 5)

        # Criando "cilindros" em vez de linhas finas
        # O Rect3f cria uma caixa (prisma) que é mais visível
        for (w_val, p_val) in zip(W, P)
            
            largura = bin_size * 0.8 # Largura do cilindro (80% do tamanho do bin)
            profundidade = 0.02
            alpha = 0.6
            if p_val > 0.0001
                mesh!(ax, 
                        Rect3f(Vec3f(w_val - largura/2, log10.(τ_val) - profundidade/2, 0), # Origem
                            Vec3f(largura, profundidade, p_val)),              # Dimensões (largura, profundidade, altura)
                        color = (colors[i], alpha),
                        shading = false
                )
            end    
        end


    end

    save(filename, fig, dpi=600)
    return fig
end


meu_reds = criar_colormap_customizado(:crimson)
meu_blues = criar_colormap_customizado(:slateblue)

omega = 1 # Energy hopping
g_0 = 0
npart = Int(14)
tau = Int(30) # total time evolution
delta = omega
g_f = -3*omega

# Formatting parameters

amt_times = 200 # How many points we calculate

tmin = 0.5
tmax = tau

fidel_interval = vcat(LinRange(tmin,tmax,amt_times))


eigen_0 = -13.000000000000002


basis = tensor(fill(SpinBasis(1/2), npart)...)


all_indexes = [1:1:4;5:5:200;]


times = fidel_interval[all_indexes]
logtimes=log10.(times)
interval = (logtimes[end]-logtimes[1])/9

# finding evenly spaced points in log scale
log_index = [findfirst(x -> abs(x - (logtimes[1]+interval*i))<0.06, logtimes) for i in 0:8]
index_to_plot = all_indexes[log_index]
push!(index_to_plot, all_indexes[end])


data_3d_even_MA = []
data_3d_odd_MA = []
data_3d_even_LR = []
data_3d_odd_LR = []

tau_list = []
W_even_list = []
P_even_list = []
W_odd_list = []
P_odd_list = []
W_even_LR_list = []
P_even_LR_list = []
W_odd_LR_list = []
P_odd_LR_list = []

@showprogress for aux in index_to_plot

    index = Int(aux)

    τ = fidel_interval[index]
    
    name_prob_even = "distrodata_MA/prob_even_$(index).jld2"
    name_prob_odd = "distrodata_MA/prob_odd_$(index).jld2"
    name_work = "distrodata_MA/work_$(index).jld2"
    name_prob_even_LR = "distrodata_LR/prob_even_$(index).jld2"
    name_prob_odd_LR = "distrodata_LR/prob_odd_$(index).jld2"
    name_work_LR = "distrodata_LR/work_$(index).jld2"

    @load name_prob_odd probabilities_odd
    @load name_prob_even probabilities_even
    @load name_work work_vals
    probs_odd_MA = probabilities_odd
    probs_even_MA = probabilities_even
    work_vals_MA = work_vals
    @load name_prob_odd_LR probabilities_odd
    @load name_prob_even_LR probabilities_even
    @load name_work_LR work_vals

    
    push!(tau_list, τ)

    push!(W_even_list, work_vals_MA)
    push!(P_even_list, probs_even_MA)
    push!(W_odd_list, work_vals_MA)
    push!(P_odd_list, probs_odd_MA)
    push!(data_3d_even_MA, (τ, work_vals_MA, probs_even_MA))
    push!(data_3d_odd_MA, (τ, work_vals_MA, probs_odd_MA))

    push!(W_even_LR_list, work_vals)
    push!(P_even_LR_list, probabilities_even)
    push!(W_odd_LR_list, work_vals)
    push!(P_odd_LR_list, probabilities_odd)
    push!(data_3d_even_LR, (τ, work_vals, probabilities_even))
    push!(data_3d_odd_LR, (τ, work_vals, probabilities_odd))
end



work_gs = data_3d_even_MA[1][2][1]
work_first_es = data_3d_even_MA[1][2][2]
work_most_es = data_3d_even_MA[1][2][end]
bin_size = abs(work_first_es - work_gs)
work_interval = abs(work_most_es - work_gs)
num_bins = Int(ceil(work_interval / bin_size))
work_vals = [work_gs + bin_size * j for j in 0:(num_bins-1)]


data_3d_even_MA_binned = []
data_3d_odd_MA_binned = []
data_3d_even_LR_binned = []
data_3d_odd_LR_binned = []

for i in 1:length(data_3d_even_MA)
    point_even = data_3d_even_MA[i]
    point_odd = data_3d_odd_MA[i]
    
    p_even = zeros(num_bins)
    p_odd = zeros(num_bins)

    for j in 1:length(point_even[2])
        aux_index_even = Int(div((point_even[2][j] - work_gs), bin_size) + 1)
        p_even[aux_index_even] += point_even[3][j]

        aux_index_odd = Int(div((point_odd[2][j] - work_gs), bin_size) + 1)
        p_odd[aux_index_odd] += point_odd[3][j]
    end


    new_point_even = (point_even[1], work_vals, p_even)
    new_point_odd = (point_odd[1], work_vals, p_odd)
    push!(data_3d_even_MA_binned, new_point_even)
    push!(data_3d_odd_MA_binned, new_point_odd)
end


fig = plot_work_evolution_3d_volumetric(data_3d_even_MA_binned, "pw_3d_even_MA_binned.png", meu_reds)
fig = plot_work_evolution_3d_volumetric(data_3d_odd_MA_binned, "pw_3d_odd_MA_binned.png", meu_blues)
fig = plot_work_evolution_3d_volumetric(data_3d_even_LR_binned, "pw_3d_even_LR_binned.png", meu_reds)
fig = plot_work_evolution_3d_volumetric(data_3d_odd_LR_binned, "pw_3d_odd_LR_binned.png", meu_blues)