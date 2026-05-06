using LinearAlgebra
using QuantumOptics
using CairoMakie
using ProgressMeter
using GeometryBasics
using JLD2
using DelimitedFiles
using Colors


npart = Int(14)



@load "results_even_MA_N14.jld2" results_even_MA
@load "results_odd_MA_N14.jld2" results_odd_MA
@load "results_even_LR_N14.jld2" results_even_LR
@load "results_odd_LR_N14.jld2" results_odd_LR

results_even_MA = results_even_MA[1:30]
results_odd_MA = results_odd_MA[1:30]
results_even_LR = results_even_LR[1:30]
results_odd_LR = results_odd_LR[1:30]


set_theme!(fonts = (; regular = "Computer Modern Roman", bold = "Computer Modern Bold", italic = "Computer Modern Italic", bold_italic = "Computer Modern Bold Italic"), fontsize = 20)

# ==============================================================================


fig1 = Figure(size=(640,480))

ax1 = Axis(fig1[1,1], xlabel=L"\omega \tau", ylabel=L"\langle W \rangle / \omega",
xticklabelsize = 27.3, yticklabelsize = 27.3, xlabelsize = 27.3, ylabelsize = 27.3, yticks = -10:2:0)

lines!(ax1,[Point2f(results_odd_MA[i].tau, results_odd_MA[i].stats.media) for i in 1:length(results_odd_MA)], color = :slateblue, linewidth = 2.5, label = L"\text{MA-odd}")
lines!(ax1,[Point2f(results_even_MA[i].tau, results_even_MA[i].stats.media) for i in 1:length(results_even_MA)], color = :crimson, linewidth = 2, label = L"\text{MA-even}")
lines!(ax1,[Point2f(results_odd_LR[i].tau, results_odd_LR[i].stats.media) for i in 1:length(results_odd_LR)], color = :slateblue, linewidth = 2.5, linestyle = :dash, label = L"\text{LR-odd}")
lines!(ax1,[Point2f(results_even_LR[i].tau, results_even_LR[i].stats.media) for i in 1:length(results_even_LR)], color = :crimson, linewidth = 2, linestyle = :dash, label = L"\text{LR-even}")

ax1_inset = Axis(fig1[1,1], 
    width = Relative(0.38), height = Relative(0.35),
    # halign = 0.42,
    halign = 0.935, valign = 0.975, # Mudado para top-right como exemplo
    backgroundcolor = :white,
    xscale = log10, xticks = [0.01, 0.1, 1.0],
    xticklabelsize = 20, yticklabelsize = 20
)

translate!(ax1_inset.blockscene, 0, 0, 150)

lines!(ax1_inset,[Point2f(results_odd_MA[i].tau, results_odd_MA[i].stats.media) for i in 1:length(results_odd_MA)], color = :slateblue, linewidth = 2.5)
lines!(ax1_inset,[Point2f(results_even_MA[i].tau, results_even_MA[i].stats.media) for i in 1:length(results_even_MA)], color = :crimson, linewidth = 2)
lines!(ax1_inset,[Point2f(results_odd_LR[i].tau, results_odd_LR[i].stats.media) for i in 1:length(results_odd_LR)], color = :slateblue, linewidth = 2.5, linestyle = :dash)
lines!(ax1_inset,[Point2f(results_even_LR[i].tau, results_even_LR[i].stats.media) for i in 1:length(results_even_LR)], color = :crimson, linewidth = 2, linestyle = :dash)

xlims!(ax1_inset, 0.01, 1.0)
ylims!(ax1_inset, -2, nothing)

# axislegend(ax1, labelsize = 27.3)

fig1


# =================================================================



fig2 = Figure(size=(640,568))

ax2 = Axis(fig2[2,1], xlabel=L"\omega \tau", ylabel=L"\langle W - \bar{W} \rangle^2 / \omega^2",
xticklabelsize = 27.3, yticklabelsize = 27.3, xlabelsize = 27.3, ylabelsize = 27.3)

lines!(ax2,[Point2f(results_odd_MA[i].tau, results_odd_MA[i].stats.variancia) for i in 1:length(results_odd_MA)], color = :slateblue, linewidth = 2.5, label = L"\text{MA-odd}")
lines!(ax2,[Point2f(results_even_MA[i].tau, results_even_MA[i].stats.variancia) for i in 1:length(results_even_MA)], color = :crimson, linewidth = 2, label = L"\text{MA-even}")
lines!(ax2,[Point2f(results_odd_LR[i].tau, results_odd_LR[i].stats.variancia) for i in 1:length(results_odd_LR)], color = :slateblue, linewidth = 2.5, linestyle = :dash, label = L"\text{LR-odd}")
lines!(ax2,[Point2f(results_even_LR[i].tau, results_even_LR[i].stats.variancia) for i in 1:length(results_even_LR)], color = :crimson, linewidth = 2, linestyle = :dash, label = L"\text{LR-even}")

# --- INSET DA FIGURA 2 ---
ax2_inset = Axis(fig2[2,1], 
    width = Relative(0.38), height = Relative(0.35),
    #halign = 0.42,
    halign = 0.935, valign = 0.975, # Mudado para top-right como exemplo
    backgroundcolor = :white,
    xscale = log10, xticks = [0.01, 0.1, 1.0], yticks = [10.0, 20.0, 30.0],
    xticklabelsize = 20, yticklabelsize = 20
)

translate!(ax2_inset.blockscene, 0, 0, 150)

lines!(ax2_inset,[Point2f(results_odd_MA[i].tau, results_odd_MA[i].stats.variancia) for i in 1:length(results_odd_MA)], color = :slateblue, linewidth = 2.5)
lines!(ax2_inset,[Point2f(results_even_MA[i].tau, results_even_MA[i].stats.variancia) for i in 1:length(results_even_MA)], color = :crimson, linewidth = 2)
lines!(ax2_inset,[Point2f(results_odd_LR[i].tau, results_odd_LR[i].stats.variancia) for i in 1:length(results_odd_LR)], color = :slateblue, linewidth = 2.5, linestyle = :dash)
lines!(ax2_inset,[Point2f(results_even_LR[i].tau, results_even_LR[i].stats.variancia) for i in 1:length(results_even_LR)], color = :crimson, linewidth = 2, linestyle = :dash)

xlims!(ax2_inset, 0.01, 1.0)
ylims!(ax2_inset, 10, 35)


Legend(fig2[1,1], ax1, orientation = :horizontal, colgap = 40, nbanks = 2, labelsize = 27.3)
# axislegend(ax2, labelsize = 27.3)

fig2

# ==================================================================================

fig3 = Figure(size=(640,480))

ax3 = Axis(fig3[1,1], xlabel=L"\omega \tau", ylabel=L"\langle W - \bar{W} \rangle^3 / \omega^3",
xticklabelsize = 27.3, yticklabelsize = 27.3, xlabelsize = 27.3, ylabelsize = 27.3)

lines!(ax3,[Point2f(results_odd_MA[i].tau, results_odd_MA[i].stats.m3) for i in 1:length(results_odd_MA)], color = :slateblue, linewidth = 2.5, label = L"\text{MA-odd}")
lines!(ax3,[Point2f(results_even_MA[i].tau, results_even_MA[i].stats.m3) for i in 1:length(results_even_MA)], color = :crimson, linewidth = 2, label = L"\text{MA-even}")
lines!(ax3,[Point2f(results_odd_LR[i].tau, results_odd_LR[i].stats.m3) for i in 1:length(results_odd_LR)], color = :slateblue, linewidth = 2.5, linestyle = :dash, label = L"\text{LR-odd}")
lines!(ax3,[Point2f(results_even_LR[i].tau, results_even_LR[i].stats.m3) for i in 1:length(results_even_LR)], color = :crimson, linewidth = 2, linestyle = :dash, label = L"\text{LR-even}")

# --- INSET DA FIGURA 3 ---
ax3_inset = Axis(fig3[1,1], 
    width = Relative(0.38), height = Relative(0.35),
    #halign = 0.42,
    halign = 0.935, valign = 0.975, # Mudado para top-right como exemplo
    backgroundcolor = :white,
    xscale = log10, xticks = [0.01, 0.1, 1.0], yticks = [50.0, 100.0],
    xticklabelsize = 20, yticklabelsize = 20
)

translate!(ax3_inset.blockscene, 0, 0, 150)

lines!(ax3_inset,[Point2f(results_odd_MA[i].tau, results_odd_MA[i].stats.m3) for i in 1:length(results_odd_MA)], color = :slateblue, linewidth = 2.5)
lines!(ax3_inset,[Point2f(results_even_MA[i].tau, results_even_MA[i].stats.m3) for i in 1:length(results_even_MA)], color = :crimson, linewidth = 2)
lines!(ax3_inset,[Point2f(results_odd_LR[i].tau, results_odd_LR[i].stats.m3) for i in 1:length(results_odd_LR)], color = :slateblue, linewidth = 2.5, linestyle = :dash)
lines!(ax3_inset,[Point2f(results_even_LR[i].tau, results_even_LR[i].stats.m3) for i in 1:length(results_even_LR)], color = :crimson, linewidth = 2, linestyle = :dash)

xlims!(ax3_inset, 0.01, 1.0)
# ylims!(ax3_inset, 10, 35)

# axislegend(ax3, labelsize = 27.3)

fig3

# ===========================================================

# ax4 = Axis(fig4[1,1], xlabel=L"\omega \tau", ylabel=L"\langle W - \bar{W} \rangle^4 / \omega^4",
# xticklabelsize = 27.3, yticklabelsize = 27.3, xlabelsize = 27.3, ylabelsize = 27.3, xscale = log10)

# lines!(ax4,[Point2f(results_odd_MA[i].tau, results_odd_MA[i].stats.m4) for i in 1:length(results_odd_MA)], color = :slateblue, linewidth = 2.5, label = L"\text{MA-odd}")
# lines!(ax4,[Point2f(results_even_MA[i].tau, results_even_MA[i].stats.m4) for i in 1:length(results_even_MA)], color = :crimson, linewidth = 2, label = L"\text{MA-even}")
# lines!(ax4,[Point2f(results_odd_LR[i].tau, results_odd_LR[i].stats.m4) for i in 1:length(results_odd_LR)], color = :slateblue, linewidth = 2.5, linestyle = :dash, label = L"\text{LR-odd}")
# lines!(ax4,[Point2f(results_even_LR[i].tau, results_even_LR[i].stats.m4) for i in 1:length(results_even_LR)], color = :crimson, linewidth = 2, linestyle = :dash, label = L"\text{LR-even}")


# fig4

fig1
fig2
fig3
save("inset_mean_work_comparison.pdf", fig1)
save("inset_upleg_variance_work_comparison.pdf", fig2)
save("inset_skewness_work_comparison.pdf", fig3)