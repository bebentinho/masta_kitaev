using CairoMakie
using DelimitedFiles
using MathTeXEngine
using LaTeXStrings



fidel_interval = readdlm("fidel_interval_comparison_n30.csv")

fidel_2nMA = readdlm("fidel_2nMA_n30.csv")
fidel_MA = readdlm("fidel_MA_n30.csv")
fidel_lin = readdlm("fidel_lin_n30.csv")

mu_2nMA = readdlm("mu_2nMA_n30.csv")
mu_MA = readdlm("mu_MA_n30.csv")
mu_lin = readdlm("mu_lin_n30.csv")

# cor20 = colorant"#11A5B1"
cor20 = colorant"#D62728FF"
# cor50 = colorant"#C71585"
cor50 = colorant"#008000FF"
# cor80 =colorant"#E69F00"
cor80 = colorant"#1F77B4FF"

fig1 = Figure(size=(640,480))
fig2 = Figure(size=(640,480))

set_theme!(fonts = (; regular = "Computer Modern Roman", bold = "Computer Modern Bold", italic = "Computer Modern Italic", bold_italic = "Computer Modern Bold Italic"), fontsize = 20)

ax1 = Axis(fig1[1,1], xlabel=L"\omega \tau", ylabel=L"\mathcal{F}", limits = ((0.0,120.0),(0.0,1.0)), yscale = Makie.PowerScale(1/2), yticks = [0.0:0.1:0.5; 1.0])
# ax1 = Axis(fig1[1,1], xlabel=L"\omega \tau", ylabel=L"\mathcal{F}", limits = ((0.0,30.0),(0.0,1.0)), xticks = [0.0:5.0:30.0;], yticks = [0.0:0.2:1.0;])
lines!(ax1,[Point2f(fidel_interval[i], fidel_2nMA[i]) for i in 1:length(fidel_2nMA)], color = :red, linewidth = 2, label = L"\text{2pMA}")
# lines!(ax1,[Point2f(fidel_interval[i], fidel_MA[i]) for i in 1:2], alpha = 0.0, label = " ")
lines!(ax1,[Point2f(fidel_interval[i], fidel_MA[i]) for i in 1:length(fidel_MA)], color = :black, linestyle = :dashdot, linewidth = 2.5, label = L"\text{MA}")
lines!(ax1,[Point2f(fidel_interval[i], fidel_lin[i]) for i in 1:length(fidel_lin)], color = :blue, linestyle = :dash, linewidth = 2, label = L"\text{LR}")

ax2 = Axis(fig2[1,1], xlabel=L"t/ \tau", ylabel=L"\mu(t)/\omega", xticks = [0.0:0.2:1.0;], yticks = [-3.0:1.0:3.0;])
lines!(ax2,[Point2f(fidel_interval[i] ./ fidel_interval[end], mu_2nMA[i]) for i in 1:length(mu_2nMA)], color = :red, linewidth = 2, label = L"\text{2pMA}")
# lines!(ax2,[Point2f(fidel_interval[i] ./ fidel_interval[end], mu_MA[i]) for i in 1:2], alpha = 0.0, label = " ")
lines!(ax2,[Point2f(fidel_interval[i] ./ fidel_interval[end], mu_MA[i]) for i in 1:length(mu_MA)], color = :black, linewidth = 2.5, linestyle = :dashdot, label = L"\text{MA}")
lines!(ax2,[Point2f(fidel_interval[i] ./ fidel_interval[end], mu_lin[i]) for i in 1:length(mu_lin)], color = :blue, linestyle = :dash, linewidth = 2, label = L"\text{LR}")

axislegend(ax1, position = :rc)
fig1
axislegend(ax2, position = :lt)
fig2

save("wleg_2n_1_MASTA_fidel_n30_tau120.pdf", fig1)
save("wleg_2n_1_MASTA_g_t_n30_tau120.pdf", fig2)
