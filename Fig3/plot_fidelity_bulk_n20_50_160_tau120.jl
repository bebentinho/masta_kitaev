using CairoMakie
using DelimitedFiles
using MathTeXEngine



fidel_interval = readdlm("fidel_interval.csv")

fidel_20 = readdlm("fidel_n20.csv")
fidel_lin_20 = readdlm("fidel_lin_n20.csv")

fidel_50 = readdlm("fidel_n50.csv")
fidel_lin_50 = readdlm("fidel_lin_n50.csv")

fidel_80 = readdlm("fidel_n80.csv")
fidel_lin_80 = readdlm("fidel_lin_n80.csv")


fig1 = Figure(size=(640,560))

set_theme!(fonts = (; regular = "Computer Modern Roman", bold = "Computer Modern Bold", italic = "Computer Modern Italic", bold_italic = "Computer Modern Bold Italic"), fontsize = 20)

# cor20 = colorant"#11A5B1"
cor20 = colorant"#D62728FF"
# cor50 = colorant"#C71585"
cor50 = colorant"#008000FF"
# cor80 =colorant"#E69F00"
cor80 = colorant"#1F77B4FF"

# ax1 = Axis(fig1[1,1], xlabel=L"\omega \tau", ylabel=L"\mathcal{F}", limits = ((0.0,120.0),(0.0,1.0)), yscale = Makie.PowerScale(1/2), yticks = [0.0:0.1:0.5; 1.0])
ax1 = Axis(fig1[2,1], xlabel=L"\omega \tau", ylabel=L"\mathcal{F}", limits = ((0.0,120.0),(0.0,1.0)), yticks = [0.0:0.2:1.0;])
lines!(ax1,[Point2f(fidel_interval[i], fidel_20[i]) for i in 1:length(fidel_20)], color = cor20, linewidth = 2, label = "2pMA, N=20")
lines!(ax1,[Point2f(fidel_interval[i], fidel_lin_20[i]) for i in 1:length(fidel_lin_20)], color = cor20, linestyle = :dash, linewidth = 2, label = "LR, N=20")
lines!(ax1,[Point2f(fidel_interval[i], fidel_50[i]) for i in 1:length(fidel_50)], color = cor50, linewidth = 2, label = "2pMA, N=50")
lines!(ax1,[Point2f(fidel_interval[i], fidel_lin_50[i]) for i in 1:length(fidel_lin_50)], color = cor50, linestyle = :dash, linewidth = 2, label = "LR, N=50")
lines!(ax1,[Point2f(fidel_interval[i], fidel_80[i]) for i in 1:length(fidel_80)], color = cor80, linewidth = 2, label = "2pMA, N=80")
lines!(ax1,[Point2f(fidel_interval[i], fidel_lin_80[i]) for i in 1:length(fidel_lin_80)], color = cor80, linestyle = :dash, linewidth = 2, label = "LR, N=80")

Legend(fig1[1,1], ax1, orientation = :horizontal, nbanks = 2)
fig1

save("wleg_2nMASTA_fidel_n205080_tau120.pdf", fig1)