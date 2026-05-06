using CairoMakie
using DelimitedFiles
using MathTeXEngine


fidel_interval = readdlm("fidel_interval_n14_30.txt")

fidelo = readdlm("fidelo_n14_30.txt")
fidel_lino = readdlm("fidel_lino_n14_30.txt")
fidele = readdlm("fidele_n14_30.txt")
fidel_line = readdlm("fidel_line_n14_30.txt")

fig1 = Figure(size=(640,480))

set_theme!(fonts = (; regular = "Computer Modern Roman", bold = "Computer Modern Bold", italic = "Computer Modern Italic", bold_italic = "Computer Modern Bold Italic"), fontsize = 20)

ax1 = Axis(fig1[1,1], xlabel=L"$\omega \tau$", ylabel=L"$\mathcal{F}$", limits = ((0.0,30.0),(0.0,1.0)), xticks = [0.0:5.0:30.0;], yticks = [0.0:0.2:1.0;])
lines!(ax1,[Point2f(fidel_interval[i], fidelo[i]) for i in 1:length(fidelo)], color = :slateblue, linewidth = 2.5, label = L"\text{MA-odd}")
lines!(ax1,[Point2f(fidel_interval[i], fidele[i]) for i in 1:length(fidele)], color = :crimson, linewidth = 2, label = L"\text{MA-even}")
lines!(ax1,[Point2f(fidel_interval[i], fidel_lino[i]) for i in 1:length(fidel_lino)], color = :slateblue, linestyle = :dash, linewidth = 2, label = L"\text{LR-odd}")
lines!(ax1,[Point2f(fidel_interval[i], fidel_line[i]) for i in 1:length(fidel_line)], color = :crimson, linestyle = :dash, linewidth = 2, label = L"\text{LR-even}")

axislegend(L"\textbf{Exact diagonalization}", position = :rb, gridshalign = :left)
fig1


save("wleg_fidel_numerical_even_odd_n14_30.pdf", fig1)