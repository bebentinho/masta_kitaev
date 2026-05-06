using CairoMakie
using DelimitedFiles
using MathTeXEngine



fidel_interval = readdlm("2ndmom_n60_120_fidel_interval.txt")

fidelo = readdlm("n60_120_fidelo.txt")
fidel_lino = readdlm("n60_120_fidel_lino.txt")
fidele = readdlm("n60_120_fidele.txt")
fidel_line = readdlm("n60_120_fidel_line.txt")

fidelo_2nd = readdlm("2ndmom_n60_120_fidelo.txt")
fidel_lino_2nd = readdlm("2ndmom_n60_120_fidel_lino.txt")
fidele_2nd = readdlm("2ndmom_n60_120_fidele.txt")
fidel_line_2nd = readdlm("2ndmom_n60_120_fidel_line.txt")







fig1 = Figure(size=(640,480))

set_theme!(fonts = (; regular = "Computer Modern Roman", bold = "Computer Modern Bold", italic = "Computer Modern Italic", bold_italic = "Computer Modern Bold Italic"), fontsize = 20)

ax1 = Axis(fig1[1,1], xlabel=L"\omega \tau", ylabel=L"\mathcal{F}", limits = ((0.0,120.0),(0.0,1.0)), yscale = Makie.PowerScale(1/2), yticks = [0.0:0.1:0.5; 1.0])

lino_ma = lines!(ax1,[Point2f(fidel_interval[i], fidelo[i]) for i in 1:length(fidelo)], color = :slateblue, linewidth = 2.5, label = L"\text{MA-odd},\ k_{target}=\pi/N")
line_ma = lines!(ax1,[Point2f(fidel_interval[i], fidele[i]) for i in 1:length(fidele)], color = :crimson, linewidth = 2, label = L"\text{MA-even},\ k_{target}=\pi/N")

lino_lr = lines!(ax1,[Point2f(fidel_interval[i], fidel_lino[i]) for i in 1:length(fidel_lino)], color = :slateblue, linestyle = :dash, linewidth = 2, label = L"\text{LR-odd}")
line_lr = lines!(ax1,[Point2f(fidel_interval[i], fidel_line[i]) for i in 1:length(fidel_line)], color = :crimson, linestyle = :dash, linewidth = 2, label = L"\text{LR-even}")

axislegend(ax1, position = :rb)



fig2 = Figure(size=(640,480))

ax2 = Axis(fig2[1,1], xlabel=L"\omega \tau", ylabel=L"\mathcal{F}", limits = ((0.0,120.0),(0.0,1.0)), yscale = Makie.PowerScale(1/2), yticks = [0.0:0.1:0.5; 1.0])

lino_2nd_ma = lines!(ax2,[Point2f(fidel_interval[i], fidelo_2nd[i]) for i in 1:length(fidelo_2nd)], color = :slateblue, linestyle = :solid, alpha = 1.0, linewidth = 2, label = L"\text{MA-odd},\ k_{target}=2\pi/N")
line_2nd_ma = lines!(ax2,[Point2f(fidel_interval[i], fidele_2nd[i]) for i in 1:length(fidele_2nd)], color = :crimson, linestyle = :solid, alpha = 1.0, linewidth = 2, label = L"\text{MA-even},\ k_{target}=2\pi/N")

lino_2nd_lr = lines!(ax2,[Point2f(fidel_interval[i], fidel_lino_2nd[i]) for i in 1:length(fidel_lino_2nd)], color = :slateblue, linestyle = :dash, linewidth = 2, label = L"\text{LR-odd}")
line_2nd_lr = lines!(ax2,[Point2f(fidel_interval[i], fidel_line_2nd[i]) for i in 1:length(fidel_line_2nd)], color = :crimson, linestyle = :dash, linewidth = 2, label = L"\text{LR-even}")

axislegend(ax2, position = :rb)



save("wleg_fidel_bulk_even_odd_n60_120.pdf", fig1)
save("wleg_fidel_bulk_even_odd_2ndmom_n60_120.pdf", fig2)