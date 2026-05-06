# using Markdown
# using InteractiveUtils
using LinearAlgebra
using QuantumOptics
using CairoMakie
using ProgressMeter
using GeometryBasics
using JLD2
using DelimitedFiles


omega = 1 # Energy hopping
g_0 = 0
npart = Int(14)
tau = Int(30) # total time evolution
delta = omega
g_f = -3*omega


# Formatting parameters

amt_times = 200 # How many points we calculate
ts = LinRange(0, tau, amt_times)


tmin = 0.5
tmax = tau

fidel_interval = vcat(LinRange(tmin,tmax,amt_times))

# eigens_0, eigsts_0 = eigenstates(dense(H_t(0, gvars, npart)), 2)
eigen_0 = -13.000000000000002

basis = tensor(fill(SpinBasis(1/2), npart)...)

@load "matriz_evol_lino_indexes.jld2" matriz_evol_lino_indexes
@load "matriz_evol_line_indexes.jld2" matriz_evol_line_indexes
# @load "matriz_evol_line.jld2" matriz_evol_line
# @load "matriz_evol_lino.jld2" matriz_evol_lino

# array_evolo = [Ket(basis, matriz_evolo_indexes[:, i]) for i in 1:size(matriz_evolo_indexes, 2)]
# array_evole = [Ket(basis, matriz_evole_indexes[:, i]) for i in 1:size(matriz_evole_indexes, 2)]
array_evol_line = [Ket(basis, matriz_evol_line_indexes[:, i]) for i in 1:size(matriz_evol_line_indexes, 2)]
array_evol_lino = [Ket(basis, matriz_evol_lino_indexes[:, i]) for i in 1:size(matriz_evol_lino_indexes, 2)]


data_3d_even_MA = []
data_3d_odd_MA = []
# data_3d_even_LR = []
# data_3d_odd_LR = []

tau_list = []
W_even_list = []
P_even_list = []
W_odd_list = []
P_odd_list = []

# index_to_plot = [Int(amt_times/40):Int(amt_times/40):Int(amt_times);]
# index = index_to_plot[parse(Int, ARGS[1])]

index = parse(Int, ARGS[1])

τ = fidel_interval[index]
# Hf_MA = H_t(τ, gvars, npart)
# Hf_LR = H_lin_t(τ, gvars_lin, npart)

name = "ma_sta_stuff_correct/estates_f_$(index).jld2"
name_en = "ma_sta_stuff_correct/evals_f_$(index).jld2"
@load name estates_f
@load name_en evals_f

indexesss = [1, 3, 5, 10, 20, 30, 45, 70, 110, 200]
index_to_ten = findfirst(x -> x == index, indexesss)
work_vals = evals_f .- eigen_0
probabilities_even = [abs2(dagger(estates_f[n]) * array_evol_line[Int(index_to_ten)]) for n in 1:length(evals_f)]
probabilities_odd = [abs2(dagger(estates_f[n]) * array_evol_lino[Int(index_to_ten)]) for n in 1:length(evals_f)]


@save "LR_work_$(index).jld2" work_vals
@save "LR_prob_even_$(index).jld2" probabilities_even
@save "LR_prob_odd_$(index).jld2" probabilities_odd