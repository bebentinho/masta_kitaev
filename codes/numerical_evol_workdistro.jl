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

alpha_gvars = (64 * abs(delta)^2 * (sin(pi/npart))^2)^(-1)
beta_gvars = -2 * omega * cos(pi/npart)
gamma_gvars = 2 * abs(delta) * sin(pi/npart)

gvars = Dict("alpha" => alpha_gvars, "beta" => beta_gvars, "gamma" => gamma_gvars, "tau" => tau)
gvars_lin = Dict("gf" => g_f, "tau" => tau)

# Formatting parameters

amt_times = 200 # How many points we calculate
ts = LinRange(0, tau, amt_times)

id = identityoperator(SpinBasis(1/2))
σˣ = sigmax(SpinBasis(1/2))
σʸ = sigmay(SpinBasis(1/2))
σᶻ = sigmaz(SpinBasis(1/2))

function g_t(t, gvars)
  # Minimal action solution
  gvars["beta"] + gvars["gamma"] * tan((1-t/gvars["tau"]) * atan((g_0 - gvars["beta"])/gvars["gamma"]) + t/gvars["tau"] * atan((g_f - gvars["beta"])/gvars["gamma"]) )
end

function g_lin_t(t, gvars)
    # Linear ramp solution
    g_0 + (gvars["gf"] - g_0) * t / gvars["tau"]
end

function ℙₑ(state)
	round(1 + expect(tensor(fill(σᶻ, npart)...), state))/2
end

function ℙₒ(state)
	round(1 - expect(tensor(fill(σᶻ, npart)...), state))/2
end

function Hstatic(N)
	# vector of operators: [σˣ, σˣ, id, ...]
    🍎 = [σˣ; σˣ; fill(id, N-2)]
    
    # vector of operators: [σˣ, id, ...]
    🍉 = [σʸ; σʸ; fill(id, N-2)]
    
    H = 0*tensor(fill(id, N)...)
    for i in 1:N-1
        # tensor multiply all operators
        H -= 0.5 * (omega+delta) * tensor(🍎...)
		H -= 0.5 * (omega-delta) * tensor(🍉...)
        # cyclic shift the operators
        🍎 = circshift(🍎,1)
		🍉 = circshift(🍉,1)
    end
	
	H
end

function Hdyn(N)
	# vector of operators
	🥭 = [σᶻ; fill(id, N-1)]
	
	H = 0*tensor(fill(id, N)...)
    for i in 1:N
        # tensor multiply all operators
        H -= 0.5 * tensor(🥭...)
        # cyclic shift the operators
        🥭 = circshift(🥭,1)
    end

	H
end

function H_t(t, gvars, N)
	Hstatic(N) + g_t(t, gvars)*Hdyn(N)
end

function H_lin_t(t, gvars, N)
	Hstatic(N) + g_lin_t(t, gvars)*Hdyn(N)
end


# Now let's calculate the fidelity

tmin = 0.5
tmax = tau

fidel_interval = vcat(LinRange(tmin,tmax,amt_times))

# eigens_0, eigsts_0 = eigenstates(dense(H_t(0, gvars, npart)), 2)
# eigen_0 = -13.000000000000002
# gs1_0 = eigsts_0[1]
# gs2_0 = eigsts_0[2]

# gse_0 = ℙₑ(gs1_0)*gs1_0 + ℙₑ(gs2_0)*gs2_0
# gso_0 = ℙₒ(gs1_0)*gs1_0 + ℙₒ(gs2_0)*gs2_0


eigenslin_0, eigstslin_0 = eigenstates(dense(H_lin_t(0, gvars_lin, npart)), 2)
gslin1_0 = eigstslin_0[1]
gslin2_0 = eigstslin_0[2]

gsline_0 = ℙₑ(gslin1_0)*gslin1_0 + ℙₑ(gslin2_0)*gslin2_0
gslino_0 = ℙₒ(gslin1_0)*gslin1_0 + ℙₒ(gslin2_0)*gslin2_0


# fidele = zeros(amt_times)
# fidelo = zeros(amt_times)
fidel_line = zeros(amt_times)
fidel_lino = zeros(amt_times)

# array_evole = []
# array_evolo = []
array_evol_line = []
array_evol_lino = []

@showprogress for index in 1:amt_times
	time = fidel_interval[index]
	
	gvarsaux = Dict("alpha" => alpha_gvars, "beta" => beta_gvars, "gamma" => gamma_gvars, "tau" => time)
	gvars_linaux = Dict("gf" => g_f, "tau" => time)
	
	# eigens_t, eigsts_t = eigenstates(H_t(time, gvarsaux, npart), 2, info=false)
	# gs_t = eigsts_t[1]
	# es_t = eigsts_t[2]
	
	eigenslin_t, eigstslin_t = eigenstates(H_lin_t(time, gvars_linaux, npart), 2, info=false)
	gs_lin_t = eigstslin_t[1]
	es_lin_t = eigstslin_t[2]
	
	# evole = timeevolution.schroedinger_dynamic([0,time], gse_0, (t, psi) -> H_t(t, gvarsaux, npart))
	# evolo = timeevolution.schroedinger_dynamic([0,time], gso_0, (t, psi) -> H_t(t, gvarsaux, npart))
	
	evol_line = timeevolution.schroedinger_dynamic([0,time], gsline_0, (t, psi) -> H_lin_t(t, gvars_linaux, npart))
	evol_lino = timeevolution.schroedinger_dynamic([0,time], gslino_0, (t, psi) -> H_lin_t(t, gvars_linaux, npart))

    # push!(array_evolo, evolo[2][2])
    # push!(array_evole, evole[2][2])
    push!(array_evol_lino, evol_lino[2][2])
    push!(array_evol_line, evol_line[2][2])

	# fidelo[index] = norm(dagger(es_t)*evolo[2][2])^2
	# fidele[index] = norm(dagger(gs_t)*evole[2][2])^2
	fidel_lino[index] = norm(dagger(es_lin_t)*evol_lino[2][2])^2
	fidel_line[index] = norm(dagger(gs_lin_t)*evol_line[2][2])^2
end

basis = tensor(fill(SpinBasis(1/2), npart)...)

@load "matriz_evolo.jld2" matriz_evolo
@load "matriz_evole.jld2" matriz_evole
# @load "matriz_evol_line.jld2" matriz_evol_line
# @load "matriz_evol_lino.jld2" matriz_evol_lino


array_evolo = [Ket(basis, matriz_evolo[:, i]) for i in 1:1:4]
array_evole = [Ket(basis, matriz_evole[:, i]) for i in 1:1:4]
# array_evolo = [Ket(basis, matriz_evolo_indexes[:, i]) for i in 1:size(matriz_evolo_indexes, 2)]
# array_evole = [Ket(basis, matriz_evole_indexes[:, i]) for i in 1:size(matriz_evole_indexes, 2)]
# array_evol_line = [Ket(basis, matriz_evol_line[:, i]) for i in 1:size(matriz_evol_line, 2)]
# array_evol_lino = [Ket(basis, matriz_evol_lino[:, i]) for i in 1:size(matriz_evol_lino, 2)]



indexesss = [1, 3, 5, 10, 20, 30, 45, 70, 110, 200]
array_evol_lino = [Ket(basis, matriz_evol_lino[:, i]) for i in indexesss]
array_evol_line = [Ket(basis, matriz_evol_line[:, i]) for i in indexesss]

matriz_evol_lino_indexes = hcat([k.data for k in array_evol_lino]...)
@save "matriz_evol_lino_indexes.jld2" matriz_evol_lino_indexes
matriz_evol_line_indexes = hcat([k.data for k in array_evol_line]...)
@save "matriz_evol_line_indexes.jld2" matriz_evol_line_indexes

