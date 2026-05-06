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


# fidel_interval = [vcat(LinRange(0.01,(0.1 - 0.01),9)); vcat(LinRange(0.1, 1.0 - 0.1, 9)); vcat(LinRange(1.0, 10.0 - 1.0, 9)); vcat(LinRange(10.0, 100.0, 10))]
fidel_interval = [vcat(LinRange(20.0, 100.0 - 10.0, 8)); vcat(LinRange(100.0, 1000.0, 10))]

eigens_0, eigsts_0 = eigenstates(dense(H_t(0, gvars, npart)), 2)
eigen_0 = eigens_0[1]
gs1_0 = eigsts_0[1]
gs2_0 = eigsts_0[2]

gse_0 = ℙₑ(gs1_0)*gs1_0 + ℙₑ(gs2_0)*gs2_0
gso_0 = ℙₒ(gs1_0)*gs1_0 + ℙₒ(gs2_0)*gs2_0


eigenslin_0, eigstslin_0 = eigenstates(dense(H_lin_t(0, gvars_lin, npart)), 2)
gslin1_0 = eigstslin_0[1]
gslin2_0 = eigstslin_0[2]

gsline_0 = ℙₑ(gslin1_0)*gslin1_0 + ℙₑ(gslin2_0)*gslin2_0
gslino_0 = ℙₒ(gslin1_0)*gslin1_0 + ℙₒ(gslin2_0)*gslin2_0



array_evole = []
array_evolo = []
array_evol_line = []
array_evol_lino = []


@showprogress for time in fidel_interval


    gvarsaux = Dict("alpha" => alpha_gvars, "beta" => beta_gvars, "gamma" => gamma_gvars, "tau" => time)
    gvars_linaux = Dict("gf" => g_f, "tau" => time)

    eigens_t, eigsts_t = eigenstates(H_t(time, gvarsaux, npart), 2, info=false)
    gs_t = eigsts_t[1]
    es_t = eigsts_t[2]

    eigenslin_t, eigstslin_t = eigenstates(H_lin_t(time, gvars_linaux, npart), 2, info=false)
    gs_lin_t = eigstslin_t[1]
    es_lin_t = eigstslin_t[2]

    evole = timeevolution.schroedinger_dynamic([0,time], gse_0, (t, psi) -> H_t(t, gvarsaux, npart),
            alg=Vern9(), reltol=1e-8, abstol=1e-10)
    evolo = timeevolution.schroedinger_dynamic([0,time], gso_0, (t, psi) -> H_t(t, gvarsaux, npart),
            alg=Vern9(), reltol=1e-8, abstol=1e-10)

    evol_line = timeevolution.schroedinger_dynamic([0,time], gsline_0, (t, psi) -> H_lin_t(t, gvars_linaux, npart),
            alg=Vern9(), reltol=1e-8, abstol=1e-10)
    evol_lino = timeevolution.schroedinger_dynamic([0,time], gslino_0, (t, psi) -> H_lin_t(t, gvars_linaux, npart),
            alg=Vern9(), reltol=1e-8, abstol=1e-10)

    push!(array_evolo, evolo[2][2])
    push!(array_evole, evole[2][2])
    push!(array_evol_lino, evol_lino[2][2])
    push!(array_evol_line, evol_line[2][2])


end

@save "array_evole_N$(npart)_2204.jld2" array_evole
@save "array_evolo_N$(npart)_2204.jld2" array_evolo
@save "array_evol_line_N$(npart)_2204.jld2" array_evol_line
@save "array_evol_lino_N$(npart)_2204.jld2" array_evol_lino