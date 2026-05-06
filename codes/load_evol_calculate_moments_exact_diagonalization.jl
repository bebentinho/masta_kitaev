using LinearAlgebra
using QuantumOptics
using JLD2
using ProgressMeter

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

I_N = tensor(fill(id, npart)...)

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

# fidel_interval = [vcat(LinRange(0.01,(0.1 - 0.01),9)); vcat(LinRange(0.1, 1.0 - 0.1, 9)); vcat(LinRange(1.0, 10.0 - 1.0, 9)); vcat(LinRange(10.0, 100.0, 10))]
fidel_interval = [vcat(LinRange(20.0, 100.0 - 10.0, 8)); vcat(LinRange(100.0, 1000.0, 10))]

@load "array_evole_N$(npart)_1704.jld2" array_evole
@load "array_evolo_N$(npart)_1704.jld2" array_evolo
@load "array_evol_line_N$(npart)_1704.jld2" array_evol_line
@load "array_evol_lino_N$(npart)_1704.jld2" array_evol_lino

# ==========================================
# 3. Calculate Moments (Memory-Safe Trick)
# ==========================================
function calculate_exact_moments(state, H_final, E_0)
    # 1. Mean Work: <W> = <H_f> - E_0
    E_final = real(expect(H_final, state))
    mean_W = E_final - E_0
    
    # 2. Build the Central Shifted Operator: W_c = H_f - <H_f> * I
    # This automatically recenters the distribution!
    W_c = H_final - E_final * I_N
    
    # 3. Successive multiplication (avoids calculating dense H^4 matrices)
    # W_c |Psi>
    state_c1 = W_c * state
    # W_c^2 |Psi>
    state_c2 = W_c * state_c1
    # W_c^3 |Psi>
    state_c3 = W_c * state_c2
    
    # 4. Central Moments via Inner Products: <Psi | W_c^n | Psi>
    # Note: dagger(state) * ket returns the complex overlap
    var_W = real(dagger(state) * state_c2)   # 2nd Central Moment
    m3_W  = real(dagger(state) * state_c3)   # 3rd Central Moment (Numerator of Skewness)
    
    # For M4, we just overlap state_c2 with itself: <Psi|W^2 * W^2|Psi>
    m4_W  = real(dagger(state_c2) * state_c2) # 4th Central Moment (Numerator of Kurtosis)
    
    return (media = mean_W, variancia = var_W, m3 = m3_W, m4 = m4_W)
end


H_i_MA = H_t(0, gvars, npart)
H_i_LR = H_lin_t(0, gvars_lin, npart)

# E_0 is just the ground state energy of the initial Hamiltonians
# E0_MA = real(eigenstates(dense(H_i_MA), 2)[1][1])
# E0_LR = real(eigenstates(dense(H_i_LR), 2)[1][1])
E0_MA = -13.000000000000002
E0_LR = -13.000000000000002

# Initialize Arrays for Results
results_even_MA = []
results_even_LR = []
results_odd_MA = []
results_odd_LR = []

@showprogress for (i, taux) in enumerate(fidel_interval)
    
    gvarsaux = Dict("alpha" => alpha_gvars, "beta" => beta_gvars, "gamma" => gamma_gvars, "tau" => taux)
    gvars_linaux = Dict("gf" => g_f, "tau" => taux)
    
    # 1. Build the Final Hamiltonians at t = tau
    Hf_MA = H_t(taux, gvarsaux, npart)
    Hf_LR = H_lin_t(taux, gvars_linaux, npart)
    
    # 2. Get the previously saved evolved states
    state_e_MA = array_evole[i]
    state_e_LR = array_evol_line[i]
    state_o_MA = array_evolo[i]
    state_o_LR = array_evol_lino[i]

    # 3. Calculate Moments
    res_e_MA = calculate_exact_moments(state_e_MA, Hf_MA, E0_MA)
    res_e_LR = calculate_exact_moments(state_e_LR, Hf_LR, E0_LR)
    res_o_MA = calculate_exact_moments(state_o_MA, Hf_MA, E0_MA)
    res_o_LR = calculate_exact_moments(state_o_LR, Hf_LR, E0_LR)

    
    push!(results_even_MA, (tau=taux, stats=res_e_MA))
    push!(results_even_LR, (tau=taux, stats=res_e_LR))
    push!(results_odd_MA, (tau=taux, stats=res_o_MA))
    push!(results_odd_LR, (tau=taux, stats=res_o_LR))
end

@save "results_even_MA_N$(npart)_1704.jld2" results_even_MA
@save "results_odd_MA_N$(npart)_1704.jld2" results_odd_MA
@save "results_even_LR_N$(npart)_1704.jld2" results_even_LR
@save "results_odd_LR_N$(npart)_1704.jld2" results_odd_LR
