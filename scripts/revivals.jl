using QuantumOptics
using Plots

#= parameters =#
δ = 0.1
g = 1.0
α = 4.0

#= basis =#
b_fock = FockBasis(40)
b_spin = SpinBasis(1 // 2)

#= operators =#
a = destroy(b_fock) ⊗ one(b_spin)
σ₋ = one(b_fock) ⊗ sigmam(b_spin)

#= Hamiltonian =#
H = -δ * dagger(a) * a + g * (dagger(a) * σ₋ + a * dagger(σ₋))

#= initial state =#
ψ₀ = coherentstate(b_fock, α) ⊗ spindown(b_spin)

#= time steps =#
T = [0:0.01:35;]

#= Schordinger's equation =#
tout, ψt = timeevolution.schroedinger(T, ψ₀, H)

#= expected value for atomic excitation =#
excitation = expect(dagger(σ₋) * σ₋, ψt)

using Plots

plot1 = plot(T, real(excitation), ylims=(0, 1), xlabel="gt", ylabel="L⟨σ₋t σ₋⟩")

plot(plot1)

savefig("../images/Revivals-JC-model.png")