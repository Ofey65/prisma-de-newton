using QuantumOptics
using Plots

N_cutoff = 10

wc = 0.1
wa = 0.1
omega = 1.

#= Bases =#
b_fock = FockBasis(N_cutoff)
b_spin = SpinBasis(1 // 2)
b = b_fock ⊗ b_spin

a = destroy(b_fock)
at = create(b_fock)
n = number(b_fock)

sm = sigmam(b_spin)
sp = sigmap(b_spin)
sz = sigmaz(b_spin)

Hatom = wa * sz / 2
Hfield = wc * n
Hint = omega * (at ⊗ sm + a ⊗ sp)
H = one(b_fock) ⊗ Hatom + Hfield ⊗ one(b_spin) + Hint


# Initial state
α = 1.
Ψ0 = coherentstate(b_fock, α) ⊗ spindown(b_spin)

# Integration time
T = [0:0.1:20;]

# Schroedinger time evolution
tout, Ψt = timeevolution.schroedinger(T, Ψ0, H)

exp_n = real(expect(n ⊗ one(b_spin), Ψt))
exp_sz = real(expect(one(b_fock) ⊗ sz, Ψt))

#= Using Plots insted of PyPlot =#
using Plots

plot1 = plot(T, exp_n, ylims=(0, 2), xlabel="Τ", ylabel="L⟨n⟩")
plot2 = plot(T, exp_sz, ylims=(-1, 1), xlabel="Τ", ylabel="L⟨ σ₂ ⟩")

plot(plot1, plot2, layout=(2, 1))

savefig("../images/JC-model.png")