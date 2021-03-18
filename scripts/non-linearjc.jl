using QuantumOptics
using SparseArrays

#= Parameters =#
N = 10
g = 1.0
α = 4.0

#= f(aat) = aat^n =#
n = .5

b_fock = FockBasis(N)
b_spin = SpinBasis(1 // 2)

σp = one(b_fock) ⊗ sigmap(b_spin)
σm = one(b_fock) ⊗ sigmam(b_spin)

function f(x)
    (x^.5) * (x - 1)^(n)
end


#= Defining A = f(a at)a =# 
function createA(b::FockBasis)
    diag = complex.(f.(b.offset + 1.:b.N))
    data = spdiagm(1 => diag)
    SparseOperator(b, data)
end

A = createA(b_fock) ⊗ one(b_spin)

V = g * (A * σp + dagger(A) * σm)

ψ0 = coherentstate(b_fock, α) ⊗ spindown(b_spin)

T = [0:0.01:35;]

tout, ψt = timeevolution.schroedinger(T, ψ0, V)
excitation = expect(dagger(σm) * σm, ψt)

using Plots

plot1 = plot(T, real(excitation), ylims=(0, 1), xlabel="gt", ylabel="⟨ σ₊ σ₋ ⟩")
plot(plot1)

savefig("../images/non-linearjc.png")