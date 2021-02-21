using QuantumOptics

#= first specifies the photon number cutoff =#
N_cutoff = 20

#= define the basis =#
b = FockBasis(N_cutoff)

#= create a state and operator respect to b =#
a = destroy(b)
at = create(b)

α = 1.5
ψ = coherentstate(b, α)

dψ = -1im * (a + at) * ψ

#= density operator with dagger function =#
ρ = ψ ⊗ dagger(ψ)

#= alternatively =#
ρ1 = tensor(ψ, dagger(ψ))

#= printing expected values =#
println("α = ", dagger(ψ) * a * ψ)

#= OR =#

println("α = ", expect(at, ψ))
