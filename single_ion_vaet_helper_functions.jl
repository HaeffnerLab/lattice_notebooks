using IonSim
using QuantumOptics
using StochasticDiffEq


νa = 0.911e6  # axial trap frequency
νr1 = 3e6  # first radial trap frequency
νr2 = 3e6  # second radial trap frequency
B = 4.21e-4  # magnitude of B-field in Tesla
Bhat = (x̂ + ẑ)/√2  # direction of B-field
dB = 2π * abs(zeeman_shift(1, 2, 5/2, -1/2) - zeeman_shift(1, 0, 1/2, -1/2))
n̄eq = 1e4
timescale = 1e-6

C = Ca40(["S-1/2", "D-1/2"])
single_mode_chain = LinearChain(
                                ions=[C], com_frequencies=(x=νr1, y=νr2, z=νa), 
                                vibrational_modes=(x=[], y=[], z=[1])
                            )

"""
    siv_three_lasers(;
                J=12e3, Δ=20e3, κ=7e3, νeff=22.5e3, N=5, rwa_cutoff=8e5, lamb_dicke_order=1,
                δB=0, δJ=0, correlated=true, δν=0, phase_drift=0, correlated_phase_drift=true
            )
"""
function siv_three_lasers(;
                            J=12e3, Δ=20e3, κ=7e3, νeff=22.5e3, N=5, 
                            rwa_cutoff=8e5, lamb_dicke_order=1, δB=0, δJ=0, correlated=true,
                            δν=0, phase_drift=0, correlated_phase_drift=true
                        )
    
    L1, L2, L3 = laser(), laser(), laser()
    T = trap(configuration=single_mode_chain, B=B, Bhat=Bhat, lasers=[L1, L2, L3])
    com = T.configuration.vibrational_modes.z[1]
    com.N = N

    # Set the laser frequencies
    δf = transition_frequency(T, 1, ("S-1/2", "D-1/2")) + dB * δB / 2π
    δϕ = phase_drift * randn()

    # L1 will be used to drive the carrier J⋅σ_ϕ = J⋅σ_x + Δ⋅σ_y
    L1.Δ = δf
    L1.k = ẑ
    L1.ϵ = x̂
    if !correlated_phase_drift && !(phase_drift == 0)
        δϕ = phase_drift * randn()
    end
    L1.ϕ = t -> atan(Δ/J) / 2π + δϕ * t

    # L2/L3 will be used to drive detuned 1ˢᵗ order sidebands
    L2.Δ = δf + com.ν - νeff
    L2.k = ẑ
    L2.ϵ = x̂
    if !correlated_phase_drift && !(phase_drift == 0)
        δϕ = phase_drift * randn()
    end
    L2.ϕ = t -> 3/4 + δϕ * t

    L3.Δ = δf - com.ν + νeff
    L3.k = ẑ
    L3.ϵ = x̂
    if !correlated_phase_drift && !(phase_drift == 0)
        δϕ = phase_drift * randn()
    end
    L3.ϕ = t -> 1/4 + δϕ * t

    η = abs(get_η(com, L2, C))
    if correlated
        δ = 1 + δJ * randn()
        Ω0 = sqrt(J^2 + Δ^2) * δ
        Ω1 = κ/η * δ  # Set κ with κ = ηΩ
    else
        Ω0 = sqrt(J^2 + Δ^2) * (1 + δJ * randn())
        Ω1 = κ/η * (1 + δJ * randn()) # Set κ with κ = ηΩ
    end

    Efield_from_rabi_frequency!(Ω0, T, 1, 1, ("S-1/2", "D-1/2"))
    E = Efield_from_rabi_frequency(Ω1, T, 2, 1, ("S-1/2", "D-1/2"))
    L2.E = E; L3.E = E
    com.ν = νa + δν * randn()

    h = hamiltonian(
            T, rwa_cutoff=rwa_cutoff, lamb_dicke_order=lamb_dicke_order, timescale=timescale
        ) 
    if phase_drift == 0
        return h
    else
        return h, t -> L1.ϕ(t) - atan(Δ/J) / 2π
    end 
end

"""
    siv_four_lasers(;
                J=12e3, Δ=20e3, κ=7e3, νeff=22.5e3, N=5, rwa_cutoff=8e5, lamb_dicke_order=1,
                δB=0, δJ=0, correlated=true, δν=0, phase_drift=0, correlated_phase_drift=true
            )
"""
function siv_four_lasers(;
                            J=12e3, Δ=20e3, κ=7e3, νeff=22.5e3, N=5, 
                            rwa_cutoff=8e5, lamb_dicke_order=1, δB=0, δJ=0, correlated=true,
                            δν=0, phase_drift=0, correlated_phase_drift=true
                        )
    
    L1, L2, L3, L4 = laser(), laser(), laser(), laser()
    T = trap(configuration=single_mode_chain, B=B, Bhat=Bhat, lasers=[L1, L2, L3, L4])
    com = T.configuration.vibrational_modes.z[1]
    com.N = N

    # Set the laser frequencies
    δf = transition_frequency(T, 1, ("S-1/2", "D-1/2")) + dB * δB / 2π
    δϕ = phase_drift * randn()

    # L1 will be used to drive the carrier: J⋅σ_x
    L1.Δ = δf
    L1.k = ẑ
    L1.ϵ = x̂
    if !correlated_phase_drift && !(phase_drift == 0)
        δϕ = phase_drift * randn()
    end
    L1.ϕ = t -> δϕ * t

    # L2 will be used to drive the carrier: Δ⋅σ_y
    L2.Δ = δf
    L2.k = ẑ
    L2.ϵ = x̂
    if !correlated_phase_drift && !(phase_drift == 0)
        δϕ = phase_drift * randn()
    end
    L2.ϕ = t -> 1/4 + δϕ * t

    # L3/L4 will be used to drive detuned 1ˢᵗ order sidebands
    L3.Δ = δf + com.ν - νeff
    L3.k = ẑ
    L3.ϵ = x̂
    if !correlated_phase_drift && !(phase_drift == 0)
        δϕ = phase_drift * randn()
    end
    L3.ϕ = t -> 3/4 + δϕ * t

    L4.Δ = δf - com.ν + νeff
    L4.k = ẑ
    L4.ϵ = x̂
    if !correlated_phase_drift && !(phase_drift == 0)
        δϕ = phase_drift * randn()
    end
    L4.ϕ = t -> 1/4 + δϕ * t

    η = abs(get_η(com, L3, C))
    if correlated
        δ = 1 + δJ * randn()
        J *= δ
        Δ *= δ
        Ω1 = κ/η * δ  # Set κ with κ = ηΩ
    else
        J *= (1 + δJ * randn())
        Δ *= (1 + δJ * randn())
        Ω1 = κ/η * (1 + δJ * randn()) # Set κ with κ = ηΩ
    end

    Efield_from_rabi_frequency!(J, T, 1, 1, ("S-1/2", "D-1/2"))
    Efield_from_rabi_frequency!(Δ, T, 2, 1, ("S-1/2", "D-1/2"))

    E = Efield_from_rabi_frequency(Ω1, T, 3, 1, ("S-1/2", "D-1/2"))
    L3.E = E; L4.E = E
    com.ν = νa + δν * randn()

    h = hamiltonian(
            T, rwa_cutoff=rwa_cutoff, lamb_dicke_order=lamb_dicke_order, timescale=timescale
        ) 
    if phase_drift == 0
        return h
    else
        return h, L1.ϕ
    end 
end


"""
    rotated_solve(tspan, h; fout=nothing, n̄=0, θ=π/2, ξ=nothing, θ2=nothing)

`n̄` is the mean initial occupation of the motional state (assumed to be a thermal state), by 
default this is the ground state. `θ` is the angle of rotation for the spin state. This is 
applied once at the beginning and once at the end. By default this is equal to π/2. If `ξ` is
not `nothing`, then the second rotation is performed about the angle atan(⟨σy⟩/⟨σx⟩)=ξ. 
Further, if `θ2` is not `nothing` then it equals the magnitude of the second rotation.
If `fout=nothing`, then the state vector (or density matrix) is returned, but after a rotation 
 that maps σy → σx.

 If `heating_rate` is not `nothing`, then its value corresponds to the approximate linear 
 heating rate (for small time) in quanta / second.
"""
function rotated_solve(
            tspan, h; fout=nothing, n̄=0, θ=π/2, ξ=nothing, θ2=nothing, heating_rate=nothing,
            phase_noise=nothing
        )
    com = single_mode_chain.vibrational_modes.z[1]
    Rx = DenseOperator(C.basis, [cos(θ/2) -im*sin(θ/2); -im*sin(θ/2) cos(θ/2)]) ⊗ one(com.basis)
    isnothing(θ2) ? θ₂ = θ : θ₂ = θ2
    function Rx2(ξ, θ₂)
        if isnothing(ξ)
            nx = 1; ny = 0
        else
            nx = sqrt(1 / (tan(ξ)^2 + 1)) * sign(cos(ξ))
            ny = (1 - nx^2) * sign(sin(ξ))
        end
        (
            DenseOperator(C.basis, [cos(θ₂/2) 0; 0 cos(θ₂/2)]) +
            nx * DenseOperator(C.basis, [0 -im*sin(θ₂/2); -im*sin(θ₂/2) 0]) +
            ny * DenseOperator(C.basis, [0 -sin(θ₂/2); sin(θ₂/2) 0])
        ) ⊗ one(com.basis)
    end
    if n̄ ≡ 0
        ψᵢ = Rx * (C["S-1/2"] ⊗ fockstate(com.basis, 0))
    else
        ψᵢ = Rx * (dm(C["S-1/2"]) ⊗ thermalstate(com.basis, n̄)) * dagger(Rx)
    end
    if isnothing(fout)
        if isnothing(heating_rate)
            _, sol = timeevolution.schroedinger_dynamic(tspan, ψᵢ, h)
        else
            J(t, ψ) = (h(t, ψ), 
                            [
                                one(C.basis) ⊗ create(com.basis), 
                                one(C.basis) ⊗ destroy(com.basis)
                        ], 
                            [
                                one(C.basis) ⊗ destroy(com.basis), 
                                one(C.basis) ⊗ create(com.basis)
                        ], 
                            [heating_rate * timescale * (1 + 1 / n̄eq), heating_rate * timescale]
                    )
            _, sol = timeevolution.master_dynamic(tspan, ψᵢ, J)
        end
        if n̄ ≡ 0 && isnothing(heating_rate)
            if typeof(ξ) <: Vector
                return [Rx2(ξ[i], θ₂) * sol[i] for i in 1:length(sol)]
            else
                return [Rx2(ξ, θ₂) * psi for psi in sol]
            end
        else
            if typeof(ξ) <: Vector
                return [(Rx2(ξ[i], θ₂) * sol[i]) * dagger(Rx2(ξ[i], θ₂)) for i in 1:length(sol)]
            else
                return [(Rx2(ξ, θ₂) * rho) * dagger(Rx2(ξ, θ₂)) for rho in sol]
            end
        end
    else
        if n̄ ≡ 0 && isnothing(heating_rate)
            if typeof(ξ) <: Vector
                sol = timeevolution.schroedinger_dynamic(tspan, ψᵢ, h)[2]
                return [fout(tspan[i], Rx2(ξ[i], θ₂) * sol[i]) for i in 1:length(sol)]
            else
                f_out(t, psi) = fout(t, Rx2(ξ, θ₂) * psi)
                return timeevolution.schroedinger_dynamic(tspan, ψᵢ, h, fout=f_out)[2]
            end
        else
            f_out1(t, rho) = fout(t, Rx2(ξ, θ₂) * rho * dagger(Rx2(ξ, θ₂)))
            if isnothing(heating_rate)
                if typeof(ξ) <: Vector
                    sol = timeevolution.schroedinger_dynamic(tspan, ψᵢ, h)[2]
                    return [fout(tspan[i], Rx2(ξ[i], θ₂) * sol[i]) for i in 1:length(sol)]
                else
                    return timeevolution.schroedinger_dynamic(tspan, ψᵢ, h, fout=f_out1)[2]
                end
            else
                J1(t, ψ) = (h(t, ψ), 
                            [
                                one(C.basis) ⊗ create(com.basis), 
                                one(C.basis) ⊗ destroy(com.basis)
                        ], 
                            [
                                one(C.basis) ⊗ destroy(com.basis), 
                                one(C.basis) ⊗ create(com.basis)
                        ], 
                            [heating_rate * timescale * (1 + 1 / n̄eq), heating_rate * timescale]
                    )
                if typeof(ξ) <: Vector
                    sol = timeevolution.master_dynamic(tspan, ψᵢ, J1)[2]
                    return [fout(tspan[i], Rx2(ξ[i], θ₂) * sol[i] * dagger(Rx2(ξ[i], θ₂))) for i in 1:length(sol)]
                else
                    return timeevolution.master_dynamic(tspan, ψᵢ, J1, fout=f_out1)[2]
                end
            end
        end
    end
end

function pfidelity(ρ, σ)
    l1 = length(ρ.basis_l.shape)
    l2 = length(σ.basis_l.shape)
    fidelity(ptrace(ρ, collect(2:l1)), ptrace(σ, collect(2:l2)))
end

function pfidelity(ρ::Vector, σ::Vector)
    v = zeros(Float64, length(ρ))
    for i in 1:length(ρ)
        v[i] = real(pfidelity(ρ[i], σ[i]))
    end
    v
end

gaussian_decay(t, γ) = (1 + exp(-(t * γ)^2/2))/2
exp_decay(t, γ) = (1 + exp(-t * γ))/2


# Ntraj = 250
# z_avg = zero(T)
# for i=1:Ntraj
#     t, z = stochastic.schroedinger(T, ψ0, H, Hs; fout=fout, alg=StochasticDiffEq.RKMil{:Stratonovich}())
#     z_avg .+= z./Ntraj
# end
