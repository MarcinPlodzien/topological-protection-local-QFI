#==============================================================================
mps_edge_qfi.jl -- edge-spin QFI in long interacting chains with MPS / TDVP

Tangent-state method: the theta-derivative of the state is evolved as a second
MPS, so no finite differences in theta are taken.  Model:

    H = -(J/2) sum_j [(1+gamma) X_j X_{j+1} + (1-gamma) Y_j Y_{j+1}]
        - h sum_j Z_j + delta sum_j Z_j Z_{j+1}        (PAULI matrices, OBC)

    |psi(theta)> = exp(-iHt) exp(-i theta sigma^a_k / 2) |dd...d>,  a = y, x

Conventions:
  * Pauli matrices in H,
  * initial state all DOWN, nearest-neighbour couplings only,
  * only the single-site QFI of site `j_read` (default 1) is evaluated,
  * both encoding axes from ONE tangent state: d|psi> = c |phi>, with
    |phi> = exp(-iHt)|flip at k>, c = -1/2 (R_y) and c = -i/2 (R_x).

Usage:
    julia mps_edge_qfi.jl L gamma h delta tmax dt maxdim [cutoff] [k] [outdir]
Output: <outdir>/mps_L{L}_g{gamma}_h{h}_d{delta}_chi{maxdim}_dt{dt}.txt with columns
    tJ  F_y  F_x  chi_psi  chi_phi  F_global  trunc_norm_phi
F_global = 4(<dpsi|dpsi> - |<psi|dpsi>|^2) must stay 1; its drift is the error monitor.
norm_phi is the accumulated norm of the tangent state lost by truncation (1 without truncation).

Method.  Two-site TDVP (ITensorMPS.tdvp, one sweep per step) for |psi(t)> = exp(-iHt)|dd...d> and for the
tangent state |phi(t)> = exp(-iHt)|flip at k>, both without conserved quantum numbers (the pairing term breaks
U(1)).  The reduced density matrix of the read-out site and the transition density matrix between |phi> and |psi>
give rho and d rho/d theta; the QFI follows from the symmetric logarithmic derivative in the eigenbasis of rho
(qfi_sld).  Site k and read-out site are 1-based here (k = 1 is the edge).

Runs used in the paper (outdir default ../data/mps):
  julia codes/mps_edge_qfi.jl 48 1.0 0.2 0.2 15.0 0.1 128 1e-9    stopped at tJ = 13 (FigR10b uses tJ <= 10)
  julia codes/mps_edge_qfi.jl 48 1.0 0.2 0.0 15.0 0.1 128 1e-9    complete to tJ = 15 (FigR10b, delta = 0)
  julia codes/mps_edge_qfi.jl 10 1.0 0.2 0.2 6.0 0.05 64 1e-12    L = 10 benchmark against full diagonalization (reply)
  julia codes/mps_edge_qfi.jl 48 1.0 0.2 0.4 15.0 0.1 128 1e-9    stopped at tJ = 3, single core (FigR12 of the reply)
  julia codes/mps_edge_qfi.jl 100 1.0 0.2 0.2 8.0 0.1 128 1e-9    complete to tJ = 8 (FigR12); OPENBLAS_NUM_THREADS=4,
                                                                  JULIA_NUM_THREADS=1, 7 h 46 min
Restart: the checkpoint <output>.ckpt is written every unit of time (every 10 steps of dt = 0.1).
Requires Julia with ITensors and ITensorMPS (versions in the README of the repository).
==============================================================================#
using ITensors, ITensorMPS
using LinearAlgebra, Printf, Serialization

# MPO of H in spin-1/2 operators S = sigma/2: XX + YY = 2(S+S- + S-S+), XX - YY = 2(S+S+ + S-S-),
# Z_i Z_{i+1} = 4 Sz Sz, Z = 2 Sz; the coefficients below reproduce the Pauli Hamiltonian of the paper exactly.
function build_H(s; L, J, gamma, h, delta)
    os = OpSum()
    for i in 1:(L - 1)
        # -(J/2)[(1+g)XX + (1-g)YY] = -J (S+S- + S-S+) - J g (S+S+ + S-S-)
        os += -J, "S+", i, "S-", i + 1
        os += -J, "S-", i, "S+", i + 1
        if gamma != 0.0
            os += -J * gamma, "S+", i, "S+", i + 1
            os += -J * gamma, "S-", i, "S-", i + 1
        end
        if delta != 0.0
            os += 4.0 * delta, "Sz", i, "Sz", i + 1                    # delta ZZ
        end
    end
    for i in 1:L
        os += -2.0 * h, "Sz", i                                         # -h Z
    end
    return MPO(os, s)
end

# tau[a,b] = <a| Tr_env(|ket><bra|) |b> on one site x (rows: ket, columns: bra)
function transition_rdm1(ket::MPS, bra::MPS, x::Int)
    s = siteinds(ket); N = length(ket)
    brap = [prime(dag(bra[j]), "Link") for j in 1:N]
    Lenv = ITensor(1.0 + 0.0im)
    for j in 1:(x - 1)
        Lenv = (Lenv * ket[j]) * brap[j]
    end
    Renv = ITensor(1.0 + 0.0im)
    for j in N:-1:(x + 1)
        Renv = (Renv * ket[j]) * brap[j]
    end
    T = ((Lenv * ket[x]) * prime(brap[x], s[x])) * Renv
    return Matrix{ComplexF64}(Array(T, s[x], prime(s[x])))
end

function qfi_sld(rho, drho; tol = 1e-12)
    rho = (rho + rho') / 2
    vals, vecs = eigen(Hermitian(rho)); vals = max.(vals, 0.0)
    P = vecs' * drho * vecs
    F = 0.0
    for i in eachindex(vals), j in eachindex(vals)
        vals[i] + vals[j] > tol && (F += 2 * abs2(P[i, j]) / (vals[i] + vals[j]))
    end
    return F
end

function tdvp_step(H, psi, dt; cutoff, maxdim)
    try
        return tdvp(H, -1im * dt, psi; nsweeps = 1, cutoff = cutoff, maxdim = maxdim, nsite = 2, outputlevel = 0)
    catch e
        e isa MethodError || e isa ErrorException || rethrow()
        return tdvp(H, -1im * dt, psi; nsweeps = 1, cutoff = cutoff, maxdim = maxdim, nsites = 2, outputlevel = 0)
    end
end

function main(args)
    L = parse(Int, args[1]); gamma, h, delta, tmax, dt = parse.(Float64, args[2:6])
    maxdim = parse(Int, args[7])
    cutoff = length(args) >= 8 ? parse(Float64, args[8]) : 1e-10
    k = length(args) >= 9 ? parse(Int, args[9]) : 1
    outdir = length(args) >= 10 ? args[10] : joinpath(@__DIR__, "..", "data", "mps")
    jread = 1; J = 1.0
    mkpath(outdir)
    fname = joinpath(outdir, @sprintf("mps_L%d_g%.2f_h%.2f_d%.2f_chi%d_dt%.3g.txt", L, gamma, h, delta, maxdim, dt))

    # Checkpoint: both MPS, the site indices and the rows written so far are saved after every unit of time
    # (file <output>.ckpt).  If it exists, the run continues from it; a finished run can be extended by
    # calling the script again with a larger tmax.
    ckpt = fname * ".ckpt"
    nsteps = round(Int, tmax / dt)
    if isfile(ckpt)
        step0, s, psi, phi, nphi, rows = deserialize(ckpt)
        println("resuming from $ckpt at t = $(step0 * dt)")
        if step0 >= nsteps
            println("nothing to do: checkpoint is at t = $(step0 * dt) >= tmax"); return
        end
    else
        s = siteinds("S=1/2", L; conserve_qns = false)          # pairing term breaks U(1)
        psi = MPS(ComplexF64, s, ["Dn" for _ in 1:L])
        phi = MPS(ComplexF64, s, [i == k ? "Up" : "Dn" for i in 1:L])
        nphi = 1.0                                               # truncation-induced norm of phi
        rows = String[]; step0 = -1
    end
    H = build_H(s; L = L, J = J, gamma = gamma, h = h, delta = delta)
    for step in (step0 + 1):nsteps
        if step > 0                                              # advance from step-1 to step
            psi = tdvp_step(H, psi, dt; cutoff = cutoff, maxdim = maxdim)
            phi = tdvp_step(H, phi, dt; cutoff = cutoff, maxdim = maxdim)
            normalize!(psi); nphi *= norm(phi); normalize!(phi)
        end
        t = step * dt
        rho = transition_rdm1(psi, psi, jread)
        tau = nphi * transition_rdm1(phi, psi, jread)       # tau = Tr_env |phi><psi|; d rho = c tau + h.c.
        Fy = qfi_sld(rho, (-0.5) * tau + ((-0.5) * tau)')
        Fx = qfi_sld(rho, (-0.5im) * tau + ((-0.5im) * tau)')
        Fg = 4 * (0.25 * nphi^2 * real(inner(phi, phi)) - abs2(0.5 * nphi * inner(psi, phi)))
        push!(rows, @sprintf("%.4f %.8e %.8e %d %d %.8f %.8f", t, Fy, Fx, maxlinkdim(psi), maxlinkdim(phi), Fg, nphi))
        if step % 10 == 0 || step == nsteps
            @printf("t=%7.2f  F_y=%.5f  F_x=%.5f  chi=%d/%d  F_global=%.6f\n", t, Fy, Fx, maxlinkdim(psi), maxlinkdim(phi), Fg)
            flush(stdout)
            open(fname, "w") do f
                println(f, "# MPS/TDVP(2-site): L=$L gamma=$gamma h=$h delta=$delta J=1 (Pauli), dt=$dt maxdim=$maxdim cutoff=$cutoff, rotation at site $k, read-out site $jread")
                println(f, "# tJ  F_y  F_x  chi_psi  chi_phi  F_global  norm_phi")
                foreach(r -> println(f, r), rows)
            end
            serialize(ckpt * ".tmp", (step, s, psi, phi, nphi, rows)); mv(ckpt * ".tmp", ckpt; force = true)
        end
    end
    println("wrote $fname")
end

main(ARGS)
