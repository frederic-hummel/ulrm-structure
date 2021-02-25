# import & setup of packages ___________________________________________________
using Interpolations, GSL, Calculus
# I/O paths ____________________________________________________________________
Path = "/home/frederic/Documents/julia/";
include("$Path/spin/Dimer\ Func.jl")

@time begin
# interactions _________________________________________________________________
    Cs = false
    swaveS = true
    swaveT = true
    pwaveS = true
    pwaveT = true
    pwaveSO = false
    fine = false
    hyperfine = false
    cfit = 1.0

# Phase shift data _____________________________________________________________
    (aS0, aS1, aP0, aP10, aP11, aP12) = scatdata(Cs, swaveS, swaveT, pwaveS, pwaveT, pwaveSO, cfit)

# Define parameters ____________________________________________________________
    n = 30 # energy level
    ρ = linspace(0.2n^2, 2.4n^2, 400);
    writecsv("$Path/pec/$(n)_grid.csv", ρ) 

    n_below = 1 # additional lower hydrogenic manifolds
    n_above = 0 # additional upper hydrogenic manifolds (one always implemented), if the closest manifold shall be neglected use L-dependent initialisation function
    m_max = 1.5
    L = 0 # orbital angular momentum
    (I1, α, α2) = Species(Cs) # nuclear spin, polarisability, ion polarizability
    I1 = 0
    Ahf = Hyperfine(Cs, hyperfine)
    δRb = Fine(Cs, fine, n)

    Basis = initializebasis(n, n_below, n_above, m_max, I1, δRb)

    offset = 1/2n^2 - eF(I1+0.5, I1, Ahf)
    kinetic(R, n) = real(sqrt(complex(2(-1/2n^2 + 1/R)))) # kinetic energy

    PES0 = Array{Float64}(length(ρ), size(bselect(Basis, 0.), 1))
    PES1 = Array{Float64}(length(ρ), size(bselect(Basis, 1.), 1))
    PES2 = Array{Float64}(length(ρ), size(bselect(Basis, 2.), 1))
end

for Ω0 in [2., 1., 0.]
# Define basis states __________________________________________________________
    basis = bselect(Basis, Ω0)
    Vbasis = Nstar(basis, δRb, α2, fine, false) # effective primary quantum number
    nlj = unique(basis[:,1:3], 1) # unperturbed basis states
    LSJ = LSJΩ(Ω0) # quantum numbers centered about perturber
    Ls = collect(0:1) # considered angular momentum
    Fs = ifelse(I1==0, 0.5, [I1-0.5, I1+0.5]) # hyperfine levels

# Set up constant matrices ___________________________________________________
    CG0 = CG(0, LSJ, basis)
    CGp = CG(1, LSJ, basis)
    CGm = CG(-1, LSJ, basis)

    V0 = spdiagm(-1./2Vbasis.^2 + FineStructure(basis)) + HF(basis, Ahf, Fs, I1)
    Delta = Hfilt(basis)

    Pnlj = Pos(nlj, basis)
    Pl = [1, 1, 2, 2, 2, 2]
    P_T = speye(size(basis)[1]) - ST(basis, 0, 0)

# Diagonalisation ______________________________________________________________
if Ω0 == 0.
    PES0 = Calc_pes(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α);
    PES0 = (PES0+offset)*n^3
    writecsv("$Path/pec/$(n)_Ω=0.csv", sort(PES0,2))
elseif Ω0 == 1.
    PES1 = Calc_pes(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α);
    PES1 = (PES1+offset)*n^3
    writecsv("$Path/pec/$(n)_Ω=1.csv", sort(PES1,2))
elseif Ω0 == 2.
    PES2 = Calc_pes(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α);
    PES2 = (PES2+offset)*n^3
    writecsv("$Path/pec/$(n)_Ω=2.csv", sort(PES2,2))
end
end

