println("prepare packages...")
using Interpolations, GSL, Calculus
Path = "/home/frederic/Documents/julia/";
include("$Path/field/nlm\ field\ Func.jl")

# Phase shift data _____________________________________________________________
(aTS, aTP) = Rubidium(true) # Stuttgart data, "true" to include P-wave scattering

@time begin
println("prepare calculation...")
# Define parameters ____________________________________________________________
  n = 30 # energy level
  m_max = 2 # Rydberg electron total angular momentum projection
  Δnp = 0 # additional above manifolds
  Δnm = 1 # additional below manifolds
  (δRb, α) = QDaP(false) # "true" for Cesium, "false" for Rubidium
  ρ = linspace(0.4, 2.4, 400)*n^2; # radial grid
  writecsv("$Path/pec/data/pec_n=$(n)_grid.csv", ρ)
  F = 1000 # field strength

# Define basis states __________________________________________________________
  b_hyd = Hydrogenic(n, Δnp, Δnm, m_max) # hydrogenic basis states
  b_qd = QDefect(n, Δnp, Δnm, m_max, δRb) # quantum defect basis states
  basis = vcat(b_qd, b_hyd)
  Vbasis = basis[:,1] - δnlj(basis[:,2], δRb) # effective primary quantum number

  nl = unique(basis[:,1:2], 1) # unperturbed basis states
  lm = unique([basis[:,2] abs.(basis[:,3])], 1)
  m_l = sort(unique(basis[:,3]))

# Set up constant matrices _____________________________________________________
  V0 = spdiagm(-1./2Vbasis.^2)
  (Rbasis, Pbasis, Abasis) = Pos(basis, nl, lm, m_l)
end

@time begin
println("calculate...")
for ω in linspace(0.,π,3)
  println("θ=$(ω/π)")
  H_E = StarkTerm(basis, ω, F, n)
  # @time PES, CLP = Calc_ves(ρ, n, nl, lm, m_l, basis, Rbasis, Pbasis, Abasis, δRb, α, V0, aTS, aTP, H_E)
  @time PES = Calc_pes(ρ, n, nl, lm, m_l, basis, Rbasis, Pbasis, Abasis, δRb, α, V0, aTS, aTP, H_E)
  # EE = (PES + 1/2n^2)*2*3.289841960355*10^6 # in GHz
  EE = (PES + 1/2n^2)
  writecsv("$Path/pec/data/pec_n=$(n)_E=$(F)_θ=$(ω/π).csv", EE)
end
end

