DataPath = "/home/frederic/Documents/julia/data";

# Phase shift data _____________________________________________________________
function scatdata(Cs, swaveS, swaveT, pwaveS, pwaveT, pwaveSO, cfit) # import scattering length
    if Cs
        aSP = readcsv("$DataPath/CsPhases_3.csv")

        a0 = extrapolate(interpolate((aSP[:,1], ), zeros(length(aSP[:,1])), Gridded(Linear())), Flat())
        aS0 = extrapolate(interpolate((aSP[:,1], ), aSP[:,2], Gridded(Linear())), Flat())
        aS1 = extrapolate(interpolate((aSP[:,1], ), aSP[:,3], Gridded(Linear())), Flat())
        aP0 = extrapolate(interpolate((aSP[:,1], ), aSP[:,4], Gridded(Linear())), Flat())
        aP10 = extrapolate(interpolate((aSP[:,1], ), aSP[:,5], Gridded(Linear())), Flat())
        aP11 = extrapolate(interpolate((aSP[:,1], ), aSP[:,6], Gridded(Linear())), Flat())
        aP12 = extrapolate(interpolate((aSP[:,1], ), aSP[:,7], Gridded(Linear())), Flat())
        aP1 = aP11
    else
        aSPm = readcsv("$DataPath/Rbscatl_matt.csv")
        if cfit == "Matt"
            a0 = extrapolate(interpolate((aSPm[:,1], ), zeros(length(aSPm[:,1])), Gridded(Linear())), Flat())
            aS0 = extrapolate(interpolate((aSPm[:,1], ), aSPm[:,2], Gridded(Linear())), Flat())
            aS1 = extrapolate(interpolate((aSPm[:,1], ), aSPm[:,3], Gridded(Linear())), Flat())
            aP0 = extrapolate(interpolate((aSPm[:,1], ), aSPm[:,4], Gridded(Linear())), Flat())
            aP10 = extrapolate(interpolate((aSPm[:,1], ), aSPm[:,5], Gridded(Linear())), Flat())
            aP11 = extrapolate(interpolate((aSPm[:,1], ), aSPm[:,6], Gridded(Linear())), Flat())
            aP12 = extrapolate(interpolate((aSPm[:,1], ), aSPm[:,7], Gridded(Linear())), Flat())
            aP1 = aP11
        else
            aSP = readcsv("$DataPath/Rbscatl_fit_c=$(cfit).csv")

            cut = 201
            aSP[1:cut,3:end] = repeat(aSP[cut,3:end]', inner=(cut,1))

            a0 = extrapolate(interpolate((aSP[:,1], ), zeros(length(aSP[:,1])), Gridded(Linear())), Flat())
            aS0 = extrapolate(interpolate((aSPm[:,1], ), aSPm[:,2], Gridded(Linear())), Flat())
            aS1 = extrapolate(interpolate((aSP[:,1], ), aSP[:,2], Gridded(Linear())), Flat())
            aP0 = extrapolate(interpolate((aSPm[:,1], ), aSPm[:,4], Gridded(Linear())), Flat())
            aP1 = extrapolate(interpolate((aSP[:,1], ), aSP[:,3], Gridded(Linear())), Flat())
            aP10 = extrapolate(interpolate((aSP[:,1], ), aSP[:,4], Gridded(Linear())), Flat())
            aP11 = extrapolate(interpolate((aSP[:,1], ), aSP[:,5], Gridded(Linear())), Flat())
            aP12 = extrapolate(interpolate((aSP[:,1], ), aSP[:,6], Gridded(Linear())), Flat())
        end
    end

    if swaveS == false
        aS0 = a0
    end
    if swaveT == false
        aS1 = a0
    end
    if pwaveS == false
        aP0 = a0
    end
    if pwaveSO == false
        aP10 = aP1
        aP11 = aP1
        aP12 = aP1
    end
    if pwaveT == false
        aP10 = a0
        aP11 = a0
        aP12 = a0
    end
    return (aS0, aS1, aP0, aP10, aP11, aP12)
end

function scatdata(λ1, λ2, λ3) # import scattering length
    aSPm = readcsv("$DataPath/Rbscatl_matt.csv")
    aSP = readcsv("$DataPath/Rbscatl_fit_c=$(cfit).csv")[1:size(aSPm,1),:]

    cut = 201
    aSP[1:cut,3:end] = repeat(aSP[cut,3:end]', inner=(cut,1))

    aS0 = extrapolate(interpolate((aSPm[:,1], ), (1-λ1)*aSP[:,2]+λ1*aSPm[:,2], Gridded(Linear())), Flat())
    aS1 = extrapolate(interpolate((aSP[:,1], ), aSP[:,2], Gridded(Linear())), Flat())
    aP0 = extrapolate(interpolate((aSPm[:,1], ), (1-λ2)*aSP[:,3]+λ2*aSPm[:,4], Gridded(Linear())), Flat())
    aP10 = extrapolate(interpolate((aSP[:,1], ), (1-λ3)*aSP[:,3]+λ3*aSP[:,4], Gridded(Linear())), Flat())
    aP11 = extrapolate(interpolate((aSP[:,1], ), (1-λ3)*aSP[:,3]+λ3*aSP[:,5], Gridded(Linear())), Flat())
    aP12 = extrapolate(interpolate((aSP[:,1], ), (1-λ3)*aSP[:,3]+λ3*aSP[:,6], Gridded(Linear())), Flat())

    return (aS0, aS1, aP0, aP10, aP11, aP12)
end

E2n(energy) = 1./sqrt.(2energy)
n2E(nstar) = 1./2nstar.^2

kinetic(R, n) = real(sqrt(complex(2(-1/2n^2 + 1/R)))) # kinetic energy

function a(L, S, J, R, n, aS0, aS1, aP0, aP10, aP11, aP12)
  if (S, L, J) == (0, 0, 0)
    aS0[kinetic(R, n)]
  elseif (S, L, J) == (1, 0, 1)
    aS1[kinetic(R, n)]
  elseif (S, L, J) == (0, 1, 1)
    aP0[kinetic(R, n)]
  elseif (S, L, J) == (1, 1, 0)
    aP10[kinetic(R, n)]
  elseif (S, L, J) == (1, 1, 1)
    aP11[kinetic(R, n)]
  elseif (S, L, J) == (1, 1, 2)
    aP12[kinetic(R, n)]
  else
    println("impossible combination of (L,S,J)")
  end # scattering length # scattering channel selection
end

# Define basis set _____________________________________________________________
function Species(Cs)
    if Cs
        I1 = 3.5 # nuclear spin
        α = 402.2 # polarisability
        α2 = 15.8 # ion polarizability
    else
        I1 = 1.5 # nuclear spin
        α = 319.2 # polarisability
        α2 = 9.11 # ion polarizability
    end
    return (I1, α, α2)
end

function Fine(Cs, fine, n) # quantum defects
    if Cs
        if fine
            δRb1 = [4.049325 4.049325; 3.591556 3.559058; 2.475365 2.466210; 0.033392 0.033537]
            δRb2 = [0.2462 0.2462; 0.3714 0.374; 0.5554 0.067; -0.191 0.191]
            return δRb1 + δRb2./(n - δRb1).^2
        else
            δRb1 = [4.049325 4.049325; 3.57 3.57; 2.477 2.47; 0.0335 0.335]
            δRb2 = [0.2462 0.2462; 0.374 0.374; 0.067 0.067; -0.191 0.191]
            return δRb1 + δRb2./(n - δRb1).^2
        end
    else
        if fine
            δRb1 = [3.1311804 3.1311804; 2.6548849 2.6416737; 1.34809171 1.34646572; 0.0165192 0.0165437]
            δRb2 = [0.1784 0.1784; 0.29 0.295; -0.60286 -0.596; -0.085 -0.086]
            return δRb1 + δRb2./(n - δRb1).^2
        else
            δRb1 = [3.1311804 3.1311804; 2.65 2.65; 1.35 1.35; 0.016 0.016]
            δRb2 = [0.1784 0.1784; 0.29 0.29; -0.6 -0.6; -0.085 -0.085]
            return δRb1 + δRb2./(n - δRb1).^2
        end
    end
end

function Nstar(basis, δRb, α2, fine, Sansonetti)
    Vbasis = basis[:,1] - δnlj(round.(Int, basis[:,2]), basis[:,3], δRb) - QDCorePol(basis, α2)
    if Sansonetti
        if fine
            δRb_exp = readcsv("$DataPath/Rb_quantum_defects.csv")
            for i in eachindex(basis[:,1]), j in eachindex(δRb_exp[:,1])
                if basis[i,1:3] == δRb_exp[j,1:3]
                    Vbasis[i] = δRb_exp[j,4]
                end
            end
        end
    end
    return Vbasis
end

function Hyperfine(Cs, hyperfine) # hyperfine splitting
  if hyperfine
      if Cs
          Ahf = 2.298/(2*3.289841960355*10^6)
      else
          Ahf = 3.417/(2*3.289841960355*10^6)
      end
  else
    Ahf = 0.
  end
end

function FineStructure(basis)
    ΔE = Float64[]
    for (i,l) in enumerate(basis[:,2])
        if l < 4
            push!(ΔE, 0.)
        else
            n = basis[i,1]
            j = basis[i,3]
            len = length(n)
            push!(ΔE, (1 ./(j+0.5) - 3./4n) ./2n.^3)
        end
    end
    return -ΔE/137^2
end

function QDCorePol(basis, α2)
    μ = Float64[]
    for (i,l) in enumerate(basis[:,2])
        if l < 4
            push!(μ, 0.)
        else
            n = basis[i,1]
            j = basis[i,3]
            push!(μ, (3n.^2 - l.*(l+1))./(4n.^2 .*(l-0.5).*l.*(l+0.5).*(l+1).*(l+1.5)))
        end
    end
    return α2*μ
end

function δnlj(l::Int, j::Float64, δRb) # quantum defect
  if l<size(δRb,1)
    if j>l
      δRb[l+1,2]
    else
      δRb[l+1,1]
    end
  else
    0.
  end
end

function δnlj(l::Array, j::Array, δRb)
  b = Float64[]
  for i in eachindex(l)
    if l[i] < size(δRb,1)
      if j[i] > l[i]
        push!(b, δRb[l[i]+1,2])
      else
        push!(b, δRb[l[i]+1,1])
      end
    else
      push!(b, 0.)
    end
  end
  return b
end

function pos(n, l, j , list)
  pos = Int[]
  for i in eachindex(list[:,1])
    if list[i,1:3] == [n,l,j]
      push!(pos, i)
    end
  end
  return pos
end

function Pos(nlj, basis)
  b = Int[]
  for i in eachindex(basis[:,1])
    append!(b, pos(basis[i,1], basis[i,2], basis[i,3], nlj))
  end
  return b
end

function Hydrogenic(n, n_below, n_above, m_max, I0) # hydrogenic states with hyperfine interaction
  b = Float64[]
  for i in n-n_below:n+n_above
    for j in 4:i-1
      for k in j-0.5:j+0.5
        for l in -min(k, m_max):min(k, m_max)
          for p in -0.5:0.5
            for q in -I0:I0
              push!(b, i, j, k, l, p, q)
            end
          end
        end
      end
    end
  end
  return reshape(b, 6, round(Int, length(b)/6))'
end

function Hydrogenic(n, L, m_max, I0) # hydrogenic states with hyperfine interaction
    b = Float64[]
    for j in 3:n-1
      for k in j-0.5:j+0.5
        for l in -min(k, m_max):min(k, m_max)
          for p in -0.5:0.5
            for q in -I0:I0
              push!(b, n, j, k, l, p, q)
            end
          end
        end
      end
    end
    return reshape(b, 6, round(Int, length(b)/6))'
end

function QDefect(n, n_below, n_above, m_max, I0, δRb) # quantum defect states with hyperfine interaction
  b = Float64[]
  for i in n-n_below:n+n_above
    for j in 0:3
      for k in j-0.5:j+0.5
        for l in -min(k, m_max):min(k, m_max)
          for p in -0.5:0.5
            for q in -I0:I0
              push!(b, i, j, k, l, p, q)
            end
          end
        end
      end
    end
  end
  b = reshape(b, 6, round(Int, length(b)/6))'
  for i in eachindex(b[:,1])
    b[i,1] += round(δnlj(round(Int, b[i,2]), b[i,3], δRb), RoundDown)
  end
  return b
end

function QDefect(n, L, m_max, I0, δRb) # quantum defect states with hyperfine interaction
    b = Float64[]
    for j in L
      for k in j-0.5:j+0.5
        for l in -min(k, m_max):min(k, m_max)
          for p in -0.5:0.5
            for q in -I0:I0
              push!(b, n, j, k, l, p, q)
            end
          end
        end
      end
    end
    b = reshape(b, 6, round(Int, length(b)/6))'
    for i in eachindex(b[:,1])
      b[i,1] += round(δnlj(round(Int, b[i,2]), b[i,3], δRb), RoundDown)
    end
    return b
end

function lowStates(maxN, maxM, I0)
    b = Float64[]
    c = Float64[]
    nljδ = readcsv("$DataPath/Rb_quantum_defects.csv")
    for i in 1:size(nljδ,1)
        if nljδ[i,4] <= maxN
            for mj in -min(nljδ[i,3], maxM):min(nljδ[i,3], maxM)
                for ms in -0.5:0.5
                    for mF in -I0:I0
                        push!(b, nljδ[i,1], nljδ[i,2], nljδ[i,3], mj, ms, mF)
                        push!(c, nljδ[i,4])
                    end
                end
            end
        end
    end
    for n in 6:maxN
        if n < 9
            for l in 4:n-1
                for j in l-0.5:l+0.5
                    for mj in -min(j, maxM):min(j, maxM)
                        for ms in -0.5:0.5
                            for mF in -I0:I0
                                push!(b, n, l, j, mj, ms, mF)
                                push!(c, n)
                            end
                        end
                    end
                end
            end
        elseif n == 9
            for l in 3:n-1
                for j in l-0.5:l+0.5
                    for mj in -min(j, maxM):min(j, maxM)
                        for ms in -0.5:0.5
                            for mF in -I0:I0
                                push!(b, n, l, j, mj, ms, mF)
                                push!(c, n)
                            end
                        end
                    end
                end
            end
        else
            println("Error: Upward basis construction mathod only valid up to n=9 hydrogenic states.")
            println("n=$n selected. Please select n=9 or lower.")
        end
    end
    return (reshape(b, 6, round(Int, length(b)/6))', c)
end

function Vselect(basis, Vbasis, Ω) # select basis set with certain Ω
  b = Float64[]
  for i in eachindex(Vbasis)
    if sum(basis[i,4:6]) .== Ω
        append!(b, Vbasis[i])
    end
  end
  return b
end

function bselect(basis, Ω) # select basis set with certain Ω
  b = Float64[]
  for i in eachindex(basis[:,1])
    if sum(basis[i,4:6]) .== Ω
      append!(b, basis[i,:])
    end
  end
  return reshape(b, 6, round(Int, length(b)/6))'
end

function initializebasis(n, L, m_j_max, I1, δRb::Array{Float64,2}, hyd::Bool)
    if hyd
        b_hyd = Hydrogenic(n, L, m_j_max, I1) # hydrogenic basis states
        b_qd = QDefect(n, L, m_j_max, I1, δRb) # quantum defect basis states
        basis = vcat(b_qd, b_hyd)
    else
        basis = QDefect(n, L, m_j_max, I1, δRb)
    end
    return basis
end

function initializebasis(n, n_below, n_above, m_j_max, I1, δRb::Array{Float64,2})
    b_hyd = Hydrogenic(n, n_below, n_above, m_j_max, I1) # hydrogenic basis states
    b_qd = QDefect(n, n_below, n_above, m_j_max, I1, δRb) # quantum defect basis states
    basis = vcat(b_qd, b_hyd)
    return basis
end

function LSJΩ() # list of scattering relevant spin states
  b = Float64[]
  for l in 0:1
    for s in 0:1
      for j in abs(l-s):abs(l+s)
        push!(b, l, s, j)
      end
    end
  end
  b = reshape(b, 3, round(Int, length(b)/3))'
end

function LSJΩ(Ω) # list of scattering relevant spin states
  b = Float64[]
  for l in 0:1
    for s in 0:1
      for j in abs(l-s):abs(l+s)
        push!(b, l, s, j, Ω)
      end
    end
  end
  b = reshape(b, 4, round(Int, length(b)/4))'
end

function Ψ_ρ(ρ, n, l::Int)
  # in terms of Whittaker functions for the quantum defect states
  f(N, L) = 1/ρ* 1/N* 1/sqrt(sf_gamma(N+L+1)* sf_gamma(N-L))* exp(-ρ/N)* (2ρ/N)^(L+1)* sf_hyperg_U(L+1-N, 2(L+1), 2ρ/N)
  # and in terms of hydrogenic functions for regular states
  g(N, L) = 2/N^2* sqrt(sf_gamma(N-L)/sf_gamma(N+L+1))* exp(-ρ/N)* (2ρ/N)^L* sf_laguerre_n(N-L-1, 2L+1, 2ρ/N)
  if l < 4
    f(n, l)
  else
    g(round(Int, n), round(Int, l))
  end # radial basis function # radial basis function
end

Ψ_d(r) = Ψ_ρ(r, nlist[tmp], 0)

function Ψ_ρ(ρ, n::Array, l::Array)
  # in terms of Whittaker functions for the quantum defect states
  f(N, L) = 1/ρ* 1/N* 1/sqrt(sf_gamma(N+L+1)* sf_gamma(N-L))* exp(-ρ/N)* (2ρ/N)^(L+1)* sf_hyperg_U(L+1-N, 2(L+1), 2ρ/N)
  # and in terms of hydrogenic functions for regular states
  g(N, L) = 2/N^2* sqrt(sf_gamma(N-L)/sf_gamma(N+L+1))* exp(-ρ/N)* (2ρ/N)^L* sf_laguerre_n(N-L-1, 2L+1, 2ρ/N)
  b = Float64[]
  for (i, l0) in enumerate(l)
    if l0 < 4
      push!(b, f(n[i], l0))
    else
      push!(b, g(round(Int, n[i]), round(Int, l0)))
    end
  end
  return b
end

function Ψ_θ(θ, l::Array, j::Array, mj::Array) # azimutal basis function
    function g(l, m)
        C = 0.
        if abs(m) <= l
            if m >= 0
                C = sf_legendre_sphPlm(l, m, cos(θ))
            else
                C = (-1)^(-m)* sf_legendre_sphPlm(l, -m, cos(θ))
            end
        end
        return C
    end
    function f(l, j, mj)
        C = 0.
        for ms in -0.5:0.5
            C += cg(l, j, mj, ms).*g(round(Int,l), round(Int,mj-ms))
        end
        return C
    end
    return map(f, l, j, mj)
end

function Q(n, l, j, L, M, R, δRb) # Matt equation (10-12)
  Ψ(ρ) = Ψ_ρ(ρ, n-δnlj(round.(Int, l), j, δRb), l)
  if (L, M) == (0, 0)
    sqrt.((2l+1)/4π).*Ψ(R) # Matt equation (10)
  elseif (L, M) == (1, 0)
    sqrt.((2l+1)/4π).*derivative(Ψ, R) # Matt equation (11)
  elseif (L, abs(M)) == (1, 1)
    sqrt.((2l+1).*l.*(l+1)/8π).*Ψ(R)/R # Matt equation (12)
  else
    0. *l
  end
end

# Set up constant matrices _____________________________________________________
function W3j(j1, m1, j2, m2, J, M) # Wigner 3-j symbols for calculation of Clebsch-Gordon coefficients
  (-1)^(j1-j2+M)* sqrt(2J+1)* sf_coupling_3j(round(Int,2j1), round(Int,2j2), round(Int,2J), round(Int,2m1), round(Int,2m2), -round(Int,2M))
end

cg(l, j, mj, ms) = W3j(l, mj-ms, 0.5, ms, j, mj) # regular Clebsch-Gordon coefficient for s=0.5, ml = mj - ms

function cg(Ml, L, S, J, Ω, n, l, j, m, mS, mI) # Matt equation (19) R-independent
  if Ω == m+mS+mI && abs(Ml) <= l && abs(m-Ml) <= 0.5 && abs(m+mS-Ml) <= S && abs(mS+m) <= J && abs(Ml) <= L
    sqrt(4π/(2L+1))* W3j(l,Ml, 0.5,m-Ml, j,m)* W3j(0.5,m-Ml, 0.5,mS, S,m+mS-Ml)* W3j(L,Ml, S,m+mS-Ml, J,mS+m)
  else
    0.
  end
end

function CG(Ml, LSJ, basis) # Matt equation (19) R-independent
  b = Array{Float64}(size(basis,1), size(LSJ,1))
  for i in eachindex(basis[:,1])
    for j in eachindex(LSJ[:,1])
      b[i,j] = cg(Ml, LSJ[j,:]'..., basis[i,:]'...)
    end
  end
  return b
end

function CG(Ml, LSJ, Ω, basis) # Matt equation (19) R-independent
  b = Array{Float64}(size(basis,1), size(LSJ,1))
  for i in eachindex(basis[:,1])
    for j in eachindex(LSJ[:,1])
      b[i,j] = cg(Ml, LSJ[j,:]'..., Ω[i], basis[i,:]'...)
    end
  end
  return b
end

eF(F, I1, Ahf) = Ahf/2* (F*(F+1)-I1*(I1+1)-0.5(0.5+1))

function hf(mS1, mI1, mS2, mI2, Ahf, Fs, I1) # hyperfine interaction
  b = Float64[]
  for F in Fs
    for mF in -F:F
      if mS1+mI1 == mF && mS2+mI2 == mF
        push!(b, (F*(F+1)-I1*(I1+1)-0.5(0.5+1))* W3j(0.5,mS1, I1,mI1, F,mF)* W3j(0.5,mS2, I1,mI2, F,mF))
      else
        push!(b, 0.)
      end
    end
  end
  return sum(b)*Ahf/2
end

function HF(basis, Ahf, Fs, I1)
  b = zeros(size(basis,1), size(basis,1))
  tic()
  for i in eachindex(basis[:,1])
    for j in eachindex(basis[:,1])
      if basis[i,1:4] == basis[j,1:4]
        b[i,j] = hf(basis[i,5:6]'..., basis[j,5:6]'..., Ahf, Fs, I1)
      end
    end
  end
  toc()
  return sparse(b)
end

function Hfilt(basis, Ω) # hyperfine coupling is diagonal
  b = zeros(Int, size(basis,1), size(basis,1))
  for i in eachindex(basis[:,1])
    for j in eachindex(basis[:,1])
      if basis[i,6] == basis[j,6]
        if Ω[i] == Ω[j]
          b[i,j] = 1
        end
      end
    end
  end
  return sparse(b)
end

function Hfilt(basis) # hyperfine coupling is diagonal
  b = zeros(Int, size(basis,1), size(basis,1))
  for i in eachindex(basis[:,1])
    for j in eachindex(basis[:,1])
      if basis[i,6] == basis[j,6]
        b[i,j] = 1
      end
    end
  end
  return sparse(b)
end

function S2filt(basis) # magnetic field couples only certain qunatum states
  b = zeros(Float64, size(basis,1), size(basis,1))
  for i in eachindex(basis[:,1]), j in eachindex(basis[:,1])
    if basis[i,1:4] == basis[j,1:4] && basis[i,6] == basis[j,6]
      b[i,j] = 1.
    end
  end
  return sparse(b)
end

function S1filt(basis) # magnetic field couples only certain qunatum states
  z = zeros(Int, size(basis,1), size(basis,1))
  p = zeros(Int, size(basis,1), size(basis,1))
  for i in eachindex(basis[:,1]), j in eachindex(basis[:,1])
    if basis[i,1:2] == basis[j,1:2] && basis[i,5:6] == basis[j,5:6]
      if basis[i,4] == basis[j,4]
        z[i,j] = 1
      end
      if basis[i,4] == basis[j,4]+1
        p[i,j] = 1
      end
    end
  end
  return (sparse(z), sparse(p))
end

function Sin(ω)
  if ω == π || ω == -π
    0.
  else
    sin.(ω)
  end
end

function Cos(ω)
  if ω == π/2 || ω == -π/2
    0.
  else
    cos.(ω)
  end
end

function magnetic(basis, θ, B, I1)
  if B != 0
    # J component
    J = basis[:,3]
    MJ = basis[:,4]
    J_p = Sin(θ)*0.5sqrt.(J.*(J+1)-MJ.*(MJ+1))
    J_m = Sin(θ)*0.5sqrt.(J.*(J+1)-MJ.*(MJ-1))
    H_J = Cos(θ)*spdiagm(MJ)

    # S1 component
    C1(l,j,mj) = cg(l,j,mj,-0.5)
    C2(l,j,mj) = cg(l,j,mj,0.5)
    cminus = map(C1, basis[:,2], basis[:,3], basis[:,4])
    cplus = map(C2, basis[:,2], basis[:,3], basis[:,4])
    (z, p) = S1filt(basis)
    S1p = cplus*cminus'.*p
    H_S1 = Cos(θ)/2*(cplus*cplus' - cminus*cminus').*z + Sin(θ)/2*(S1p + S1p')
    # H_S1 = Cos(θ)*diagm((-1).^(L+0.5-J).*MJ./(1+2L)) + Sin(θ)/2*(S1p + S1p')

    # S2 component
    ms = basis[:,5]
    H_S2 = Cos(θ)*spdiagm(ms) - Sin(θ)*(ms*ms'-1/4).*S2filt(basis)

    # I2 component
    mi = basis[:,6]
    I_p = Sin(θ)*0.5sqrt.(I1*(I1+1)-mi.*(mi+1))
    I_m = Sin(θ)*0.5sqrt.(I1*(I1+1)-mi.*(mi-1))
    H_I = Cos(θ)*spdiagm(mi)

    for i in eachindex(I_p), j in eachindex(I_p)
      if basis[i,1:3] == basis[j,1:3]
        if basis[i,5:6] == basis[j,5:6]
          if MJ[i]+1 == MJ[j]
            H_J[i,j] = J_p[i]
          elseif MJ[i]-1 == MJ[j]
            H_J[i,j] = J_m[i]
          end
        elseif basis[i,4:5] == basis[j,4:5]
          if mi[i]+1 == mi[j]
            H_I[i,j] = I_p[i]
          elseif mi[i]-1 == mi[j]
            H_I[i,j] = I_m[i]
          end
        end
      end
    end

    H_B = B/2.35051746439E9*(H_J + H_S1 + 2H_S2)/2
    return H_B, H_J, H_S2 + H_I, H_S1 + H_S2 + H_I
  else
    H_B = spzeros(size(basis,1), size(basis,1))
    return H_B, H_B, H_B, H_B
  end
end

# Define diagonalisation procedure _____________________________________________
function qred(Ml, R, nlj, Ls, δRb)
  b = zeros(size(nlj,1), size(Ls,1))
  for (i,L) in enumerate(Ls)
    if abs(Ml) <= L
      b[:,i] = Q(nlj[:,1], nlj[:,2], nlj[:,3], L, Ml, R, δRb)
    end
  end
  return b
end

function qmat(Ml, R, nlj, Ls, δRb, Pnlj, Pl)
  q = qred(Ml, R, nlj, Ls, δRb)
  b = Float64[]
  for i in Pnlj
    for j in Pl
      push!(b, q[i,j])
    end
  end
  return reshape(b, length(Pl), length(Pnlj))'
end

function A(R, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm) # Matt equation (19)
  q0 = qmat(0, R, nlj, Ls, δRb, Pnlj, Pl)
  qp = qmat(1, R, nlj, Ls, δRb, Pnlj, Pl)
  qm = qmat(-1, R, nlj, Ls, δRb, Pnlj, Pl)

  return CG0.*q0 + CGp.*qp + CGm.*qm
end

function U(R, n, LSJ, aS0, aS1, aP0, aP10, aP11, aP12) # Matt equation (18)
  b = Float64[]
  for i in eachindex(LSJ[:,1])
    push!(b, (2LSJ[i,1]+1)^2/2* a(LSJ[i,1:3]'..., R, n, aS0, aS1, aP0, aP10, aP11, aP12))
  end
  return diagm(b)
end

function U(R, n, J, aP10, aP11, aP12) # Matt equation (18)
    b = zeros(6)
    if J == 0
        b[4] = 9/2*aP10[kinetic(R, n)]
    elseif J == 1
        b[5] = 9/2*aP11[kinetic(R, n)]
    elseif J == 2
        b[6] = 9/2*aP12[kinetic(R, n)]
    else
        nothing
    end
    return diagm(b)
end

function U(J) # for J_Projector
    b = zeros(6)
    if J == 0
        b[4] = 1
    elseif J == 1
        b[5] = 1
    elseif J == 2
        b[6] = 1
    else
        nothing
    end
    return diagm(b)
end

Pol(R, α) = -α./2R.^4 # polarisation potential

function H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12) # Matt equation (20)
  a = A(R, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
  u = U(R, n, LSJ, aS0, aS1, aP0, aP10, aP11, aP12)
  return V0 + (a*u*a').*Delta
end

function J_Projector(J, R, Delta, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm) # Matt equation (20)
  a = A(R, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
  u0 = (a*U(0)*a').*Delta
  u1 = (a*U(1)*a').*Delta
  u2 = (a*U(2)*a').*Delta
  u = u0 + u1 + u2
  if J == 0
      return u0/norm(Array(u))
  elseif J == 1
      return u1/norm(Array(u))
  elseif J == 2
      return u2/norm(Array(u))
  else
      return u
  end
end

function J_Contribution(J, R, V0, Delta, n, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aP10, aP11, aP12) # Matt equation (20)
    a = A(R, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
    u = U(R, n, J, aP10, aP11, aP12)
    return V0 + (a*u*a').*Delta
end

# Define diagonalisation procedure _____________________________________________
function Calc_pes(ρ, V0, Delta, H_B, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α)
  pes = Array{Float64}(size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H) and Diagonalisation (D)):")
  for (i,R) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed pes[:,i] = real(eigvals(full(hamiltonian) + H_B)) + P(R, α)
    println("$i of $(length(ρ)): H: $(round(Htime,5))s; D: $(round(Dtime,5))s.")
  end
  return pes'
end

function Calc_ang_pes_FULL(θ, B, ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α, I1)
    pes = Array{Float64}(size(basis,1), length(θ))
    println("Starting calculation for $(size(basis,1)) basis states at $(length(θ)) grid points.")
    println("Timing of the Hamiltonian setup (H) and Diagonalisation (D):")
    for (i,t) in enumerate(θ)
      H_B, P_MJ, P_MF, P_Ω = magnetic(basis, t, B, I1)
      Htime = @elapsed hamiltonian = H(ρ, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
      Dtime = @elapsed  begin
          EigenSystem = eigfact(full(hamiltonian) + H_B)
          pes[:,i] = real(EigenSystem[:values]) + Pol(ρ, α)
      end
      println("$i of $(length(θ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s.")
    end
    return pes'
end

function Calc_ang_pes_ARNOLDI(σ, Num, θ, B, ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α, I1)
  pes = Array{Float64}(Num, length(θ))
  v_pes = zeros((0, ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(θ)) grid points.")
  println("Timing of the Hamiltonian setup (H) and Diagonalisation (D):")
  for (i,t) in enumerate(θ)
    H_B, P_MJ, P_MF, P_Ω = magnetic(basis, t, B, I1)
    Htime = @elapsed hamiltonian = H(ρ, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed  begin
        EigenSystem = eigs(full(hamiltonian) + H_B, nev=Num, sigma=σ, ritzvec=false, v0=v_pes)
        v_pes = vcat(real(EigenSystem[1]), zeros(size(basis,1)-Num))
        pes[:,i] = real(EigenSystem[1]) + Pol(ρ, α)
    end
    println("$i of $(length(θ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s.")
  end
  return pes'
end

function Calc_rad_pes(θ, B, ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α, I1)
  pes = Array{Float64}(size(basis,1), length(ρ))
  c1 = Array{Float64}(size(basis,1), length(ρ))
  c2 = Array{Float64}(size(basis,1), length(ρ))
  H_B, P_MJ, P_MF, P_Ω = magnetic(basis, θ, B, I1)
  V_mF = eigfact(FmF(basis, θ, I1), 1.99, 2.01)[:vectors]
  V_mJu = eigfact(JmJ(basis, θ), 0.49, 0.51)[:vectors]
  V_mJd = eigfact(JmJ(basis, θ), -0.51, -0.49)[:vectors]
  P_mF = V_mF*V_mF'
  P_mJu = V_mJu*V_mJu'
  P_mJd = V_mJd*V_mJd'
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and chosen Projection (P):")
  for (i,r) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(r, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed  begin
        EigenSystem = eigfact(full(hamiltonian) + H_B)
        pes[:,i] = real(EigenSystem[:values]) + Pol(r, α)
    end
    Ctime = @elapsed begin
        v = real(EigenSystem[:vectors])
        c1[:,i] = diag(v'*(P_mF*P_mJu)*v)
        c2[:,i] = diag(v'*(P_mF*P_mJd)*v)
    end
    println("$i of $(length(ρ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ctime,5))s.")
  end
  return pes', c1', c2'
end

function Calc_ang_pes(θ, B, ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α, I1)
  pes = Array{Float64}(size(basis,1), length(θ))
  c1 = Array{Float64}(size(basis,1), length(θ))
  c2 = Array{Float64}(size(basis,1), length(θ))
  P_l = spdiagm(ifelse.(basis[:,2].==0, 1, 0))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(θ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and chosen Projection (P):")
  for (i,t) in enumerate(θ)
    H_B, P_MJ, P_MF, P_Ω = magnetic(basis, t, B, I1)
    Htime = @elapsed hamiltonian = H(ρ, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed  begin
        EigenSystem = eigfact(full(hamiltonian) + H_B)
        pes[:,i] = real(EigenSystem[:values]) + Pol(ρ, α)
    end
    Ctime = @elapsed begin
        v = real(EigenSystem[:vectors])
        V_mF = eigfact(FmF(basis, t, I1), 1.99, 2.01)[:vectors]
        V_mJu = eigfact(JmJ(basis, t), 0.49, 0.51)[:vectors]
        V_mJd = eigfact(JmJ(basis, t), -0.51, -0.49)[:vectors]
        P_mF = V_mF*V_mF'
        P_mJu = V_mJu*V_mJu'
        P_mJd = V_mJd*V_mJd'
        c1[:,i] = diag(v'*(P_l*P_mF*P_mJu)*v)
        c2[:,i] = diag(v'*(P_l*P_mF*P_mJd)*v)
    end
    println("$i of $(length(θ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ctime,5))s.")
  end
  return pes', c1', c2'
end

function Calc_ang_pes_pj(θ, B, ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α, I1)
  pes = Array{Float64}(size(basis,1), length(θ))
  c0 = Array{Float64}(size(basis,1), length(θ))
  c1 = Array{Float64}(size(basis,1), length(θ))
  c2 = Array{Float64}(size(basis,1), length(θ))
  jprojector0 = J_Projector(0, ρ, Delta, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
  jprojector1 = J_Projector(1, ρ, Delta, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
  jprojector2 = J_Projector(2, ρ, Delta, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
  println("Starting calculation for $(size(basis,1)) basis states at $(length(θ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and J Projections (P):")
  for (i,t) in enumerate(θ)
    H_B, P_MJ, P_MF, P_Ω = magnetic(basis, t, B, I1)
    Htime = @elapsed hamiltonian = H(ρ, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed begin
        EigenSystem = eigfact(full(hamiltonian) + H_B)
        pes[:,i] = real(EigenSystem[:values]) + Pol(ρ, α)
    end
    Ptime = @elapsed begin
        v = real(EigenSystem[:vectors])
        c0[:,i] = diag(v'*jprojector0*v)
        c1[:,i] = diag(v'*jprojector1*v)
        c2[:,i] = diag(v'*jprojector2*v)
    end
    println("$i of $(length(θ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ptime,5))s.")
  end
  return pes', c0', c1', c2'
end

function Calc_ang_pes_Jcontr(θ, B, ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α, I1)
  pes = Array{Float64}(size(basis,1), length(θ))
  c0 = Array{Float64}(size(basis,1), length(θ))
  c1 = Array{Float64}(size(basis,1), length(θ))
  c2 = Array{Float64}(size(basis,1), length(θ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(θ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and J Contributions (C):")
  for (i,t) in enumerate(θ)
    H_B, P_MJ, P_MF, P_Ω = magnetic(basis, t, B, I1)
    Htime = @elapsed begin
        hamiltonian = H(ρ, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
        J0 = J_Contribution(0, ρ, V0, Delta, n, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aP10, aP11, aP12)
        J1 = J_Contribution(1, ρ, V0, Delta, n, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aP10, aP11, aP12)
        J2 = J_Contribution(2, ρ, V0, Delta, n, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aP10, aP11, aP12)
    end
    Dtime = @elapsed begin
        EigenSystem = eigfact(full(hamiltonian) + H_B)
        pes[:,i] = real(EigenSystem[:values]) + Pol(ρ, α)
    end
    Ctime = @elapsed begin
        c0[:,i] = real(eigfact(full(J0) + H_B)[:values]) + Pol(ρ, α)
        c1[:,i] = real(eigfact(full(J1) + H_B)[:values]) + Pol(ρ, α)
        c2[:,i] = real(eigfact(full(J2) + H_B)[:values]) + Pol(ρ, α)
    end
    println("$i of $(length(θ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ctime,5))s.")
  end
  return pes', c0', c1', c2'
end

# Define diagonalisation procedure without Bfield ______________________________
function Calc_pes(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α)
  pes = Array{Float64}(size(basis)[1], length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H) and Diagonalisation (D)):")
  for (i,R) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed pes[:,i] = real(eigvals(full(hamiltonian))) + Pol(R, α)
    println("$i of $(length(ρ)): H: $(round(Htime,5))s; D: $(round(Dtime,5))s.")
  end
  return pes'
end

function Calc_ves(ρ, Proj, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α)
  pes = Array{Float64}(size(basis,1), length(ρ))
  c = Array{Float64}(size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and chosen Projection (P):")
  for (i,R) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed tmp = eigfact(full(hamiltonian))
    pes[:,i] = real(tmp[:values]) + Pol(R, α)
    v = real(tmp[:vectors])
    Ctime = @elapsed c[:,i] = diag(v'*Proj*v)
    println("$i of $(length(ρ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ctime,5))s.")
  end
  return pes', c'
end

function Calc_pes_pS(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α)
  pes = Array{Float64}(size(basis,1), length(ρ))
  p_s = Array{Float64}(size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and Singlet Projection (P):")
  P_S = ST(basis, 0, 0)
  for (i,R) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed begin
        EigenSystem = eigfact(full(hamiltonian))
        pes[:,i] = real(EigenSystem[:values]) + Pol(R, α)
    end
    Ctime = @elapsed begin
        v = real(EigenSystem[:vectors])
        p_s[:,i] = diag(v'*P_S*v)
    end
    println("$i of $(length(ρ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ctime,5))s.")
  end
  return pes', p_s'
end

function Calc_pes_pF(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, I1, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α)
  pes = Array{Float64}(size(basis,1), length(ρ))
  p_f = Array{Float64}(size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and Hyperfine Projection (P):")
  P_F = Fsquare(basis, I1)
  for (i,R) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed begin
        EigenSystem = eigfact(full(hamiltonian))
        pes[:,i] = real(EigenSystem[:values]) + Pol(R, α)
    end
    Ctime = @elapsed begin
        v = real(EigenSystem[:vectors])
        p_f[:,i] = diag(v'*P_F*v)
    end
    println("$i of $(length(ρ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ctime,5))s.")
  end
  return pes', p_f'
end

function Calc_pes_pL(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α)
  pes = Array{Float64}(size(basis,1), length(ρ))
  p_l = Array{Float64}(size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and Orbital angular momentum Projection (P):")
  P_L = spzeros(Float64, size(basis,1))
  P_L[find(basis[:,2].==3)] += 1
  P_L = dropzeros(spdiagm(P_L))
  # P_F = spzeros(Float64, size(basis,1))
  # P_F[find(sum(basis[:,5:6],2).==2)] += 1
  # P_F = dropzeros(spdiagm(P_F))
  for (i,R) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed begin
        EigenSystem = eigfact(full(hamiltonian))
        pes[:,i] = real(EigenSystem[:values]) + Pol(R, α)
    end
    Ctime = @elapsed begin
        v = real(EigenSystem[:vectors])
        p_l[:,i] = diag(v'*P_L*v)
    end
    println("$i of $(length(ρ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ctime,5))s.")
  end
  return pes', p_l'
end

function Calc_pes_pj(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α)
  pes = Array{Float64}(size(basis,1), length(ρ))
  c0 = Array{Float64}(size(basis,1), length(ρ))
  c1 = Array{Float64}(size(basis,1), length(ρ))
  c2 = Array{Float64}(size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D), and chosen Projection (P):")
  for (i,R) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Ctime = @elapsed begin
        jprojector0 = J_Projector(0, R, Delta, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
        jprojector1 = J_Projector(1, R, Delta, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
        jprojector2 = J_Projector(2, R, Delta, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm)
    end
    Dtime = @elapsed EigenSystem = eigfact(full(hamiltonian))
    pes[:,i] = real(EigenSystem[:values]) + Pol(R, α)
    v = real(EigenSystem[:vectors])
    c0[:,i] = diag(v'*jprojector0*v)
    c1[:,i] = diag(v'*jprojector1*v)
    c2[:,i] = diag(v'*jprojector2*v)
    println("$i of $(length(ρ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s. P: $(round(Ctime,5))s.")
  end
  return pes', c0', c1', c2'
end

function Calc_pes_ev(ρ, V0, Delta, n, LSJ, basis, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12, α)
  pes = Array{Float64}(size(basis,1), length(ρ))
  ev = Array{Float64}(size(basis,1), size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H) and Diagonalisation (D):")
  for (i,R) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(R, V0, Delta, n, LSJ, nlj, Ls, δRb, Pnlj, Pl, CG0, CGp, CGm, aS0, aS1, aP0, aP10, aP11, aP12)
    Dtime = @elapsed begin
        EigenSystem = eigfact(full(hamiltonian))
        pes[:,i] = real(EigenSystem[:values]) + Pol(R, α)
        ev[:,:,i] = real(EigenSystem[:vectors])
    end
    println("$i of $(length(ρ)): H: $(round(Htime,5))s. D: $(round(Dtime,5))s.")
  end
  return pes', ev
end

# Projectors ___________________________________________________________________
function Gsquare(basis, I1) # Total spin projector needs to be implemented, this is just a copy of Ksquare
    if I1 != 0
        mj = basis[:,4]
        m2 = basis[:,5]
        mI = basis[:,6]
        b = zeros(size(basis,1), size(basis,1))
        for i in eachindex(mI), j in eachindex(mI)
            if basis[i,1:3] == basis[j,1:3] && mI[i] + mj[i] + m2[i] == mI[j] + mj[j] + m2[j]
                for F in I1-0.5:I1+0.5, G in abs(I1-1):I1+1
                    b[i,j] += K*(K+1)* W3j(0.5,m1[i], 0.5,m2[i], S,m1[i]+m2[i])* W3j(S,m1[i]+m2[i], I1,mI[i], K,m1[i]+m2[i]+mI[i])* W3j(0.5,m1[j], 0.5,m2[j], S,m1[j]+m2[j])* W3j(S,m1[j]+m2[j], I1,mI[j], K,m1[j]+m2[j]+mI[j])
                end
            end
        end
        return sparse(b)
    else
        return speye(size(basis)[1]) - ST(basis, 0, 0)
    end
end

function Ksquare(basis, I1) # Total spin projector excluding Rydberg orbital angular momentum not tested yet
    if I1 != 0
        m1 = basis[:,4]
        m2 = basis[:,5]
        mI = basis[:,6]
        b = zeros(size(basis,1), size(basis,1))
        for i in eachindex(mI), j in eachindex(mI)
            if basis[i,1:3] == basis[j,1:3] && mI[i] + m1[i] + m2[i] == mI[j] + m1[j] + m2[j]
                for S in 0:1, K in abs(I1-1):I1+1
                    b[i,j] += K*(K+1)* W3j(0.5,m1[i], 0.5,m2[i], S,m1[i]+m2[i])* W3j(S,m1[i]+m2[i], I1,mI[i], K,m1[i]+m2[i]+mI[i])* W3j(0.5,m1[j], 0.5,m2[j], S,m1[j]+m2[j])* W3j(S,m1[j]+m2[j], I1,mI[j], K,m1[j]+m2[j]+mI[j])
                end
            end
        end
        return sparse(b)
    else
        return speye(size(basis)[1]) - ST(basis, 0, 0)
    end
end

function Fsquare(basis, I1) # Total nulcear spin projector
    if I1 != 0
        mS = basis[:,5]
        mI = basis[:,6]
        b = zeros(size(basis,1), size(basis,1))
        for i in eachindex(mI), j in eachindex(mI)
            if basis[i,1:4] == basis[j,1:4] && mI[i] + mS[i] == mI[j] + mS[j]
                for F in I1-0.5:I1+0.5
                    b[i,j] += F*(F+1)* W3j(0.5,mS[i], I1,mI[i], F,mI[i]+mS[i])* W3j(0.5,mS[j], I1,mI[j], F,mI[j]+mS[j])
                end
            end
        end
        return sparse(b)
    else
        return spdiagm(abs.(basis[:,5]))
    end
end

function FmF(basis, θ, I1)
    # S2 component
    ms = basis[:,5]
    H_S2 = Cos(θ)*spdiagm(ms) - Sin(θ)*(ms*ms'-1/4).*S2filt(basis)

    # I2 component
    mi = basis[:,6]
    I_p = Sin(θ)*0.5sqrt.(I1*(I1+1)-mi.*(mi+1))
    I_m = Sin(θ)*0.5sqrt.(I1*(I1+1)-mi.*(mi-1))
    H_I = Cos(θ)*spdiagm(mi)

    for i in eachindex(I_p), j in eachindex(I_p)
        if basis[i,1:5] == basis[j,1:5]
            if mi[i]+1 == mi[j]
                H_I[i,j] = I_p[i]
            elseif mi[i]-1 == mi[j]
                H_I[i,j] = I_m[i]
            end
        end
    end
    return Symmetric(full(H_S2 + H_I))
end

function JmJ(basis, θ)
    # J component
    J = basis[:,3]
    MJ = basis[:,4]
    J_p = Sin(θ)*0.5sqrt.(J.*(J+1)-MJ.*(MJ+1))
    J_m = Sin(θ)*0.5sqrt.(J.*(J+1)-MJ.*(MJ-1))
    H_J = Cos(θ)*spdiagm(MJ)

    for i in eachindex(J_p), j in eachindex(J_p)
      if basis[i,1:3] == basis[j,1:3] && basis[i,5:6] == basis[j,5:6]
          if MJ[i]+1 == MJ[j]
            H_J[i,j] = J_p[i]
          elseif MJ[i]-1 == MJ[j]
            H_J[i,j] = J_m[i]
          end
      end
    end
    return Symmetric(full(H_J))
end

function Omega(basis, θ, I1)
    # J component
    J = basis[:,3]
    MJ = basis[:,4]
    J_p = Sin(θ)*0.5sqrt.(J.*(J+1)-MJ.*(MJ+1))
    J_m = Sin(θ)*0.5sqrt.(J.*(J+1)-MJ.*(MJ-1))
    H_J = Cos(θ)*spdiagm(MJ)

    # S1 component
    C1(l,j,mj) = cg(l,j,mj,-0.5)
    C2(l,j,mj) = cg(l,j,mj,0.5)
    cminus = map(C1, basis[:,2], basis[:,3], basis[:,4])
    cplus = map(C2, basis[:,2], basis[:,3], basis[:,4])
    (z, p) = S1filt(basis)
    S1p = cplus*cminus'.*p
    H_S1 = Cos(θ)/2*(cplus*cplus' - cminus*cminus').*z + Sin(θ)/2*(S1p + S1p')
    # H_S1 = Cos(θ)*diagm((-1).^(L+0.5-J).*MJ./(1+2L)) + Sin(θ)/2*(S1p + S1p')

    # S2 component
    ms = basis[:,5]
    H_S2 = Cos(θ)*spdiagm(ms) - Sin(θ)*(ms*ms'-1/4).*S2filt(basis)

    # I2 component
    mi = basis[:,6]
    I_p = Sin(θ)*0.5sqrt.(I1*(I1+1)-mi.*(mi+1))
    I_m = Sin(θ)*0.5sqrt.(I1*(I1+1)-mi.*(mi-1))
    H_I = Cos(θ)*spdiagm(mi)

    for i in eachindex(I_p), j in eachindex(I_p)
        if basis[i,1:3] == basis[j,1:3]
            if basis[i,5:6] == basis[j,5:6]
                if MJ[i]+1 == MJ[j]
                    H_J[i,j] = J_p[i]
                elseif MJ[i]-1 == MJ[j]
                    H_J[i,j] = J_m[i]
                end
            elseif basis[i,4:5] == basis[j,4:5]
                if mi[i]+1 == mi[j]
                    H_I[i,j] = I_p[i]
                elseif mi[i]-1 == mi[j]
                    H_I[i,j] = I_m[i]
                end
            end
        end
    end
    return Symmetric(full(H_J + H_S2 + H_I))
end

function ST(basis, S, ms) # Projector |S, MS⟩⟨S, MS|
  L = basis[:,2]
  J = basis[:,3]
  MJ = basis[:,4]
  M2 = basis[:,5]

  P = zeros(Float64, length(J), length(J))
  for i in eachindex(J), j in eachindex(J)
    if basis[i,1:2] == basis[j,1:2] && basis[i,6] == basis[j,6]
      for m in -0.5:0.5, n in -0.5:0.5
        P[i,j] += W3j(L[i], MJ[i]-m, 0.5, m, J[i], MJ[i])*W3j(L[i], MJ[i]-m, 0.5, n, J[j], MJ[j])*W3j(0.5, m, 0.5, M2[i], S, ms)*W3j(0.5, n, 0.5, M2[j], S, ms)
      end
    end
  end

  return sparse(P)
end

function STθ(basis)
  L = basis[:,2]
  J = basis[:,3]
  MJ = basis[:,4]
  M2 = basis[:,5]

  P = zeros(Float64, length(J), length(J))
  for i in eachindex(J), j in eachindex(J)
    if basis[i,1:2] == basis[j,1:2] && basis[i,6] == basis[j,6]
      for m in -0.5:0.5, n in -0.5:0.5, ms in -1:1
        P[i,j] += ms*W3j(L[i], MJ[i]-m, 0.5, m, J[i], MJ[i])*W3j(L[i], MJ[i]-m, 0.5, n, J[j], MJ[j])*W3j(0.5, m, 0.5, M2[i], 1, ms)*W3j(0.5, n, 0.5, M2[j], 1, ms)
      end
    end
  end

  return sparse(P)
end

function STθ(basis, θ) # Projector |MS⟩⟨MS|
  L = basis[:,2]
  J = basis[:,3]
  MJ = basis[:,4]
  M2 = basis[:,5]

  P = zeros(Float64, length(J), length(J))
  Pp = zeros(Float64, length(J), length(J))
  Pm = zeros(Float64, length(J), length(J))
  for i in eachindex(J), j in eachindex(J)
    if basis[i,1:2] == basis[j,1:2] && basis[i,6] == basis[j,6]
      for m in -0.5:0.5, n in -0.5:0.5, ms in -1:1
        P[i,j] += ms*W3j(L[i], MJ[i]-m, 0.5, m, J[i], MJ[i])*W3j(L[i], MJ[i]-m, 0.5, n, J[j], MJ[j])*W3j(0.5, m, 0.5, M2[i], 1, ms)*W3j(0.5, n, 0.5, M2[j], 1, ms)
        Pp[i,j] += sqrt(2-ms*(ms+1))*W3j(L[i], MJ[i]-m, 0.5, m, J[i], MJ[i])*W3j(L[i], MJ[i]-m, 0.5, n, J[j], MJ[j])*W3j(0.5, m, 0.5, M2[i], 1, ms+1)*W3j(0.5, n, 0.5, M2[j], 1, ms)
        Pm[i,j] += sqrt(2-ms*(ms-1))*W3j(L[i], MJ[i]-m, 0.5, m, J[i], MJ[i])*W3j(L[i], MJ[i]-m, 0.5, n, J[j], MJ[j])*W3j(0.5, m, 0.5, M2[i], 1, ms-1)*W3j(0.5, n, 0.5, M2[j], 1, ms)
      end
    end
  end

  return Cos(θ)*sparse(P) + Sin(θ)/2*(Pp + Pm)
end

# Misc _________________________________________________________________________
function Colorize(P, ρ, v, basis)
  c = Array{Float64}(size(basis,1), length(ρ))
  for i in eachindex(ρ)
    c[:,i] = diag(v[:,:,i]'*P*v[:,:,i])
  end
  return c'
end

#CC3 = ColorGradient([RGB(1,1,0),RGB(0,1,1),RGB(1,0,1)]); # converging
#CC4 = ColorGradient([RGB(1,0,0),RGB(1,1,0),RGB(0,1,1),RGB(1,0,1)]); # converging
#CC5 = ColorGradient([RGB(1,0,0),RGB(1,0.5,0),RGB(1,1,0),RGB(0,1,1),RGB(0,0,1)]); # converging
#CC5r = ColorGradient([RGB(0,0,1),RGB(0,1,1),RGB(1,1,0),RGB(1,0.5,0),RGB(1,0,0)]); # converging
#CC6 = ColorGradient([RGB(0,0,0),RGB(1,0,0),RGB(1,1,0),RGB(0,1,0),RGB(0,0,1),RGB(0.5,0,0.5)]); # converging
#CC10 = ColorGradient([RGB(0,0,0),RGB(1,0,0),RGB(1,0.5,0),RGB(1,1,0),RGB(0,0.5,0),RGB(0,1,0.5),RGB(0,1,1),RGB(0,0,1),RGB(0.5,0,0.5),RGB(1,0,1)]);

function rel(ρ, PES, lower, upper) # find relevant values for plotting
  r = length(ρ)
  A = Float64[]
  for (i, p) in enumerate(PES)
    if lower < p < upper
      push!(A, ρ[ifelse(i%r!=0, i%r, r)], p)
    end
  end
  return reshape(A, 2, round(Int, length(A)/2))'
end

function findlocalmaxima(signal::Vector, threshold::Real)
   inds = Int[]
   if length(signal) > 1
       if signal[1] > signal[2] && signal[1] >= threshold
           push!(inds,1)
       end
       for i = 2:length(signal)-1
           if signal[i-1] < signal[i] > signal[i+1] && signal[i] >= threshold
               push!(inds,i)
           end
       end
       if signal[end] > signal[end-1] && signal[end] >= threshold
           push!(inds,length(signal))
       end
   end
   return inds
end

function findlocalminima(signal::Vector, threshold::Real)
   inds = Int[]
   if length(signal) > 1
       if signal[1] < signal[2] && signal[1] >= threshold
           push!(inds,1)
       end
       for i = 2:length(signal)-1
           if signal[i-1] > signal[i] < signal[i+1] && signal[i] >= threshold
               push!(inds,i)
           end
       end
       if signal[end] < signal[end-1] && signal[end] >= threshold
           push!(inds,length(signal))
       end
   end
   return inds
end

function Base.sortperm(A::AbstractMatrix, dim::Integer)
    P = mapslices(sortperm, A, dim)
    if dim == 1
        for j = 1:size(P,2)
            offset = (j-1) * size(P,1)
            for i = 1:size(P,1)
                P[i,j] += offset
            end
        end
    else # if dim == 2
        for j = 1:size(P,2)
            for i = 1:size(P,1)
                P[i,j] = (P[i,j] - 1) * size(P,1) + i
            end
        end
    end
    return P
end

function Gauss(x, E, FC, Δ)
    pdf = zeros(size(x,1))
    for j in eachindex(x), i in eachindex(E)
        if x[1] < E[i] < x[end]
            pdf[j] +=  FC[i]/(Δ*sqrt(2π))*exp(-1/2*((x[j]-E[i])/Δ)^2)
        else
            nothing
        end
    end
    return pdf
end

# ______________________________________________________________________________
Mod(i, s) = ifelse(i%s!=0, i%s, s)

function findPES(Tol, ρ, PES, P_Ω)
    b = Float64[]
    for i in eachindex(PES[:,1]), j in eachindex(PES[1,:])
        if Tol <= P_Ω[i,j] <= 1
            push!(b, ρ[Mod.(i, size(PES,1))], PES[i,j])
        end
    end
    return reshape(b, 2, round(Int, length(b)/2))'
end

function selectPES(Tol, ρ, PES, P_Ω)
    sPES = Array{Float64}(length(ρ));
    tmpPES = findPES(Tol, ρ, PES, P_Ω)
    double = find(x -> x .== 0, tmpPES[2:end,1]-tmpPES[1:end-1,1])
    for j in reverse(double)
        if tmpPES[j,2] < tmpPES[j+1,2]
            tmpPES = vcat(tmpPES[1:j,:], tmpPES[j+2:end,:])
        elseif tmpPES[j,2] >= tmpPES[j+1,2]
            tmpPES = vcat(tmpPES[1:j-1,:], tmpPES[j+1:end,:])
        end
    end
    return sPES
end

# ______________________________________________________________________________
function Vionpair(ρ)
    αp = 9.11 # Rb^{+} polarizability
    αm = 526 # Rb^{-} polarizability
    EA = -0.01786 # Rubidium electron affinity
    b = Float64[]
    f(r) = -1/r -(αp + αm)/2r^4
    for r in ρ
        push!(b, f(r))
    end
    return b + EA
end

function Vionpair(ρ,l)
    αp = 9.11 # Rb^{+} polarizability
    αm = 526 # Rb^{-} polarizability
    EA = -0.01786 # Rubidium electron affinity
    b = Float64[]
    f(r) = -1/r +l*(l+1)/158426r^2 -(αp + αm)/2r^4
    for r in ρ
        push!(b, f(r))
    end
    return b + EA
end

function lowStates(maxN::Int64)
    b = Float64[]
    c = Float64[]
    δRb = [3.1311804; 2.65; 1.35; 0.016; 0.; 0.; 0.; 0.; 0.]
    for n in 4:maxN+3
        if n == 4
            for l in 2:3
                for m in -l:l
                    push!(b, n, l, m)
                    push!(c, n-δRb[l+1])
                end
            end
        elseif 4 < n <= maxN
            for l in 0:n-1
                for m in -l:l
                    push!(b, n, l, m)
                    push!(c, n-δRb[l+1])
                end
            end
        elseif n > maxN
            for l in 0:maxN+3-n
                for m in -l:l
                    push!(b, n, l, m)
                    push!(c, n-δRb[l+1])
                end
            end
        else
            println("Error: Upward basis construction mathod only valid up to n=9 hydrogenic states.")
            println("n=$n selected. Please select n=9 or lower.")
        end
    end
    p = sortperm(c)
    bb = reshape(b, 3, round(Int, length(b)/3))'
    return (bb[p,:], c[p])
end

function boundEckart(ρ, R)
    ρ0 = 9.004786
    V0 = 0.061675
    λ0 = (sqrt(8*V0*ρ0^2+1)-1)/4
    b = Float64[]
    f(r) = 1/cosh(r/ρ0)^(2λ0)*sinh(r/ρ0)#*sf_hyperg_2F1(0, -2λ0+1, 1.5, -sinh(ρ/ρ0)^2)
    for r in ρ
        push!(b, f(r-R))
    end
    # normEckart = 5.47358356076463
    return b/sqrt(sum(b.^2))
end
