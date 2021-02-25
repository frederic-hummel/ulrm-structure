DataPath = "/home/frederic/Documents/julia/data";

CT(Ψ, Ξ) = sparse(real(Ψ*ctranspose(Ξ)))

function Cesium(Pwave, J) # phaseshift to scattering length
    aSP = readcsv("$DataPath/CsPhases_3.csv")

    a0 = extrapolate(interpolate((aSP[:,1], ), zeros(length(aSP[:,1])), Gridded(Linear())), Flat())
    aS = extrapolate(interpolate((aSP[:,1], ), aSP[:,3], Gridded(Linear())), Flat())
    aP10 = extrapolate(interpolate((aSP[:,1], ), aSP[:,5], Gridded(Linear())), Flat())
    aP11 = extrapolate(interpolate((aSP[:,1], ), aSP[:,6], Gridded(Linear())), Flat())
    aP12 = extrapolate(interpolate((aSP[:,1], ), aSP[:,7], Gridded(Linear())), Flat())
    if Pwave
        if J == 0
            aP = aP10
        elseif J == 1
            aP = aP11
        elseif J == 2
            aP = aP12
        else
            aP1 = aP11
            println("J=1 autoselected. Choose J=0,1,2 to select a phaseshift for p-wave scattering.")
        end
    else
        aP = a0
    end
    return (aS, aP)
end

function Rubidium(Pwave)
    aSP = readcsv("$DataPath/Rbscatl_fit_c=1.0.csv")

    cut = 201
    aSP[1:cut,3:end] = repeat(aSP[cut,3:end]', inner=(cut,1))

    a0 = extrapolate(interpolate((aSP[:,1], ), zeros(length(aSP[:,1])), Gridded(Linear())), Flat())
    aS = extrapolate(interpolate((aSP[:,1], ), aSP[:,2], Gridded(Linear())), Flat())
    aP = extrapolate(interpolate((aSP[:,1], ), aSP[:,3], Gridded(Linear())), Flat())

    if Pwave == false
        aP = a0
    end

    return (aS, aP)
end

function Strontium() # phaseshift to scattering length
    aS0 = -13.2
    aP0 = 8.4
    αk = 186
    K = collect(0:0.001:0.13)
    aSk = aS0 + K*αk*π/3
    aS = extrapolate(interpolate((K, ), aSk, Gridded(Linear())), Linear())
    aPk = fill(aP0, length(K))
    aP = extrapolate(interpolate((K, ), aPk, Gridded(Linear())), Flat())
    δ = Float64[]
    push!(δ, 3.371)
    return (aS, aP, δ)
end

k(R,n) = real(sqrt.(complex(2(-1/2n^2 + 1/R)))) # kinetic energy
Pol(α, R) = -α./2R.^4 # polarisation potential

function QDaP(Cs)
    if Cs
        δRb = [4.049325; 3.57; 2.46; 0.0334] # quantum defects
        α = 402.2 # polarisability
    else
        δRb = [3.1311804; 2.65; 1.34; 0.0165192] # quantum defects
        α = 319.5 # polarisability
    end
    return (δRb, α)
end

# Define basis set _____________________________________________________________
function δnlj(l::Int, δRb) # quantum defect
  if l<size(δRb)[1]
    δRb[l+1,1]
  else
    0.
  end
end

function δnlj(l::Array, δRb)
  b = Float64[]
  for i in eachindex(l)
    if l[i] < size(δRb)[1]
      push!(b, δRb[l[i]+1,1])
    else
      push!(b, 0.)
    end
  end
  return b
end

function pos(n, l, list)
  pos = Int[]
  for i in eachindex(list[:,1])
    if list[i,1:2] == [n,l]
      push!(pos, i)
    end
  end
  return pos
end

function Pos(basis::Array, nl::Array, lm::Array, m_l::Array)

  pnl = Int[]
  for i in eachindex(basis[:,1])
    append!(pnl, pos(basis[i,1], basis[i,2], nl))
  end

  plm = Int[]
  for i in eachindex(basis[:,1])
    append!(plm, pos(basis[i,2], abs.(basis[i,3]), lm))
  end

  pm_l = Int[]
  for i in eachindex(basis[:,1])
    append!(pm_l, find(m_l.==basis[i,3]))
  end

  Rbasis = sparse(pnl, eachindex(pnl), ones(Int, length(pnl)))
  Pbasis = sparse(plm, eachindex(plm), ones(Int, length(plm)))
  Abasis = sparse(pm_l, eachindex(pm_l), ones(Int, length(pm_l)))

  return (Rbasis, Pbasis, Abasis)
end

function Hydrogenic(n, Δnp, Δnm, m_max) # hydrogenic states with hyperfine interaction
  b = Int64[]
  for i in n-Δnm:n+Δnp
    for j in 4:i-1
      for k in -min(j, m_max):min(j, m_max)
        push!(b, i, j, k)
      end
    end
  end
  return reshape(b, 3, round.(Int, length(b)/3))'
end

function QDefect(n, Δnp, Δnm, m_max, δRb) # quantum defect states with hyperfine interaction
  b = Int64[]
  for i in n-Δnm:n+Δnp
    for j in 0:3
      for k in -min(j, m_max):min(j, m_max)
        push!(b, i, j, k)
      end
    end
  end
  b = reshape(b, 3, round.(Int, length(b)/3))'
  for i in eachindex(b[:,1])
    b[i,1] += round.(δnlj(b[i,2], δRb), RoundDown)
  end
  return b
end

# Basis functions ______________________________________________________________
function Ψ(ρ, nl, lm, m_l, δRb) # accumulated basis function calculation

  Ψ_ρ_intern(ρ) = Ψ_ρ(ρ, nl[:,1]-δnlj(nl[:,2], δRb), nl[:,2])

  R = zeros(size(nl)[1], 3)
  R[:,1] = Ψ_ρ_intern(ρ)
  R[:,2] = derivative(Ψ_ρ_intern, ρ)
  R[:,3] = R[:,1]/ρ

  P = zeros(size(lm)[1], 3)
  P[:,1] = Ψ_θ(lm)
  P[:,2] = DΨ_θ(lm)
  P[:,3] = θΨ_θ(lm)

  A = complex(zeros(size(m_l)[1], 2))
  A[:,1] = Ψ_ϕ(m_l)
  A[:,2] = DΨ_ϕ(m_l)

  return (R, P, A)
end

function Ψ_ρ(ρ, n, l) # radial basis function
  # in terms of Whittaker functions for the quantum defect states
  f(N, L) = 1/ρ* 1/N* 1/sqrt.(sf_gamma(N+L+1)* sf_gamma(N-L))* exp(-ρ/N)* (2ρ/N)^(L+1)* sf_hyperg_U(L+1-N, 2(L+1), 2ρ/N)
  # and in terms of hydrogenic functions for regular states
  g(N, L) = 2/N^2* sqrt.(sf_gamma(N-L)/sf_gamma(N+L+1))* exp(-ρ/N)* (2ρ/N)^L* sf_laguerre_n(N-L-1, 2L+1, 2ρ/N)
  b = Float64[]
  for (i, l0) in enumerate(l)
    if l0 < 4
      push!(b, f(n[i], l0))
    else
      push!(b, g(round.(Int, n[i]), round.(Int, l0)))
    end
  end
  return b
end

function Ψ_θ(lm) # azimutal basis function
  function g(l, m)
      if m >= 0
          sf_legendre_sphPlm(l, m, 1)
      else
          (-1)^(-m)* sf_legendre_sphPlm(l, -m, 1)
      end
  end
  return round.(map(g, lm[:,1], lm[:,2]), 15)
end

Ψ_ϕ(m) = ones(Float64, size(m,1))
DΨ_ϕ(m) = im* m

function DΨ_θ(lm) # azimutal gradient basis function
  function g(l, m)
    if l==0
      0.
    elseif m==0
      sqrt.((2l+1)/4π)* sf_legendre_Plm(l, 1, 1)
    elseif m==l
      -0.5* sqrt.((2l+1)/4π* sf_gamma(2l+1))* l* sf_legendre_Plm(l, l-1, 1)
    else
      -0.5* sqrt.((2l+1)/4π* sf_gamma(l-m+1)/sf_gamma(l+m+1))* ((l+m)*(l-m+1)* sf_legendre_Plm(l, m-1, 1)- sf_legendre_Plm(l, m+1, 1))
    end
  end
  map(g, lm[:,1], lm[:,2])
end

function θΨ_θ(lm) # polar gradient azimutal basis function
  function g(l, m)
    if m==0
      0.
    else
      -sqrt.((2l+1)/4π)* sqrt.(sf_gamma(l-m+1)/sf_gamma(l+m+1))* 1/2m* ((l-m+1)*(l-m+2)* sf_legendre_Plm(l+1, m-1, 1)+ sf_legendre_Plm(l+1, m+1, 1))
    end
  end
  map(g, lm[:,1], lm[:,2])
end

# Set up B field Hamiltonian __________________________________________________
function Sin(ω)
    if ω == π || ω == -π
        0.
    else
        sin(ω)
    end
end

function Cos(ω)
    if ω == π/2 || ω == -π/2
        0.
    else
        cos(ω)
    end
end

function ZeemanTerm(basis, ω, B)
    if B != 0
        l = basis[:,2]
        ml = basis[:,3]
        L_p = Sin(ω)*0.5sqrt.(l.*(l+1)-ml.*(ml+1))
        L_m = Sin(ω)*0.5sqrt.(l.*(l+1)-ml.*(ml-1))
        H_L = Cos(ω)*spdiagm(ml)

        for i in eachindex(L_p), j in eachindex(L_p)
            if basis[i,1:2] == basis[j,1:2]
                if ml[i]+1 == ml[j]
                    H_L[i,j] = L_p[i]
                elseif ml[i]-1 == ml[j]
                    H_L[i,j] = L_m[i]
                end
            end
        end

        out = B/2.35051756758E9/2*H_L # in units auf Gauß = 10^-4 Tesla
    else
        out = spzeros(size(basis)[1], size(basis)[1])
    end
    return out
end

function StarkTerm(basis, ω, F, n)
    Radial = spzeros(size(basis)[1], size(basis)[1])
    AngularZ = spzeros(size(basis)[1], size(basis)[1])
    AngularX = spzeros(size(basis)[1], size(basis)[1])
    if F != 0
        Radial = sparse(readcsv("$DataPath/Radial$n.csv"))
        AngularZ = sparse(readcsv("$DataPath/AngularZ$n.csv"))
        AngularX = sparse(readcsv("$DataPath/AngularX$n.csv"))
    end
    return F/5.14220674763E11*Radial.*(Cos(ω)*AngularZ + Sin(ω)*AngularX) # in units of V/m
end

function StarkTermDebug(calc, basis, ω, F)
    Radial = spzeros(size(basis)[1], size(basis)[1])
    AngularZ = spzeros(size(basis)[1], size(basis)[1])
    AngularX = spzeros(size(basis)[1], size(basis)[1])
    if calc
        n = basis[:,1]
        l = basis[:,2]
        m = basis[:,3]
        for i in eachindex(n), j in eachindex(n)
            Radial[i,j] = StarkRadialIntegrals(n[i], n[j], l[i], l[j])
            AngularZ[i,j] = StarkAngularZIntegrals(l[i], l[j], m[i], m[j])
            AngularX[i,j] = StarkAngularXIntegrals(l[i], l[j], m[i], m[j])
        end
    else
        Radial = sparse(readcsv("$DataPath/Radial.csv"))
        AngularZ = sparse(readcsv("$DataPath/AngularZ.csv"))
        AngularX = sparse(readcsv("$DataPath/AngularX.csv"))
    end
    return Radial, AngularZ, AngularX
end

function CalcStarkTerm(basis, ω, F) # Does not converge for n>10
    Radial = spzeros(size(basis)[1], size(basis)[1])
    AngularZ = spzeros(size(basis)[1], size(basis)[1])
    AngularX = spzeros(size(basis)[1], size(basis)[1])
    if F != 0
        n = basis[:,1]
        l = basis[:,2]
        m = basis[:,3]
        for i in eachindex(n), j in eachindex(n)
            Radial[i,j] = StarkRadialIntegrals(n[i], n[j], l[i], l[j])
            AngularZ[i,j] = StarkAngularZIntegrals(l[i], l[j], m[i], m[j])
            AngularX[i,j] = StarkAngularXIntegrals(l[i], l[j], m[i], m[j])
        end
    end
    return F/5.14220674763E11*Radial.*(Cos(ω)*AngularZ + Sin(ω)*AngularX)
end

function StarkAngularZIntegrals(li, lj, mi, mj)
    if mi == mj && abs(li-lj) == 1
        l = min(li, lj)
        return sqrt((l+mi+1)*(l-mi+1)/(2l+1)/(2l+3))
    else
        return 0.
    end
end

function StarkAngularXIntegrals(li, lj, mi, mj)
    if mi+1 == mj && li-1 == lj
        sqrt((li-mi)*(li-mi-1)/(4li^2-1))/2
    elseif mi-1 == mj && li-1 == lj
        -sqrt((li+mi)*(li+mi-1)/(4li^2-1))/2
    elseif mi+1 == mj && li+1 == lj
        sqrt((2li+1)*(li+mi+2)/(2li+3)/(li+mi+1))*((li-mi)/(2li+1)-1)/2
    elseif mi-1 == mj && li+1 == lj
        -sqrt((2li+1)*(li-mi+2)/(2li+3)/(li-mi+1))*((li+mi)/(2li+1)-1)/2
    else
        0.
    end
end

function StarkRadialIntegrals(ni, nj, li, lj) # Does not converge for n>10
    if abs(li-lj) == 1
        Prefactor = 4* 2^(li+lj)/ ni^(li+2)/ nj^(lj+2)* sqrt(sf_fact(ni-li-1)* sf_fact(nj-lj-1)/ sf_fact(ni+li)/ sf_fact(nj+lj))
        Sum = 0.
        for mi in 0:(ni-li-1), mj in 0:(nj-lj-1)
            Sum += (-2)^(mi+mj)/ sf_fact(mi)/ sf_fact(mj)/ ni^mi/ nj^mj* sf_fact(li+lj+mi+mj+3)/ ((ni+nj)/ni/nj)^(4+li+lj+mi+mj)* binomial(ni+li, ni-li-mi-1)* binomial(nj+lj, nj-lj-mj-1)
        end
        return Prefactor*Sum
    else
        return 0.
    end
end

# Set coordinate dependant matrices ____________________________________________
function Ψmat(r, nl, lm, m_l, Rbasis, Pbasis, Abasis, δRb)
  (R, P, A) = Ψ(r, nl, lm, m_l, δRb)

  BG_R = Rbasis'*R
  BG_P = Pbasis'*P
  BG_A = Abasis'*A

  BG = complex(zeros(size(BG_R)[1], 4))

  BG[:,1] = BG_R[:,1].*BG_P[:,1].*BG_A[:,1]
  BG[:,2] = BG_R[:,2].*BG_P[:,1].*BG_A[:,1]
  BG[:,3] = BG_R[:,3].*BG_P[:,2].*BG_A[:,1]
  BG[:,4] = BG_R[:,3].*BG_P[:,3].*BG_A[:,2]

  return BG
end

function H(r, n, nl, lm, m_l, Rbasis, Pbasis, Abasis, δRb, V0, aS, aP, H_B)
  BG = Ψmat(r, nl, lm, m_l, Rbasis, Pbasis, Abasis, δRb)

  VS = CT(BG[:,1], BG[:,1])
  VP = CT(BG[:,2], BG[:,2]) + CT(BG[:,3], BG[:,3]) + CT(BG[:,4], BG[:,4])
  v = k(r,n)

  return V0 + 2π*aS[v]*VS + 6π*aP[v]*VP + H_B
end

function Calc_pes(ρ, n, nl, lm, m_l, basis, Rbasis, Pbasis, Abasis, δRb, α, V0, aS, aP, H_B)
  pes = Array{Float64}(size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  #println("Timing of the Hamiltonian setup (H) and Diagonalisation (D):")
  for (i,r) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(r, n, nl, lm, m_l, Rbasis, Pbasis, Abasis, δRb, V0, aS, aP, H_B)
    Dtime = @elapsed pes[:,i] = eigvals(full(hamiltonian)) + Pol(α, r)
    #println("$i of $(length(ρ)): H: $(round.(Htime,5))s. D: $(round.(Dtime,5))s.")
  end
  return pes'
end

function Calc_ves(ρ, n, nl, lm, m_l, basis, Rbasis, Pbasis, Abasis, δRb, α, V0, aS, aP, H_B)
  pes = Array{Float64}(size(basis,1), length(ρ))
  cev = Array{Float64}(size(basis,1), size(basis,1), length(ρ))
  println("Starting calculation for $(size(basis,1)) basis states at $(length(ρ)) grid points.")
  println("Timing of the Hamiltonian setup (H), Diagonalisation (D):")
  for (i,r) in enumerate(ρ)
    Htime = @elapsed hamiltonian = H(r, n, nl, lm, m_l, Rbasis, Pbasis, Abasis, δRb, V0, aS, aP, H_B)
    Dtime = @elapsed EigSys = eigfact(full(hamiltonian))
    pes[:,i] = EigSys[:values] + Pol(α, r)
    cev[:,:,i] = EigSys[:vectors]
    println("$i of $(length(ρ)): H: $(round.(1000Htime,2))μs. D: $(round.(1000Dtime,2))μs.")
  end
  return pes', cev
end

#CC5 = ColorGradient([RGB(1,0,0),RGB(1,0.5,0),RGB(1,1,0),RGB(0,1,1),RGB(0,0,1)]); # converging
#CC5r = ColorGradient([RGB(0,0,1),RGB(0,1,1),RGB(1,1,0),RGB(1,0.5,0),RGB(1,0,0)]); # converging
#CC6 = ColorGradient([RGB(0,0,0),RGB(1,0,0),RGB(1,1,0),RGB(0,1,0),RGB(0,0,1),RGB(0.5,0,0.5)]); # converging
