using BenchmarkTools, Dierckx, LinearAlgebra, NLsolve, Roots

"""    ndgrid(x)
return the grid vectors x1, x2 to produce an 2-dimensional full grid """
function ndgrid(A_vec::AbstractArray{T}) where T <: Float64
    X = [i for i in A_vec, j in A_vec]
    Y = [j for i in A_vec, j in A_vec]
    return X, Y
end

"""    soleq3(x)
return the solution of the third degree equation from 'x' """
function soleq3(a::AbstractArray{T}) where T <: Complex
    Res_R = Array{Complex}(zeros(N, N, 3))
    p = -a[:, :, 1] .^ 2 / 3 + a[:, :, 2]
    q = a[:, :, 1] .^ 3 / 13.5 - a[:, :, 1] .* a[:, :, 2] / 3 + a[:, :, 3]
    Q = p .^ 3 / 27 + q .^ 2 / 4
    A = (-q / 2 + sqrt.(Q)) .^ (1 / 3)
    B = -p ./ A / 3
    ApB = A + B
    AmB = im * (A - B) * sqrt(3) / 2
    Res_R[:, :, 1] = +ApB - a[:, :, 1] / 3
    Res_R[:, :, 2] = -ApB / 2 + AmB - a[:, :, 1] / 3
    Res_R[:, :, 3] = -ApB / 2 - AmB - a[:, :, 1] / 3
    return Res_R
end

"""    fermi(Spectrum, temperature)
return the Furmi function from 'Spectrum' and 'temperature' """
function fermi(Spectrum, Temperature)
    return 1 ./ (exp.(Spectrum / Temperature) .+ 1)
end

"""    aEk(kx, ky)
return three bands of Fermi spectrum in Spin-Fermion Modeд in k space """
function aEk(kx::AbstractArray{T}, ky::AbstractArray{T}) where T <: Float64
    Res_eq_3 = Array{Complex}(zeros(N, N, 3))
    I1  = (1 - p) * I
    csx = cos.(kx .- 2 * α_x); csy = cos.(ky)
    snx = sin.(kx / 2 .- α_x); sny = sin.(ky / 2)
    psi = snx .* sny; sx2 = snx .^ 2; sy2 = sny .^ 2
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(2 * kx .- 4 * α_x) + cos.(2 * ky)) / 2
    Kk  = 3 / 4 .- C1 * gm1
    D11 = ε_p .+ 2 * τ_ * sx2
    D12 = 2 * (τ_ - 2 * t) * psi
    D13 = J * Kk .* snx
    D23 = J * Kk .* sny
    D22 = ε_p .+ 2 * τ_ * sy2
    D33 = (ε_p + 4 * (τ_ - t)) * Kk .- 3 * (3 / 4 * τ_ - t) / 2 +
        + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
        + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    Res_eq_3[:, :, 1] = - D11 +
                        - D22 +
                        - D33 ./ Kk;
    Res_eq_3[:, :, 2] = + D11 .* D22 - D12 .* D12 +
                        + D11 .* D33 ./ Kk - D23 .* D23 ./ Kk +
                        - D13 .* D13 ./ Kk + D22 .* D33 ./ Kk;
    Res_eq_3[:, :, 3] = - D11 .* D22 .* D33 ./ Kk + D11 .* D23 .* D23 ./ Kk +
                        - D12 .* D13 .* D23 ./ Kk - D12 .* D13 .* D23 ./ Kk +
                        + D12 .* D12 .* D33 ./ Kk + D13 .* D13 .* D22 ./ Kk;
    return sort(real(soleq3(Res_eq_3)), dims=3)
end

function aEk_(kx::AbstractArray{T}, ky::AbstractArray{T}) where T <: Float64
    Res_eq_3 = Array{Complex}(zeros(N, N, 3))
    I1 = (1 - p) * I
    csx = cos.(-kx .- 2 * α_x); csy = cos.(-ky)
    snx = sin.(-kx / 2 .- α_x); sny = sin.(-ky / 2)
    sx2 = snx .^ 2; sy2 = sny .^ 2; psi = snx .* sny
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(-2 * kx .- 4 * α_x) + cos(-2 * ky)) / 2
    Kk_ = 3 / 4 .- C1 * gm1
    D11_ = ε_p .+ 2 * τ_ * sx2
    D12_ = 2 * (τ_ - 2 * t) * psi
    D13_ = J * Kk_ .* snx
    D23_ = J * Kk_ .* sny
    D22_ = ε_p .+ 2 * τ_ * sy2
    D33_ = (ε_p + 4 * (τ_ - t)) * Kk_ .- 9 / 8 * τ_ .+ 3 / 2 * t +
         + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
         + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    Res_eq_3[:, :, 1]= D11_ +
                     + D22_ +
                     + D33_ ./ Kk_
    Res_eq_3[:, :, 2]= D11_ .* D22_ - D12_ .* D12_ +
                     + D11_ .* D33_ ./ Kk_ - D23_ .* D23_ ./ Kk_ +
                     - D13_ .* D13_ ./ Kk_ + D22_ .* D33_ ./ Kk_
    Res_eq_3[:, :, 3]= D11_ .* D22_ .* D33_ ./ Kk_ - D11_ .* D23_ .* D23_ ./ Kk_ +
                     + D12_ .* D13_ .* D23_ ./ Kk_ + D12_ .* D13_ .* D23_ ./ Kk_ +
                     - D12_ .* D12_ .* D33_ ./ Kk_ - D13_ .* D13_ .* D22_ ./ Kk_
    return sort(real(soleq3(Res_eq_3)), dims=3)
end

function anpN(float)
    I1  = (1 - p) * I
    csx = cos.(kx .- 2 * α_x); csy = cos.(ky)
    snx = sin.(kx / 2 .- α_x); sny = sin.(ky / 2)
    psi = snx .* sny; sx2 = snx .^ 2; sy2 = sny .^ 2
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(2 * kx .- 4 * α_x) + cos.(2 * ky)) / 2
    Kk  = 3 / 4 .- C1 * gm1
    D11 = ε_p .+ 2 * τ_ * sx2
    D12 = 2 * (τ_ - 2 * t) * psi
    D13 = J * Kk .* snx
    D23 = J * Kk .* sny
    D22 = ε_p .+ 2 * τ_ * sy2
    D33 = (ε_p + 4 * (τ_ - t)) * Kk .- 3 * (3 / 4 * τ_ - t) / 2 +
        + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
        + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    E = aEk(kx, ky)
    Fx = Array{Float64}(zeros(N, N))
    for ll in 1:3
        om = E[:, :, ll]
        S1m = (om - D11)
        S2m = (om - D22)
        S3m = (om - D33 ./ Kk)
        Q11 = S2m .* S3m - D23 .^2 ./ Kk
        Q22 = S1m .* S3m - D13 .^2 ./ Kk
        denom = 1
        for jj in 1:3
            if jj != ll
                denom = denom .* (om - E[:, :, jj])
            end
        end
        denom[findall(denom -> denom==0, denom)] .= 1e-10

        Jh = Q11 + Q22
        Fx = Fx .+ fermi(om .- mu, T) .* Jh ./ denom
    end
    return 2 * sum(real(Fx)) / N ^ 2
end

function amuN(muu::Float64)
    global mu, x
    mu = muu
    return anpN(1) - x
end

function anpSC(X)

    global kx, ky, T, mu, Dlt1, Dlt2, Dlt3, Dlt4, Dlt5, Dlt6
    global ε_p, τ_, t, α_x, x, p, C1, C2, C3, J, I, Up, V2

    """ ====== Finding Deltas: =============================================="""
    nlsolve(aEqB, [.01;.01;.01;.01;.01;.01;])
    """ ====================================================================="""
    d_sym = (cos.(kx)-cos.(ky))
    D14k = d_sym * Dlt1
    D25k = d_sym * Dlt2
    D36k = cos.(kx) * Dlt3 + cos.(ky) * Dlt4 + d_sym * (Dlt5 + Dlt6)
    """ --- """
    I1 = (1 - p) * I; εpp = ε_p - mu
    """ ====== Normal elements of Matrix Dij and Kij at +k =================="""
    csx = cos.(kx .- 2 * α_x); csy = cos.(ky)
    snx = sin.(kx / 2 .- α_x); sny = sin.(ky / 2)
    psi = snx .* sny; sx2 = snx .^ 2; sy2 = sny .^ 2
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(2 * kx .- 4 * α_x) + cos.(2 * ky)) / 2
    Kk  = 3 / 4 .- C1 * gm1
    D11 = εpp .+ 2 * τ_ * sx2
    D12 = 2 * (τ_ - 2 * t) * psi
    D13 = J * Kk .* snx
    D23 = J * Kk .* sny
    D22 = εpp .+ 2 * τ_ * sy2
    D33 = (εpp + 4 * (τ_ - t)) * Kk .- 3 * (3 / 4 * τ_ - t) / 2 +
        + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
        + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    """ ====== Normal elements of Matrix Dij and Kij at -k =================="""
    csx = cos.(-kx .- 2 * α_x); csy = cos.(-ky)
    snx = sin.(-kx / 2 .- α_x); sny = sin.(-ky / 2)
    sx2 = snx .^ 2; sy2 = sny .^ 2; psi = snx .* sny
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(-2 * kx .- 4 * α_x) + cos(-2 * ky)) / 2
    Kk_  = 3 / 4 .- C1 * gm1
    D11_ = εpp .+ 2 * τ_ * sx2
    D12_ = 2 * (τ_ - 2 * t) * psi
    D13_ = J * Kk_ .* snx
    D23_ = J * Kk_ .* sny
    D22_ = εpp .+ 2 * τ_ * sy2
    D33_ = (εpp + 4 * (τ_ - t)) * Kk_ .- 9 / 8 * τ_ .+ 3 / 2 * t +
         + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
         + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    """ ====== Calculating spectrum in N-phase and SC-phase ================="""
    ESC = Array{Float64}(zeros(N, N, 6))
    E = aEk(kx, ky) .- mu; E1 = E[:, :, 1]; E2 = E[:, :, 2]; E3 = E[:, :, 3]
    E_=aEk_(kx, ky) .+ mu; E1_=E_[:, :, 3]; E2_=E_[:, :, 2]; E3_=E_[:, :, 1]
    Ek = sqrt.(((E1 - E1_) / 2) .^ 2 + D14k .^ 2 + D25k .^ 2 + D36k .^ 2 ./ Kk ./ Kk_)
    ESC[:, :, 1] = (E1 + E1_) / 2 + Ek; ESC[:, :, 2] = E2 ; ESC[:, :, 3] = E3
    ESC[:, :, 4] = (E1 + E1_) / 2 - Ek; ESC[:, :, 5] = E2_; ESC[:, :, 6] = E3_
    """ ====================================================================="""
    Fx = Array{Float64}(zeros(N, N))
    for ll in [1 4]
        om = ESC[:, :, ll]; DT3_ = (om - E1_) .* (om - E2_) .* (om - E3_)
        S1m = (om - D11);  S2m = (om - D22);  S3m = (om - D33 ./ Kk )
        S1p = (om + D11_); S2p = (om + D22_); S3p = (om + D33_./ Kk_)
        Q11 = S2m .* S3m - D23 .^ 2 ./ Kk
        Q22 = S1m .* S3m - D13 .^ 2 ./ Kk
        Q11_= S2p .* S3p - D23_ .^ 2 ./ Kk_
        Q22_= S1p .* S3p - D13_ .^ 2 ./ Kk_
        Q13_=-S2p .* D13_+ D12_ .* D23_
        Q33_= S1p .* S2p - D12_ .^ 2
        Q23_= S1p .* D23_- D12_ .* D13_
        """ ====== Denominator =============================================="""
        denom = 1
        for jj in 1:6
            if jj != ll
                denom = denom .* (om - ESC[:, :, jj])
            end
        end
        denom[findall(denom -> denom==0, denom)] .= 1e-10
        """=================================================================="""
        Qaa = Q11 .* DT3_
            - S3m .* Q22_ .* D25k .* D25k +
            + D23 .* Q23_ .* D25k .* D36k ./ Kk ./ Kk_ +
            + D23 .* Q23_ .* D25k .* D36k ./ Kk ./ Kk_ +
            - S2m .* Q33_ .* D36k .* D36k ./ Kk ./ Kk_ +
            + S1p .* D25k .* D25k .* D36k .* D36k ./ Kk ./ Kk_
        Qbb = Q22 .* DT3_
            - S3m .* Q11_ .* D14k .* D14k +
            - D13 .* Q13_ .* D14k .* D36k ./ Kk ./ Kk_ +
            - D13 .* Q13_ .* D14k .* D36k ./ Kk ./ Kk_ +
            - S1m .* Q33_ .* D36k .* D36k ./ Kk ./ Kk_ +
            + S2p .* D14k .* D14k .* D36k .* D36k ./ Kk ./ Kk_
    """ ====================================================================="""
    Fx = Fx .+ fermi(om, T) .* (Qaa + Qbb) ./ denom
    end
    return 2 * sum(Fx) / N ^ 2
end

function amuSC(vmu::Float64)
    global x, mu
    mu = vmu
     return anpSC(1) - x
end

function aEqB(Dlt::AbstractArray{T}) where T <: Float64
    global kx, ky, T, mu, Dlt1, Dlt2, Dlt3, Dlt4, Dlt5, Dlt6
    global ε_p, τ_, t, α_x, x, p, C1, C2, C3, J, I, Up, V2

    Dlt1 = Dlt[1]; Dlt2 = Dlt[2]; Dlt3 = Dlt[3]
    Dlt4 = Dlt[4]; Dlt5 = Dlt[5]; Dlt6 = Dlt[6]
    d_sym = (cos.(kx)-cos.(ky))
    D14k = d_sym * Dlt1
    D25k = d_sym * Dlt2
    D36k = cos.(kx) * Dlt3 + cos.(ky) * Dlt4 + d_sym * (Dlt5 + Dlt6)

    I1 = (1 - p) * I; εpp = ε_p - mu
    """ ====== Normal elements of Matrix Dij and Kij at +k =================="""
    csx = cos.(kx .- 2 * α_x); csy = cos.(ky)
    snx = sin.(kx / 2 .- α_x); sny = sin.(ky / 2)
    psi = snx .* sny; sx2 = snx .^ 2; sy2 = sny .^ 2
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(2 * kx .- 4 * α_x) + cos.(2 * ky)) / 2
    Kk  = 3 / 4 .- C1 * gm1
    D11 = εpp .+ 2 * τ_ * sx2
    D12 = 2 * (τ_ - 2 * t) * psi
    D13 = J * Kk .* snx
    D23 = J * Kk .* sny
    D22 = εpp .+ 2 * τ_ * sy2
    D33 = (εpp + 4 * (τ_ - t)) * Kk .- 3 * (3 / 4 * τ_ - t) / 2 +
        + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
        + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    """ ====== Normal elements of Matrix Dij and Kij at -k =================="""
    csx = cos.(-kx .- 2 * α_x); csy = cos.(-ky)
    snx = sin.(-kx / 2 .- α_x); sny = sin.(-ky / 2)
    sx2 = snx .^ 2; sy2 = sny .^ 2; psi = snx .* sny
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(-2 * kx .- 4 * α_x) + cos(-2 * ky)) / 2
    Kk_  = 3 / 4 .- C1 * gm1
    D11_ = εpp .+ 2 * τ_ * sx2
    D12_ = 2 * (τ_ - 2 * t) * psi
    D13_ = J * Kk_ .* snx
    D23_ = J * Kk_ .* sny
    D22_ = εpp .+ 2 * τ_ * sy2
    D33_ = (εpp + 4 * (τ_ - t)) * Kk_ .- 9 / 8 * τ_ .+ 3 / 2 * t +
         + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
         + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    """ ====== Calculating spectrum in N-phase and SC-phase ================="""
    ESC = Array{Float64}(zeros(N, N, 6))
    E = aEk(kx, ky) .- mu; E1 = E[:, :, 1]; E2 = E[:, :, 2]; E3 = E[:, :, 3]
    E_=aEk_(kx, ky) .+ mu; E1_=E_[:, :, 3]; E2_=E_[:, :, 2]; E3_=E_[:, :, 1]
    Ek = sqrt.(((E1 - E1_) / 2) .^ 2 + D14k .^ 2 + D25k .^ 2 + D36k .^ 2 ./ Kk ./ Kk_)
    ESC[:, :, 1] = (E1 + E1_) / 2 + Ek; ESC[:, :, 2] = E2 ; ESC[:, :, 3] = E3
    ESC[:, :, 4] = (E1 + E1_) / 2 - Ek; ESC[:, :, 5] = E2_; ESC[:, :, 6] = E3_
    """======================================================================"""
    Fx1 = Fx2 = Fx3 = Fx4 = Fx5 = Fx6 = Array{Float64}(zeros(N, N))
    px = sin.(+kx / 2 .- α_x); py = sin.(+ky / 2)
    mx = sin.(-kx / 2 .- α_x); my = sin.(-ky / 2)
    for ll in 1:6
        om = ESC[:, :, ll]
        S1m  = (om - D11 ); S2m = (om - D22 ); S3m = (om - D33  ./ Kk )
        S1p  = (om + D11_); S2p = (om + D22_); S3p = (om + D33_ ./ Kk_)
        Q11  =  S2m .* S3m - D23 .^ 2 ./ Kk
        Q12  = -D12 .* S3m - D13 .* D23 ./ Kk
        Q22  =  S1m .* S3m - D13 .^ 2 ./ Kk
        Q13  =  S2m .* D13 + D12 .* D23
        Q23  = -S1m .* D23 - D12 .* D13
        Q33  =  S1m .* S2m - D12 .^ 2
        Q11_ =  S2p .* S3p - D23_ .^ 2 ./ Kk_
        Q12_ = D12_ .* S3p - D13_ .* D23_ ./ Kk_
        Q22_ =  S1p .* S3p - D13_ .^ 2 ./ Kk_
        Q13_ = -S2p .* D13_ + D12_ .* D23_
        Q23_ =  S1p .* D23_ - D12_ .* D13_
        Q33_ =  S1p .* S2p - D12_ .^ 2
        """ ====== Denominator =============================================="""
        denom = 1
        for jj in 1:6
            if jj != ll
                denom = denom .* (om - ESC[:, :, jj])
            end
        end
        denom[findall(denom -> denom==0, denom)] .= 1e-10
        """=================================================================="""
        Qaa_= Q11 .* Q11_ .* D14k + Q12 .* Q12_ .* D25k +
            + Q13 .* Q13_ .* D36k ./ Kk ./ Kk_ +
            - S3m .* S3p .* D14k .* D25k .* D25k +
            - S2m .* S2p .* D14k .* D36k .* D36k ./ Kk ./ Kk_ +
            + D23 .* D23_ .* D14k .* D25k .* D36k ./ Kk ./ Kk_ +
            + D23 .* D23_ .* D14k .* D25k .* D36k ./ Kk ./ Kk_ +
            + D13 .* D13_ .* D25k .* D25k .* D36k ./ Kk ./ Kk_ +
            + D12 .* D12_ .* D25k .* D36k .* D36k ./ Kk ./ Kk_ +
            + D14k .* D25k .* D25k .* D36k .* D36k ./ Kk ./ Kk_
        Qab_= -Q11 .* Q12_ .* D14k - Q12 .* Q22_ .* D25k
            - Q13 .* Q23_ .* D36k ./ Kk ./ Kk_ +
            - D12 .* S1p .* D25k .* D36k .* D36k ./ Kk ./ Kk_ +
            - D23 .* D13_ .* D14k .* D25k .* D36k ./ Kk ./ Kk_ +
            + S2m .* D12_ .* D14k .* D36k .* D36k ./ Kk ./ Kk_
        Qba_= -Q12 .* Q11_ .* D14k - Q22 .* Q12_ .* D25k
            - Q23 .* Q13_ .* D36k ./ Kk ./ Kk_ +
            - D12 .* S2p .* D14k .* D36k .* D36k ./ Kk ./ Kk_ +
            + S1m .* D12_ .* D25k .* D36k .* D36k ./ Kk ./ Kk_+
            - D13 .* D23_ .* D14k .* D25k .* D36k ./ Kk ./ Kk_
        Qbb_= Q12 .* Q12_ .* D14k + Q22 .* Q22_ .* D25k +
            + Q23 .* Q23_ .* D36k ./ Kk ./ Kk_ +
            - S3m .* S3p .* D14k .* D14k .* D25k +
            - S1m .* S1p .* D25k .* D36k .* D36k ./ Kk ./ Kk_ +
            + D13 .* D13_ .* D14k .* D25k .* D36k ./ Kk ./ Kk_ +
            + D13 .* D13_ .* D14k .* D25k .* D36k ./ Kk ./ Kk_ +
            + D23 .* D23_ .* D14k .* D14k .* D36k ./ Kk ./ Kk_ +
            + D12 .* D12_ .* D14k .* D36k .* D36k ./ Kk ./ Kk_ +
            + D14k .* D14k .* D25k .* D36k .* D36k ./ Kk ./ Kk_
        QLL_= Q13 .* Q13_ .* D14k + Q23 .* Q23_ .* D25k + Q33 .* Q33_ .* D36k +
            - S2m .* S2p .* D14k .* D14k .* D36k +
            - S1m .* S1p .* D25k .* D25k .* D36k +
            + D12 .* D12_ .* D14k .* D25k .* D36k +
            + D12 .* D12_ .* D14k .* D25k .* D36k +
            + D23 .* D23_ .* D14k .* D14k .* D25k +
            + D13 .* D13_ .* D14k .* D25k .* D25k +
            + D14k .* D14k .* D25k .* D25k .* D36k
        Quu_= px .* mx .* Qaa_ + py .* my .* Qbb_ +
            + px .* my .* Qab_ + py .* mx .* Qba_
        """=================================================================="""
        Fermi = fermi(-om, T)
        Fx1 = Fx1 .+ Fermi .* (-2 * V2 * cos.(kx) .* Qaa_) ./ denom
        Fx2 = Fx2 .+ Fermi .* (-2 * V2 * cos.(kx) .* Qbb_) ./ denom
        Fx3 = Fx3 .+ Fermi .* (Up * C1 / 2 * cos(2*α_x) * Qaa_) ./ denom
        Fx4 = Fx4 .+ Fermi .* (Up * C1 / 2 * Qbb_) ./ denom
        Fx5 = Fx5 .+ Fermi .* (-V2 * C1 * cos.(kx) .* (Qaa_ + Qbb_)) ./ denom
        Fx6 = Fx6 .+ Fermi .* (I1 * d_sym .* (QLL_ - C1 * Quu_)) ./ denom
    end
    return [sum(Fx1) / N ^ 2 - Dlt1;
            sum(Fx2) / N ^ 2 - Dlt2;
            sum(Fx3) / N ^ 2 - Dlt3;
            sum(Fx4) / N ^ 2 - Dlt4;
            sum(Fx5) / N ^ 2 - Dlt5;
            sum(Fx6) / N ^ 2 - Dlt6];
end

function acurSC(float)
    global kx, ky, T, mu, Dlt1, Dlt2, Dlt3, Dlt4, Dlt5, Dlt6
    global ε_p, τ_, t, α_x, x, p, C1, C2, C3, J, I, Up, V2

    mu = fzero(amuSC, [-10, -0], xatol=1e-10)
    """ ====================================================================="""
    d_sym = (cos.(kx)-cos.(ky))
    D14k = d_sym * Dlt1
    D25k = d_sym * Dlt2
    D36k = cos.(kx) * Dlt3 + cos.(ky) * Dlt4 + d_sym * (Dlt5 + Dlt6)

    I1 = (1 - p) * I; εpp = ε_p - mu
    """ ====== Normal elements of Matrix Dij and Kij at +k =================="""
    csx = cos.(kx .- 2 * α_x); csy = cos.(ky)
    snx = sin.(kx / 2 .- α_x); sny = sin.(ky / 2)
    psi = snx .* sny; sx2 = snx .^ 2; sy2 = sny .^ 2
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(2 * kx .- 4 * α_x) + cos.(2 * ky)) / 2
    Kk  = 3 / 4 .- C1 * gm1
    D11 = εpp .+ 2 * τ_ * sx2
    D12 = 2 * (τ_ - 2 * t) * psi
    D13 = J * Kk .* snx
    D22 = εpp .+ 2 * τ_ * sy2
    D23 = J * Kk .* sny
    D33 = (εpp + 4 * (τ_ - t)) * Kk .- 3 * (3 / 4 * τ_ - t) / 2 +
        + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
        + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    """ ====== Normal elements of Matrix Dij and Kij at -k =================="""
    csx = cos.(-kx .- 2 * α_x); csy = cos.(-ky)
    snx = sin.(-kx / 2 .- α_x); sny = sin.(-ky / 2)
    sx2 = snx .^ 2; sy2 = sny .^ 2; psi = snx .* sny
    gm1 = (csx + csy) / 2
    gm2 = csx .* csy
    gm3 = (cos.(-2 * kx .- 4 * α_x) + cos(-2 * ky)) / 2
    Kk_  = 3 / 4 .- C1 * gm1
    D11_ = εpp .+ 2 * τ_ * sx2
    D12_ = 2 * (τ_ - 2 * t) * psi
    D13_ = J * Kk_ .* snx
    D22_ = εpp .+ 2 * τ_ * sy2
    D23_ = J * Kk_ .* sny
    D33_ = (εpp + 4 * (τ_ - t)) * Kk_ .- 9 / 8 * τ_ .+ 3 / 2 * t +
         + (τ_ - 2 * t) * C2 * gm2 + τ_ / 2 * C3 * gm3 .- 3 / 4 * J +
         + J * C1 * (1 / 4 .+ 2 * gm1) - I1 * C1 * (gm1 .+ 4)
    """ ====== Calculating spectrum in N-phase and SC-phase ================="""
    ESC = Array{Float64}(zeros(N, N, 6))
    E = aEk(kx, ky) .- mu; E1 = E[:, :, 1]; E2 = E[:, :, 2]; E3 = E[:, :, 3]
    E_=aEk_(kx, ky) .+ mu; E4= E_[:, :, 3]; E5 =E_[:, :, 2]; E6 =E_[:, :, 1]
    Ek = sqrt.(((E1 - E4) / 2) .^ 2 + D14k .^ 2 + D25k .^ 2 + D36k .^ 2 ./ Kk ./ Kk_)
    ESC[:, :, 1] = (E1 + E4) / 2 + Ek; ESC[:, :, 2] = E2; ESC[:, :, 3] = E3
    ESC[:, :, 4] = (E1 + E4) / 2 - Ek; ESC[:, :, 5] = E5; ESC[:, :, 6] = E6
    """ ====================================================================="""
    Fx = Array{Float64}(zeros(N, N))
    for ll in [1 4]
        om = ESC[:, :, ll]
        S1p = (om + D11_); DT3_= (om - E4) .* (om - E5) .* (om - E6)
        S2m = (om - D22 ); S3m = (om - D33 ./ Kk )
        S2p = (om + D22_); S3p = (om + D33_./ Kk_)
        Q11 = S2m .* S3m - D23 .^ 2 ./ Kk
        Q12 =-D12 .* S3m - D13 .* D23 ./ Kk
        Q13 = S2m .* D13 + D12 .* D23
    	Q12_= D12_.* S3p - D13_.* D23_./ Kk_
    	Q22_= S1p .* S3p - D13_.^ 2 ./ Kk_
    	Q13_=-S2p .* D13_+ D12_.* D23_
    	Q33_= S1p .* S2p - D12_.^ 2
    	Q23_= S1p .* D23_- D12_.* D13_
        """ ================================================================="""
	    Qaa = Q11 .* DT3_ +
            - S3m .* Q22_ .* D25k .* D25k +
            + D23 .* Q23_ .* D25k .* D36k ./ Kk ./ Kk_ +
            + D23 .* Q23_ .* D25k .* D36k ./ Kk ./ Kk_ +
            - S2m .* Q33_ .* D36k .* D36k ./ Kk ./ Kk_ +
            + S1p .* D25k .* D25k .* D36k .* D36k ./ Kk ./ Kk_
	    Qab =-Q12 .* DT3_ - S3m .* Q12_ .* D14k .* D25k +
            + D23 .* Q13_ .* D14k .* D36k ./ Kk ./ Kk_ +
            - D13 .* Q23_ .* D25k .* D36k ./ Kk ./ Kk_ +
            - D12 .* Q33_ .* D36k .* D36k ./ Kk ./ Kk_ +
            + D12_ .* D14k .* D25k .* D36k .* D36k ./ Kk ./ Kk_
	    QaL = Q13 .* DT3_ +
            - D23 .* Q12_ .* D14k .* D25k +
            - D13 .* Q22_ .* D25k .* D25k +
            + S2m .* Q13_ .* D14k .* D36k ./ Kk_ +
            - D12 .* Q23_ .* D25k .* D36k ./ Kk_ +
            + D13_ .* D14k .* D25k .* D25k .* D36k ./ Kk_
        """ ====== Denominator =============================================="""
        denom = 1
        for jj in 1:6
            if jj != ll
                denom = denom .* (om - ESC[:, :, jj])
            end
        end
        denom[findall(denom -> denom==0, denom)] .= 1e-10
        """=================================================================="""
        Faa = 4 * τ_ * sin.(kx / 2 .- α_x) .* Qaa
	    Fab = 4 * (τ_ - 2 * t) * sin.(ky / 2) .* Qab
	    FaL = 2 * J * QaL
	    JhQ = cos.(kx / 2 .- α_x) .* (Faa + Fab + FaL)
	    Fx = Fx .+ fermi(om, T) .* JhQ ./ denom
    end
    return sum(Fx) / N ^ 2
end

begin
    global kx, ky, T, mu, Dlt1, Dlt2, Dlt3, Dlt4, Dlt5, Dlt6
    global ε_p, τ_, t, α_x, x, p, C1, C2, C3, J, I, Up, V2

    ε_p = 0.0;   # энергия связи p - дырки (eV)
    τ_  = 0.1;   # effective hole hoppings (O-subsystem)
    J   = 3.4;   # 2.391695906432749"""
    I   = 0.136; # Cu-Cu exchange interaction (for LSCO I = 0.136 eV)
    t   = 0.11;  # integral of direct O-O hoppings (eV)
    Up  = 4.0;   # on-site Coulomb repulsion Up = [4 4.7 6]
    V2  = 0.12;  # Coulomb between the next-nearest-neighbor
    α_x = 2e-3;  # External field correction
    TK  = 2;     # T=8.61734e-5*TK; Temperature in Kelvin and conversion to eV
    x   = 0.17;  # x -- concentration of holes per cell
    N   = 10;   # number of discretization of the Brillouin zone
    """ ====== Tabl of vulues: x, p, C1, C2, C3 ============================="""
                   vX =[ 0.030,  0.070,  0.150,  0.220,  0.3000,  0.400]
    p_x = Spline1D(vX, [ 0.150,  0.210,  0.250,  0.275,  0.3000,  0.320])
    C1x = Spline1D(vX, [-0.287, -0.255, -0.231, -0.214, -0.1940, -0.170])
    C2x = Spline1D(vX, [ 0.124,  0.075,  0.036,  0.009, -0.0222, -0.022])
    C3x = Spline1D(vX, [ 0.095,  0.064,  0.051,  0.045,  0.0400,  0.020])
    """ ====== Discretization of the Brillouin zone ========================="""
    k0 = range(-π+1e-6, π-1e-6, length = N); (kx, ky) = ndgrid(k0) # range(a, b, step=0.01)

    T  = 8.61734e-5 * TK; p = p_x(x); C1 = C1x(x); C2 = C2x(x); C3 = C3x(x);
end

@benchmark acurSC(1)

    # E = aEk(kx, ky); surface(kx, ky, E[:, :, 1] .- mu)
    # fzero(amuN, [-8, 0], xatol=1e-10)

    # @show mu = find_zero(amuN, (-8, 0), Bisection())
    # mu1 = find_zero(amuN, (-8, 0), Roots.Brent())
    # vTK = [2 5 10 15 20 25 30]
    # Cur = Array{Float64}(zeros(length(vTK), 2))
    # for j in 1:length(vTK)
    #     TK = vTK[j];  T = TK * 8.61734e-5
    #     Cur[j, :] = [TK, acurSC(1)]
    # end
    # plot(Cur[:, 1], -87.42 * 2 * Cur[:, 2] / α_x)



# save([path 'Cur.T..x_',num2str(x),'.J_',num2str(J),'.I_',num2str(I),'.t_',num2str(t),'.tau_',num2str(taum),...
# '.Up_',num2str(Up),'.V2_',num2str(V2),'.Nk_',num2str(N),'.T_',num2str(TK),'.mat'],'Cur','-ascii','-double');
# Cur=load([path 'Cur.T..x_0.17.J_3.40.I_0.136.t_0.11.tau_0.100.Up_4.0.V2_0.12.Nk_1100_2.0.mat'],'-ascii');
