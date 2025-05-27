"""
    sld_operator(rho::AbstractMatrix, drho::AbstractMatrix; tol=1e-12)

Compute the symmetric logarithmic derivative (SLD) `L` for a given
density matrix `rho` and its parameter derivative `drho`, using the
eigenbasis method. Returns `L` in the original basis.

# Arguments
- `rho`: Hermitian density matrix (N×N)
- `drho`: Derivative of the density matrix with respect to some parameter
- `tol`: Threshold for eigenvalues considered nonzero (default: 1e-5)

# Returns
- `L`: Symmetric logarithmic derivative (N×N Hermitian matrix)
"""
function sld_operator(rho, drho; eps=1e-5)
    @assert size(rho) == size(drho) "rho and drho must be the same size"
    @assert issymmetric(rho) || ishermitian(rho) "rho must be Hermitian"
    @assert ishermitian(drho) "drho must be Hermitian"

    vals, vecs = eigen(Hermitian(rho))
    dim = size(rho, 1)
    L_eig = zeros(ComplexF64, dim, dim)

    for i in 1:dim
        for j in 1:dim
            pi, pj = vals[i], vals[j]
            denom = pi + pj
            if abs(denom) > eps
                vi = vecs[:, i]
                vj = vecs[:, j]
                L_eig[i, j] = 2 * (vi' * drho * vj)[1] / denom
            end
        end
    end

    # Reconstruct L in the original basis
    L = vecs * L_eig * vecs'
    return Hermitian(L)
end


function compute_qfi(ρ::AbstractMatrix, ρ_dot::AbstractMatrix; eps=1e-5)
    # Check Hermiticity
    @assert ishermitian(ρ) "ρ is not Hermitian"
    @assert ishermitian(ρ_dot) "ρ_dot is not Hermitian"

    evals, evecs = eigen(Hermitian(ρ))  # eigen-decomposition of ρ
    FQ = zero(eltype(ρ))

    for j in eachindex(evals)
        λ_j = evals[j]
        ψ_j = evecs[:, j]
        for k in eachindex(evals)
            λ_k = evals[k]
            ψ_k = evecs[:, k]
            denom = λ_j + λ_k
            if abs(denom) > eps && j ≠ k
                matrix_element = dot(ψ_k', ρ_dot* ψ_j)
                FQ += 4 * λ_j * abs2(matrix_element) / denom^2
            end
        end
    end

    return FQ 
end

function compute_qfi_alt(ρ::AbstractMatrix, ρ_dot::AbstractMatrix; eps=1e-5)
    # Check Hermiticity
    @assert ishermitian(ρ) "ρ is not Hermitian"
    @assert ishermitian(ρ_dot) "ρ_dot is not Hermitian"

    evals, evecs = eigen(Hermitian(ρ))  # eigen-decomposition of ρ
    FQ = zero(eltype(ρ))

    for j in eachindex(evals)
        λ_j = evals[j]
        ψ_j = evecs[:, j]
        for k in eachindex(evals)
            λ_k = evals[k]
            ψ_k = evecs[:, k]
            denom = λ_j + λ_k

            if abs(denom) > eps && j ≠ k
                matrix_element = dot(ψ_k', ρ_dot* ψ_j)
                FQ += 2 * abs2(matrix_element) / denom
            end
        end
    end

    return FQ 
end

function compute_qfi_fdm(H, L, ps, Ncutoff, Ψ₀, T, T_f, p, param_index, dh; eps=1e-12)
    p1 = Base.setindex(p, p[param_index] - dh/2, param_index)
    p2 = Base.setindex(p, p[param_index] + dh/2, param_index)

    ρ1 = simulate_density_matrix(H, L, ps, Ncutoff, Ψ₀, T, T_f, p1)
    ρ2 = simulate_density_matrix(H, L, ps, Ncutoff, Ψ₀, T, T_f, p2)
    ρ_dot = (ρ2 - ρ1) / dh
    ρ = (ρ1 + ρ2) / 2

    L_op = sld_operator(ρ, ρ_dot, eps=eps)
    return tr(ρ * L_op * L_op)
end

function simulate_density_matrix(H, L, ps, Ncutoff, Ψ₀, T, T_f, p)
    QTH = convert_to_QT([H], Dict(ps .=> p), Ncutoff)[1]
    QTL = convert_to_QT(L, Dict(ps .=> p), Ncutoff)
    sol = mesolve(QTH, Ψ₀, T, QTL)
    return data(sol.states[T_f])
end