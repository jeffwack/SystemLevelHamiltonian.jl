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


function compute_qfi_old(ρ::AbstractMatrix, ρ_dot::AbstractMatrix; eps=1e-5)
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

function compute_qfi(sys, Ncutoff, T,params, param,backend) 
    symb = collect(keys(params))
    p0 = [params[sy] for sy in symb]
    idx = findfirst(isequal(param),symb)

    function closure(para)
        p = copy(p0)
        p[idx] = para
        return simulate_density_matrix(sys, Ncutoff, T,symb, p)
    end
    ρ = closure(p0[idx])
    ρ_dot = derivative(closure,backend,p0[idx])

    L_op = sld_operator(ρ, ρ_dot, eps=1e-4)
    return tr(ρ * L_op * L_op)
end


function simulate_density_matrix(sys,Ncutoff,T,params,p)
    paramrules = Dict([ps => pn for (ps,pn) in zip(params,p) ])
    (numsys,Ψ0) = convert_to_QT(sys,Ncutoff,paramrules)
    #=
    println(numsys)
    println(Ψ0)
    =#
    sol = mesolve(numsys.H, Ψ0, T, numsys.L)
    return hermitian_data(sol.states[end])
end

"""
hermitian_data(rho)

Where rho is a density matrix which can be accessed as rho.data. This function returns the average of rho and the adjoint of rho to guarantee Hermitian output.
"""
function hermitian_data(rho)
    return 1/2*(rho.data + adjoint(rho.data))    
end