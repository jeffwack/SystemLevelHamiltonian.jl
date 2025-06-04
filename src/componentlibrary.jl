function cavity(name)
    #Define the filter cavity 
    hilb = FockSpace(name)
    a = Destroy(hilb,Symbol(:a,:_,name))

    (κ, Δ) = rnumbers(Symbol(:κ,:_,name),Symbol(:Δ,:_,name))

    return SLH(name,
                [:In],
                [:Out],
                [1],
                [sqrt(κ)*a],
                Δ*adjoint(a)*a)
end

function qed_cavity(name)
    hcav = FockSpace(Symbol(:cav,:_,name))
    hspin = NLevelSpace(Symbol(:cav,:_,name),2)
    hilb = QuantumCumulants.tensor(hcav,hspin)

    # Define the operators
    a = Destroy(hilb,Symbol(:a,:_,name))

    σ(i,j) = Transition(hilb, Symbol(:σ,:_,name), i, j, 2)

    (Δ, g, κ, h,) = cnumbers(Symbol(:Δ,:_,name),Symbol(:g,:_,name),Symbol(:κ,:_,name),Symbol(:h,:_,name))

    # Hamiltonian
    H = Δ*a'*a + g*(a'*σ(1,2) + a*σ(2,1)) + h*(a + a')
    L = [sqrt(κ)*a]

    return SLH(name,
                [:In],
                [:Out],
                [1],
                L,
                H)
end


function squeezing_cavity(name)
    hilb = FockSpace(:squeezer)
    a = Destroy(hilb, Symbol(:a,:_,name))

    (κ,ϵ) = rnumbers(Symbol(:κ,:_,name),Symbol(:ϵ,:_,name))

    return SLH(name,
                [:In],
                [:Out],
                [1],
                [sqrt(κ)*a],
                1im*ϵ*(adjoint(a)^2- a^2))
end

function radiation_pressure_cavity(name)
    hilb = QuantumCumulants.tensor(FockSpace(Symbol(name,:_,:optical)), FockSpace(Symbol(name,:_,:mechanical)))
    a = Destroy(hilb,Symbol(:a,:_,name),1) 
    b = Destroy(hilb,Symbol(:b,:_,name),2) 

    (κ, Δ, Ω, g) = rnumbers(Symbol(:κ,:_,name),Symbol(:Δ,:_,name),Symbol(:Ω,:_,name),Symbol(:g,:_,name))

    return SLH(name,
                [:In],
                [:Out],
                [1],
                [sqrt(κ)*a],
                Δ*a'*a+Ω*b'*b - g*a'*a*(b'+b))
end
