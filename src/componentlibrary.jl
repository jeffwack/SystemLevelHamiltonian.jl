function cavity(name)
    #Define the filter cavity 
    hilb = FockSpace(name)
    a = Destroy(hilb,Symbol(:a,:_,name))

    (κ, Δ) = rnumbers(Symbol(:κ,:_,name),Symbol(:Δ,:_,name))

    return SLH(name,
                [:Input],
                [:Output],
                [1],
                [sqrt(κ)*a],
                Δ*adjoint(a)*a)
end

function squeezing_cavity(name)
    hilb = FockSpace(:squeezer)
    a = Destroy(hilb, Symbol(:a,:_,name))

    (κ,ϵ) = rnumbers(Symbol(:κ,:_,name),Symbol(:ϵ,:_,name))

    return SLH(name,
                ["Input"],
                ["Output"],
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
                [:Input],
                [:Output],
                [1],
                [sqrt(κ)*a],
                Δ*a'*a+Ω*b'*b - g*a'*a*(b'+b))
end
