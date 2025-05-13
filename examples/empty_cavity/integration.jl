using Integrals
using GLMakie

function wavelet(t,p)
    return exp(-(t-p[1])^2/p[2]^2)*cos(p[3]*(t-p[1]))
end

function wavelet_squared(t,p)
    return (wavelet(t,p))^2
end
p = [21,4,8]

prob = IntegralProblem(wavelet_squared,(-Inf,Inf),p)
sol = solve(prob,QuadGKJL())
println(sol.u)
#=
t = LinRange(-10,10,1000)
vals = [wavelet(time,p) for time in t]

lines(t,vals)
=#