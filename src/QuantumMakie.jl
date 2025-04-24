
function plotlist(sol,idxs,names=idxs)

    curves = [[real(x[idx]) for x in sol.u] for idx in idxs]

    fig = Figure()
    ax = Axis(fig[1,1])
    for (name, curve) in zip(names,curves)
        lines!(ax,sol.t,curve,label = name)
    end 
    axislegend(ax)
    return fig
    
end

function plotlist(sol,ops)

    strings = []

    #this loop is for removing parentheses which are only sometimes added around the operator names
    #they have to be removed because they never appear in sys.unknowns
    for op in ops
        S = string(op)
        if S[1] == '('
            S = S[2:end-1]
        end
        push!(strings,S)
    end

    #adds the expected value brackets
    namestoplot = ["⟨"*name*"⟩" for name in strings]
    
    #this creates a list of strings sharing the same position as the variables
    unknown_names = [match(r"""(["])(.*?)\1""",string(n))[2] for n in sol.prob.f.sys.unknowns]

    #gets the indices of the variables we want to plot by matching up the strings
    idxs = [findfirst(x->x==name , unknown_names) for name in namestoplot]

    return plotlist(sol,idxs,strings)
end

function plotall(sol)
    idxs = eachindex(sol.prob.f.sys.unknowns)
    names = [match(r"""(["])(.*?)\1""",string(n))[2] for n in sol.prob.f.sys.unknowns]
    return plotlist(sol,idxs,names)
end