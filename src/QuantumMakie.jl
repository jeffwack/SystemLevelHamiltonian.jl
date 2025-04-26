function paramwidget(sys,ops)
    fig = Figure()
    plotax = Axis(fig[1,1])

    ps = parameters(sys)

    sliders = [(label = string(symb), range = 0:0.1:10, startvalue = 1) for symb in ps]

    sg = SliderGrid(fig[2,1],
                    sliders...)
    
    u0 = zeros(ComplexF64,length(unknowns(sys)))

    for changedslider in sg.sliders
        on(changedslider.value) do update
            empty!(plotax)
            p0 = [sl.value[] for sl in sg.sliders]
            prob = ODEProblem(sys,u0,(0.0,100),ps.=>p0)
            sol = solve(prob,RK4())
            plotops!(plotax,sol,ops)
        end
    end   

    return fig
end