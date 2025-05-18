using QuantumOptics
using GLMakie

function wignerfock!(ax,coeff)
    b = FockBasis(50)
    psi = normalize!(sum([a*fockstate(b, ii-1) for (ii,a) in enumerate(coeff)]))
    x = [-5:0.1:5;]
    W = wigner(psi, x, x)
    heatmap!(ax,W)
end

function fockax(ax, max_n)

    #hidedecorations!(ax)
    #hidespines!(ax)
    deactivate_interaction!(ax, :rectanglezoom)
    deactivate_interaction!(ax, :dragpan)
    #deactivate_interaction!(ax, :scrollzoom)

    alpha = 2

    sliders = [(label = string(n), range = 0:0.1:100, startvalue = alpha^n/sqrt(factorial(n))) for n in 0:max_n]

    sg = SliderGrid(fig[2,1],
                    sliders...)
    for changedslider in sg.sliders
        on(changedslider.value) do update

            coeffs = [sl.value[] for sl in sg.sliders]
            wignerfock!(ax,coeffs)
        end
    end
    coeffs = [sl.value[] for sl in sg.sliders]
    wignerfock!(ax,coeffs)
end


function newfockax(ax)

    #hidedecorations!(ax)
    #hidespines!(ax)
    deactivate_interaction!(ax, :rectanglezoom)
    deactivate_interaction!(ax, :dragpan)
    #deactivate_interaction!(ax, :scrollzoom)

    alpha = 2
    max_n = 4

    sliders = []
    for ii in 1:max_n
        push!(sliders, Slider2(fig[1,ii+1],xrange = -1:0.1:1, yrange = -1:0.1:1, startvalue = (alpha^ii/sqrt(factorial(ii)),0)))
    end

    for changedslider in sliders
        on(changedslider.value) do update
            coeffs = [sl.value[][1]+1.0im*sl.value[][2] for sl in sliders]
            wignerfock!(ax,coeffs)
        end
    end
    coeffs = [sl.value[][1]+1.0im*sl.value[][2] for sl in sliders]
    println(coeffs)
    wignerfock!(ax,coeffs)
end


fig = Figure()
ax = Axis(fig[1,1],aspect = 1)
sg = newfockax(ax)
rowsize!(fig.layout,  1, Auto(3))
fig