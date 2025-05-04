using QuantumOptics
using GLMakie

function wignerfock!(ax,coeff_obs)
    b = FockBasis(50)
    psi = normalize!(sum([a[]*fockstate(b, ii) for (ii,a) in enumerate(coeff_obs)]))
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

            #calculate new normalization
            for slider in sg.sliders
            end

            coeffs = [sl.value[] for sl in sg.sliders]
            wignerfock!(ax,coeffs)
        end
    end
    coeffs = [sl.value[] for sl in sg.sliders]
    wignerfock!(ax,coeffs)
end

fig = Figure()
ax = Axis(fig[1,1],aspect = 1)
sg = fockax(ax,10)
rowsize!(fig.layout,  1, Auto(3))
fig