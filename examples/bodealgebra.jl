#This is supposed to render a bode diagram algerbra derivation of 
#a closed loop transfer function as was done in group meeting

using ControlSystems
using GLMakie

G = tf([1], [1, 1])      # Example: low-pass filter
H = tf([10], [1, 10])    # Example: high-pass filter

T = feedback(G * H)      # Closed-loop transfer function

function bode_image(sys, filename)
    p = bodeplot(sys)
    savefig(p, filename)
end

bode_image(G, "G.png")
bode_image(H, "H.png")
bode_image(T, "T.png")

function render_diagram()
    fig = Figure()
    ax = Axis(fig[1, 1], aspect = DataAspect())

    # Insert Bode diagram images as boxes
    image!(ax, load("G.png"), position = Point2f(0, 0), scale = 0.3)
    image!(ax, load("H.png"), position = Point2f(2, 0), scale = 0.3)
    image!(ax, load("T.png"), position = Point2f(4, 0), scale = 0.3)

    # Draw arrows and labels manually as needed
    lines!(ax, [0.5, 1.5], [0, 0], color = :black, linewidth = 2)
    # etc.

    fig
end
