using GLMakie
using QuantumToolbox

#scene = mesh(Sphere(Point3f(0), 1f0))
#display(scene)

fig = plot_wigner(basis(5, 0))


fig,ax = plot_fock_distribution(basis(5,0))