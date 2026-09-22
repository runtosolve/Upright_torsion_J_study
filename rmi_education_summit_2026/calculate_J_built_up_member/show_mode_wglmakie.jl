# show_mode_wglmakie.jl — open the global FT buckling mode interactively in the browser with WGLMakie.
# Run:  julia --project=. show_mode_wglmakie.jl [mode_global_1.csv]     (Ctrl-C to stop the server)
using DelimitedFiles, Statistics
using WGLMakie, Bonito
using WGLMakie.Makie.GeometryBasics: Point3f, GLTriangleFace, Mesh
include(joinpath(@__DIR__, "plot_mode_functions.jl"))

path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "mode_global_1.csv")
res = readdlm(joinpath(@__DIR__, "buckling_results.csv"), ','; skipstart = 1)
P = Float64(res[1, 4])
title = rich("Global flexural-torsional buckling mode, two-C welded upright, L = 44 in, pinned warping-free, 3 in welds at 18 in:  P", subscript("cre"), " = $(round(P, digits = 1)) kips")

WGLMakie.activate!()
Bonito.browser_display()
fig = mode_figure(path, title)
display(fig)
println("WGLMakie scene served to the browser; drag to rotate, scroll to zoom. Press Ctrl-C to stop.")
while true
    sleep(1)
end
