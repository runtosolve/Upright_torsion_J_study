# plot_mode_wglmakie.jl — 3D buckling mode shape of the two-C welded upright from mode_*.csv written by
# buckling_built_up.jl. Interactive WGLMakie scene exported to a standalone HTML (Bonito), plus a CairoMakie PNG.
# Run:  julia --project=. plot_mode_wglmakie.jl [mode_global_1.csv]
using DelimitedFiles, Statistics, LinearAlgebra
using WGLMakie, Bonito
import CairoMakie
using WGLMakie.Makie.GeometryBasics: Point3f, GLTriangleFace, Mesh

include(joinpath(@__DIR__, "plot_mode_functions.jl"))

path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "mode_global_1.csv")
res = readdlm(joinpath(@__DIR__, "buckling_results.csv"), ','; skipstart = 1)
P = Float64(res[1, 4])
title = rich("Global flexural-torsional buckling mode, two-C welded upright, L = 44 in, pinned warping-free, 3 in welds at 18 in:  P", subscript("cre"), " = $(round(P, digits = 1)) kips")
base = splitext(basename(path))[1]

# interactive WGLMakie scene → standalone HTML
WGLMakie.activate!()
app = App() do
    mode_figure(path, title)
end
# export to a local temp file first (file close on the Google Drive path can time out), then copy over
html = joinpath(@__DIR__, base * "_wglmakie.html")
tmp = joinpath(mktempdir(), base * "_wglmakie.html")
Bonito.export_static(tmp, app)
cp(tmp, html; force = true)
println("wrote ", html, "  (", round(filesize(html) / 1e6, digits = 1), " MB)")

# static PNG with the same figure code
CairoMakie.activate!()
fig = mode_figure(path, title)
png = joinpath(@__DIR__, base * ".png")
CairoMakie.save(png, fig; px_per_unit = 2)
println("wrote ", png)
