# make_figure_section_centers.jl — 2D cross section of the two-C welded upright with the composite centroid and the
# shear center found by statics (ft_analytical_check.jl).  Run:  julia --project=. make_figure_section_centers.jl
using CairoMakie, DelimitedFiles
include(joinpath(@__DIR__, "..", "calculate_J_single_member", "calculate_J_single_member.jl"))
X1, Y1 = centerline(shape); sp = section_properties(X1, Y1, shape.t)
xc = sp.xc + shape.B / 2; yc = sp.yc                                    # composite centroid (two equal C's, C2 offset by B)
res = readdlm(joinpath(@__DIR__, "ft_analytical_check.csv"), ','; skipstart = 1)
x_o = Float64(res[findfirst(==("x_o"), String.(res[:, 1])), 2])
xs = xc + x_o; ys = yc                                                  # shear center by statics (on the symmetry axis)
gray = "#8a8a8a"; ink = "#0b0b0b"; blue = "#2a78d6"; orange = "#eb6834"; grid = "#e6e5e1"

fig = Figure(size = (1100, 720), fontsize = 16, backgroundcolor = :white)
ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "X (in)", ylabel = "Y (in)", xgridcolor = grid, ygridcolor = grid,
          xticks = -3:1:6, yticks = 0:1:3, title = "Two-C welded upright: centroid and shear center from the Ferrite.jl shell model with 3 in welds at 18 in", titlealign = :left)
lines!(ax, X1, Y1; color = gray, linewidth = 5)
lines!(ax, X1 .+ shape.B, Y1; color = gray, linewidth = 5)
scatter!(ax, [xc], [yc]; color = blue, marker = :circle, markersize = 22, strokecolor = :white, strokewidth = 2)
scatter!(ax, [xs], [ys]; color = orange, marker = :diamond, markersize = 24, strokecolor = :white, strokewidth = 2)
text!(ax, xc, yc + 0.28; text = "centroid\n(x = $(round(xc, digits = 2)), y = $(round(yc, digits = 2)))", align = (:center, :bottom), fontsize = 14, color = ink)
text!(ax, xs + 0.9, ys - 0.28; text = "shear center, Ferrite shell model with welds\n(zero twist under a mid-length transverse force)\n(x = $(round(xs, digits = 2)), y = $(round(ys, digits = 2)))", align = (:center, :top), fontsize = 14, color = ink)
text!(ax, 3.02, -0.55; text = "C1 (web at X = 0, lips welded to the C2 web)          C2 (offset by B = 3 in)", align = (:center, :top), fontsize = 13, color = "#52514e")
xlims!(ax, -1.0, 6.6); ylims!(ax, -1.2, 3.9)
save(joinpath(@__DIR__, "section_centers.png"), fig; px_per_unit = 2)
println("wrote section_centers.png   centroid ($(round(xc, digits = 3)), $(round(yc, digits = 3)))  shear center ($(round(xs, digits = 3)), $(round(ys, digits = 3)))")
