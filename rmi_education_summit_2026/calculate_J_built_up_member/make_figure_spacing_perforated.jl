# make_figure_spacing_perforated.jl — J_eff vs weld spacing (3 in welds): Ferrite shell model with and without the
# perforation pattern (Gmsh mesh), and the Tlumak equation.   Run:  julia --project=. make_figure_spacing_perforated.jl
using CairoMakie, DelimitedFiles
blue = "#2a78d6"; orange = "#eb6834"; ink = "#0b0b0b"; ink2 = "#52514e"; grid = "#e6e5e1"
d = readdlm(joinpath(@__DIR__, "perforated_weld_spacing_results.csv"), ','; skipstart = 1)
ws = Float64.(d[:, 2]); J0 = Float64.(d[:, 3]); J1 = Float64.(d[:, 4])
idx = sortperm(ws); ws = ws[idx]; J0 = J0[idx]; J1 = J1[idx]
r = readdlm(joinpath(@__DIR__, "weld_spacing_results.csv"), ','; skipstart = 1)
Jtube = Float64(r[1, 12]); JC = Float64(r[1, 13])
p = readdlm(joinpath(@__DIR__, "perforated_J_gmsh_results.csv"), ','; skipstart = 1)
JC_perf = Float64(p[findfirst(==("single C, perforated"), String.(strip.(p[:, 1]))), 2])
equation(w, s) = Jtube * w / s + JC + JC * (1 - w / s)
ss = collect(range(3.0, 110.0, length = 300)); w = 3.0

fig = Figure(size = (1100, 820), fontsize = 16, backgroundcolor = :white)
ax = Axis(fig[1, 1]; xlabel = rich("weld spacing, w", subscript("spacing"), " (in)"), ylabel = rich("J", subscript("eff"), " (in⁴)"),
          title = "Ferrite shell model vs. Tlumak equation, 3 in welds: gross and perforated uprights", titlealign = :left, titlesize = 16,
          xgridcolor = grid, ygridcolor = grid, xticks = [0, 18, 36, 54, 72, 90, 108])
lines!(ax, ss, equation.(w, ss); color = blue, linewidth = 1.5, linestyle = :dash, label = "Tlumak equation, 3 in welds (gross section)")
scatterlines!(ax, ws, J0; color = blue, marker = :circle, markersize = 12, linewidth = 2, label = rich("Ferrite shell model, no holes (w", subscript("spacing"), " = 3 in: fully welded)"))
scatterlines!(ax, ws, J1; color = orange, marker = :rect, markersize = 12, linewidth = 2, label = "Ferrite shell model, perforated (teardrop web holes, square flange holes)")
hlines!(ax, [Jtube]; color = ink2, linewidth = 1, linestyle = :dot)
text!(ax, 40.0, Jtube - 0.04; text = rich("J", subscript("tube"), " (Bredt closed cell) = $(round(Jtube, digits = 3)) in⁴"), fontsize = 12, color = ink2, align = (:left, :top))
hlines!(ax, [2JC]; color = ink2, linewidth = 1, linestyle = :dot)
text!(ax, 70.0, 2JC + 0.12; text = rich("2 J", subscript("C"), " (no welds) = $(round(2JC, sigdigits = 3)) in⁴ gross, $(round(2JC_perf, sigdigits = 3)) perforated"), fontsize = 12, color = ink2, align = (:left, :bottom))
for (x, a, b) in zip(ws, J0, J1)
    x >= 6 && text!(ax, x + 1.0, (a + b) / 2; text = "−$(round(Int, 100 * (1 - b / a)))%", fontsize = 11, color = orange, align = (:left, :center))
end
axislegend(ax; position = :rt, framevisible = false, labelsize = 12)
xlims!(ax, 0, 112); ylims!(ax, 0, 2.3)
save(joinpath(@__DIR__, "weld_spacing_comparison_perforated.png"), fig; px_per_unit = 2)
println("wrote weld_spacing_comparison_perforated.png")
