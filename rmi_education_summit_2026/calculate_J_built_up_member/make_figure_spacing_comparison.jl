# make_figure_spacing_comparison.jl — Makie plot: Ferrite J_eff vs weld spacing and the Tlumak equation (3 in welds).
using CairoMakie, DelimitedFiles
blue = "#2a78d6"; orange = "#eb6834"; aqua = "#1baf7a"; violet = "#4a3aa7"; ink = "#0b0b0b"; ink2 = "#52514e"; grid = "#e6e5e1"

res = readdlm(joinpath(@__DIR__, "weld_spacing_results.csv"), ','; skipstart = 1)
ser = String.(strip.(res[:, 1])); L = Float64.(res[:, 2]); wl = Float64.(res[:, 3]); ws = Float64.(res[:, 5])
Jf = Float64.(res[:, 7]); Jtube = Float64(res[1, 12]); JC = Float64(res[1, 13])
equation(w, s) = Jtube * w / s + JC + JC * (1 - w / s)

fig = Figure(size = (1500, 640), fontsize = 15, backgroundcolor = :white)
xt = [0, 18, 36, 54, 72, 90, 108]

ax1 = Axis(fig[1, 1]; xlabel = rich("weld spacing, w", subscript("spacing"), " (in)"), ylabel = rich("J", subscript("eff"), " (in⁴)"),
           title = "Ferrite shell model vs. Tlumak equation, 3 in welds", titlealign = :left, titlesize = 15,
           xgridcolor = grid, ygridcolor = grid, xticks = (xt, string.(xt)))
ss = collect(range(3.0, 110.0, length = 300))
hlines!(ax1, [Jtube]; color = ink2, linewidth = 1, linestyle = :dot)
text!(ax1, 40.0, Jtube - 0.04; text = rich("J", subscript("tube"), " (Bredt closed cell) = $(round(Jtube, digits = 3)) in⁴"), fontsize = 12, color = ink2, align = (:left, :top))
hlines!(ax1, [2JC]; color = ink2, linewidth = 1, linestyle = :dot)
text!(ax1, 70.0, 2JC + 0.12; text = rich("2 J", subscript("C"), " (no welds) = $(round(2JC, sigdigits = 3)) in⁴"), fontsize = 12, color = ink2, align = (:left, :bottom))
# the w_spacing = 3 in point (37 abutting 3 in welds) is the fully welded member
for (name, w, col, mk) in (("A", 3.0, blue, :circle),)
    m = ser .== name; idx = sortperm(ws[m])
    xs_ = ws[m][idx]; ys_ = Jf[m][idx]
    lines!(ax1, ss[ss .>= w], equation.(w, ss[ss .>= w]); color = col, linewidth = 1.5, linestyle = :dash, label = "Tlumak equation, $(Int(w)) in welds")
    scatterlines!(ax1, xs_, ys_; color = col, marker = mk, markersize = 12, linewidth = 2, label = rich("Ferrite shell model, $(Int(w)) in welds (w", subscript("spacing"), " = 3 in: fully welded)"))
end
axislegend(ax1; position = :rt, framevisible = false, labelsize = 11)

ax2 = Axis(fig[1, 2]; xlabel = rich("weld spacing, w", subscript("spacing"), " (in)"), ylabel = rich("J", subscript("eff"), " / J", subscript("tube")),
           title = rich("Normalized by the Bredt closed-cell J", subscript("tube")), titlealign = :left, titlesize = 15,
           xgridcolor = grid, ygridcolor = grid, xticks = (xt, string.(xt)), yticks = 0:0.2:1.2)
hlines!(ax2, [1.0]; color = ink2, linewidth = 1, linestyle = :dot)
for (name, w, col, mk) in (("A", 3.0, blue, :circle),)
    m = ser .== name; idx = sortperm(ws[m])
    xs_ = ws[m][idx]; ys_ = Jf[m][idx]
    lines!(ax2, ss[ss .>= w], equation.(w, ss[ss .>= w]) ./ Jtube; color = col, linewidth = 1.5, linestyle = :dash, label = "Tlumak equation")
    scatterlines!(ax2, xs_, ys_ ./ Jtube; color = col, marker = mk, markersize = 12, linewidth = 2, label = "Ferrite shell model, $(Int(w)) in welds")
end
axislegend(ax2; position = :rt, framevisible = false, labelsize = 12)
ylims!(ax2, 0, 1.25); xlims!(ax1, 0, 112); xlims!(ax2, 0, 112); ylims!(ax1, 0, 2.3)

Label(fig[0, :], rich("Effective J of the two-C welded upright vs. weld spacing — Ferrite.jl shell model and Tlumak J", subscript("combined"), " equation"),
      fontsize = 17, font = :bold, halign = :left, color = ink)
save(joinpath(@__DIR__, "weld_spacing_comparison.png"), fig; px_per_unit = 2)
println("wrote weld_spacing_comparison.png")
