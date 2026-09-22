# make_figure_weld_spacing.jl — Ferrite J_eff vs. the proposed prediction equation as weld spacing varies
using CairoMakie, DelimitedFiles
blue = "#2a78d6"; orange = "#eb6834"; aqua = "#1baf7a"; violet = "#4a3aa7"; ink = "#0b0b0b"; ink2 = "#52514e"; grid = "#e6e5e1"

res = readdlm(joinpath(@__DIR__, "weld_spacing_results.csv"), ','; skipstart = 1)
ser = String.(strip.(res[:, 1])); L = Float64.(res[:, 2]); wl = Float64.(res[:, 3]); ws = Float64.(res[:, 5])
ratio = Float64.(res[:, 6]); Jf = Float64.(res[:, 7]); Jeq = Float64.(res[:, 8]); Jtube = Float64(res[1, 12]); JC = Float64(res[1, 13])
eq(r) = Jtube * r + JC + JC * (1 - r)

fig = Figure(size = (1500, 620), fontsize = 15, backgroundcolor = :white)
ax1 = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "w_length / w_spacing", ylabel = "J_eff (in⁴)",
           title = "Effective J vs. weld length / spacing (L = 111 in)", titlealign = :left, titlesize = 15,
           xgridcolor = grid, ygridcolor = grid, xticks = ([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0], ["0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1"]))
rr = 10 .^ range(log10(0.012), 0, length = 200)
lines!(ax1, rr, eq.(rr); color = ink2, linewidth = 2, linestyle = :dash, label = "equation: J_tube·(w/s) + J_back + J_front·(1 − w/s)")
hlines!(ax1, [Jtube]; color = ink2, linewidth = 1, linestyle = :dot)
text!(ax1, 0.013, Jtube * 1.12; text = "J_tube (Bredt closed cell) = $(round(Jtube, digits = 3)) in⁴", fontsize = 12, color = ink2, align = (:left, :bottom))
hlines!(ax1, [2JC]; color = ink2, linewidth = 1, linestyle = :dot)
text!(ax1, 0.013, 2JC * 1.15; text = "2 J_C (no welds)", fontsize = 12, color = ink2, align = (:left, :bottom))
for (name, col, mk, lab) in (("A", blue, :circle, "Ferrite, 3 in welds"), ("B", orange, :rect, "Ferrite, 6 in welds"), ("C", aqua, :utriangle, "Ferrite, 1.5 in welds"))
    m = ser .== name
    idx = sortperm(ratio[m])
    scatterlines!(ax1, ratio[m][idx], Jf[m][idx]; color = col, marker = mk, markersize = 11, linewidth = 1.5, label = lab)
end
mD = ser .== "D"
scatter!(ax1, ratio[mD], Jf[mD]; color = violet, marker = :diamond, markersize = 12, label = "Ferrite, 3 in welds at 18 in, L = 57 and 219 in")
for i in findall(mD)
    text!(ax1, ratio[i] * 1.06, Jf[i]; text = "L = $(Int(L[i]))", fontsize = 11, color = violet, align = (:left, :center))
end
axislegend(ax1; position = :rb, framevisible = false, labelsize = 12)

ax2 = Axis(fig[1, 2]; xscale = log10, xlabel = "weld spacing w_spacing (in)", ylabel = "J_eff (Ferrite) ÷ J_combined (equation)",
           title = "Ratio of the shell-model J to the equation", titlealign = :left, titlesize = 15,
           xgridcolor = grid, ygridcolor = grid, xticks = ([3, 6, 9, 12, 18, 27, 36, 54, 108], string.([3, 6, 9, 12, 18, 27, 36, 54, 108])))
hlines!(ax2, [1.0]; color = ink, linewidth = 1)
for (name, col, mk, lab) in (("A", blue, :circle, "3 in welds"), ("B", orange, :rect, "6 in welds"), ("C", aqua, :utriangle, "1.5 in welds"))
    m = ser .== name
    idx = sortperm(ws[m])
    scatterlines!(ax2, ws[m][idx], (Jf ./ Jeq)[m][idx]; color = col, marker = mk, markersize = 11, linewidth = 1.5, label = lab)
end
scatter!(ax2, ws[mD], (Jf ./ Jeq)[mD]; color = violet, marker = :diamond, markersize = 12, label = "L = 57, 219 in")
axislegend(ax2; position = :rt, framevisible = false, labelsize = 12)
Label(fig[0, :], "Two-C welded OneRack column: Ferrite.jl J_eff = T L / (G βo) vs. proposed J_combined equation (J_tube = Bredt, J_front = J_back = single-C J)",
      fontsize = 17, font = :bold, halign = :left, color = ink)
save(joinpath(@__DIR__, "weld_spacing_summary.png"), fig; px_per_unit = 2)
println("wrote weld_spacing_summary.png")
