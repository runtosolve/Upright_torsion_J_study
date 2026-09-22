# make_figure_built_up.jl — twist profiles of the two C shapes along the welded built-up member and the
# J_eff comparison, from the CSVs written by calculate_J_built_up_member.jl.
using CairoMakie, DelimitedFiles
blue = "#2a78d6"; orange = "#eb6834"; aqua = "#1baf7a"; ink = "#0b0b0b"; ink2 = "#52514e"; grid = "#e6e5e1"; band = "#f0efec"

prof = readdlm(joinpath(@__DIR__, "twist_profile_built_up.csv"), ','; skipstart = 1)
z = Float64.(prof[:, 1]); β1 = Float64.(prof[:, 2]); β2 = Float64.(prof[:, 3]); welded = Bool.(prof[:, 4] .== 1)
res = readdlm(joinpath(@__DIR__, "J_built_up_results.csv"), ','; skipstart = 1)
labels = String.(strip.(res[:, 1])); Jeff = Float64.(res[:, 6]); Jbox = Float64(res[1, 8]); twoJC = Float64(res[1, 9]); Jab = Float64(res[1, 10])

fig = Figure(size = (1500, 600), fontsize = 15, backgroundcolor = :white)
ax1 = Axis(fig[1, 1]; xlabel = "distance along member z (in)", ylabel = "β / βo, normalized angle of twist",
           title = "Twist of each C along the member (welded zones shaded)", titlealign = :left, titlesize = 15,
           xgridcolor = grid, ygridcolor = grid, xticks = 0:20:120)
# weld bands
function weld_bands!(ax, z, welded, color)
    j = 1
    while j <= length(z)
        if welded[j]
            k = j
            while k < length(z) && welded[k+1]; k += 1; end
            vspan!(ax, z[j] - 0.25, z[k] + 0.25; color)
            j = k + 1
        else
            j += 1
        end
    end
end
weld_bands!(ax1, z, welded, band)
lines!(ax1, [0, z[end]], [0, 1]; color = ink2, linewidth = 1, linestyle = :dash, label = "uniform twist βo z / L")
lines!(ax1, z, β1; color = blue, linewidth = 2, label = "C1 (web at X = 0, lips welded to C2 web)")
lines!(ax1, z, β2; color = orange, linewidth = 2, linestyle = :dot, label = "C2 (offset by B = 3 in)")
axislegend(ax1; position = :lt, framevisible = false, labelsize = 13)
text!(ax1, 62, 0.22; text = "J_eff = T L / (G βo) = $(round(Jeff[1], sigdigits = 4)) in⁴\nwelds: 7 × 3 in at 18 in spacing", color = ink2, fontsize = 13, align = (:left, :top))

ax2 = Axis(fig[1, 2]; xlabel = "J_eff (in⁴), log scale", xscale = log10,
           title = "Effective J: Ferrite shell model vs references", titlealign = :left, titlesize = 15,
           yticks = (1:length(labels), reverse(labels)), xgridcolor = grid, ygridcolor = grid, yticklabelsize = 12)
n = length(labels)
for (k, J) in enumerate(Jeff)
    y = n - k + 1
    scatter!(ax2, [J], [y]; color = blue, markersize = 11)
    text!(ax2, J * 1.15, y; text = string(round(J, sigdigits = 4)), align = (:left, :center), fontsize = 12, color = ink2)
end
vlines!(ax2, [twoJC]; color = ink2, linewidth = 1, linestyle = :dot)
vlines!(ax2, [Jbox]; color = ink2, linewidth = 1, linestyle = :dot)
vlines!(ax2, [Jab]; color = orange, linewidth = 1.5, linestyle = :dash)
text!(ax2, twoJC * 1.1, 0.55; text = "2 J_C (no welds)", fontsize = 12, color = ink2, align = (:left, :center))
text!(ax2, Jbox * 1.1, 0.55; text = "Bredt closed cell", fontsize = 12, color = ink2, align = (:left, :center))
text!(ax2, Jab * 1.1, n + 0.45; text = "Abaqus r5 (S4R, contact)", fontsize = 12, color = orange, align = (:left, :center))
xlims!(ax2, 1e-3, 10); ylims!(ax2, 0.2, n + 0.8)
Label(fig[0, :], "Effective torsion constant of the two-C welded OneRack column (L = 111 in) by the Moen (2008) §4.2.7.3.2.3 static twist method, Ferrite.jl",
      fontsize = 17, font = :bold, halign = :left, color = ink)
colsize!(fig.layout, 1, Relative(0.5))
save(joinpath(@__DIR__, "J_built_up_summary.png"), fig; px_per_unit = 2)
println("wrote J_built_up_summary.png")
