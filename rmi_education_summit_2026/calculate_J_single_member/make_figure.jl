# make_figure.jl — twist profile (Moen 2008 Fig. 4.42 analog) and J summary from the CSVs written by
# calculate_J_single_member.jl.   Run:  julia --project=. make_figure.jl
using CairoMakie, DelimitedFiles
blue = "#2a78d6"; orange = "#eb6834"; ink = "#0b0b0b"; ink2 = "#52514e"; grid = "#e6e5e1"

prof = readdlm(joinpath(@__DIR__, "twist_profile_abaqus_mesh.csv"), ','; skipstart = 1)
z = Float64.(prof[:, 1]); βls = Float64.(prof[:, 2]); βweb = Float64.(prof[:, 3])
res = readdlm(joinpath(@__DIR__, "J_results.csv"), ','; skipstart = 1)
drilling = String.(strip.(res[:, 2])); ratio = Float64.(res[:, 16]); Jthin_ratio = Float64.(res[:, 15]) ./ Float64.(res[:, 14])
labels = ["Abaqus mesh 4/flat 5/corner 222 along, Hughes–Brezzi", "same, Cs = 0", "same, twist about centroid", "same, L = 55.5 in",
          "coarse 2/flat 2/corner 111 along, Hughes–Brezzi", "fine 8/flat 10/corner 444 along, Hughes–Brezzi",
          "Abaqus mesh, original drilling penalty 1/100", "fine mesh, original drilling penalty 1/100"]
@assert length(labels) == size(res, 1)

fig = Figure(size = (1500, 560), fontsize = 15, backgroundcolor = :white)
ax1 = Axis(fig[1, 1]; xlabel = "distance along upright z (in)", ylabel = "β / βo, normalized angle of twist",
           title = "Twist is linear along the warping-free upright",
           titlealign = :left, titlesize = 15, xgridcolor = grid, ygridcolor = grid, xticks = 0:20:120)
lines!(ax1, z, βweb; color = blue, linewidth = 2, label = "relative rotation of web/flange corners")
scatter!(ax1, z[1:12:end], βls[1:12:end]; color = :white, strokecolor = blue, strokewidth = 2, markersize = 9, label = "least-squares rigid rotation of section")
lines!(ax1, [0, z[end]], [0, 1]; color = ink2, linewidth = 1, linestyle = :dash, label = "β = βo z / L")
axislegend(ax1; position = :lt, framevisible = false, labelsize = 13)
text!(ax1, 60, 0.25; text = "R² = 1.00000000\nmax |β − fit| / βo = 1.7e-5\nJ = T L / (G βo) = 1.3256e-3 in⁴", color = ink2, fontsize = 13, align = (:left, :top))

ax2 = Axis(fig[1, 2]; xlabel = "J from Ferrite shell model ÷ exact Saint-Venant J (1.3349e-3 in⁴)",
           title = "Hughes–Brezzi drilling (now default): J within 1 % of exact; original penalty: 2 to 3.4 × too stiff", titlealign = :left, titlesize = 15,
           yticks = (1:length(labels), reverse(labels)), xgridcolor = grid, ygridcolor = grid, yticklabelsize = 12,
           xticks = 0:0.5:3.5)
vlines!(ax2, [1.0]; color = ink, linewidth = 1)
vlines!(ax2, [1.0 / Jthin_ratio[1]]; color = ink2, linewidth = 1, linestyle = :dot)
n = length(labels)
for (k, (r, d)) in enumerate(zip(ratio, drilling))
    y = n - k + 1
    c = d == "hughes_brezzi" ? blue : orange
    lines!(ax2, [1.0, r], [y, y]; color = c, linewidth = 2)
    scatter!(ax2, [r], [y]; color = c, markersize = 11)
    text!(ax2, r + (r < 0.95 ? -0.04 : 0.04), y; text = string(round(r, digits = 3)), align = (r < 0.95 ? :right : :left, :center), fontsize = 12, color = ink2)
end
text!(ax2, 1.0 / Jthin_ratio[1] + 0.02, 0.45; text = "thin-walled Σbt³/3", fontsize = 12, color = ink2, align = (:left, :center))
xlims!(ax2, 0.7, 3.75); ylims!(ax2, 0.3, n + 0.7)
Label(fig[0, :], "Saint-Venant J of one 3 × 3 × 0.074 in lipped C upright (L = 111 in) by the Moen (2008) §4.2.7.3.2.3 static twist method, Ferrite.jl + QuadShellFiniteElement.jl",
      fontsize = 17, font = :bold, halign = :left, color = ink)
colsize!(fig.layout, 1, Relative(0.42))
save(joinpath(@__DIR__, "J_single_member_summary.png"), fig; px_per_unit = 2)
println("wrote J_single_member_summary.png")
