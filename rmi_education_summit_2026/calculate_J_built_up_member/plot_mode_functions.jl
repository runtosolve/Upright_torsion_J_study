# plot_mode_functions.jl — shared mode-shape figure code (included by plot_mode_wglmakie.jl and show_mode_wglmakie.jl)
using DelimitedFiles

function load_mode(path)
    d = readdlm(path, ','; skipstart = 1)
    s = Int.(d[:, 1]); i = Int.(d[:, 2]); j = Int.(d[:, 3])
    xyz = Float64.(d[:, 4:6]); u = Float64.(d[:, 7:9])
    nn = maximum(i); nz = maximum(j)
    idx = Dict((s[k], i[k], j[k]) => k for k in eachindex(s))
    faces = Vector{Int}[]
    for jj in 1:nz-1, ss in 1:2, ii in 1:nn-1
        a = idx[(ss, ii, jj)]; b = idx[(ss, ii, jj + 1)]; c = idx[(ss, ii + 1, jj + 1)]; e = idx[(ss, ii + 1, jj)]
        push!(faces, [a, b, c]); push!(faces, [a, c, e])
    end
    return xyz, u, reduce(vcat, permutedims.(faces)), nz
end

function mode_figure(path, title; scale_frac = 0.035, weld_length = 3.0, weld_spacing = 18.0, B = 3.0, D = 3.0 - 0.074, R = 0.199)
    xyz, u, F, nz = load_mode(path)
    L = maximum(xyz[:, 3])
    umag = vec(sqrt.(sum(u .^ 2; dims = 2)))
    sc = scale_frac * L / maximum(umag)
    def = xyz .+ sc .* u
    tri = [GLTriangleFace(F[k, 1], F[k, 2], F[k, 3]) for k in 1:size(F, 1)]
    m0 = Mesh([Point3f(xyz[k, 1], xyz[k, 2], xyz[k, 3]) for k in 1:size(xyz, 1)], tri)
    m1 = Mesh([Point3f(def[k, 1], def[k, 2], def[k, 3]) for k in 1:size(def, 1)], tri)

    fig = Figure(size = (1500, 850), backgroundcolor = :white, fontsize = 14)
    Label(fig[0, 1:3], title; fontsize = 16, font = :bold, halign = :left)
    ax = Axis3(fig[1, 1]; aspect = :data, azimuth = 1.3π, elevation = 0.18, perspectiveness = 0.3,
               xlabel = "X (in)", ylabel = "Y (in)", zlabel = "z (in)")
    mesh!(ax, m1; color = umag ./ maximum(umag), colormap = :viridis, shading = NoShading)
    wireframe!(ax, m1; color = (:black, 0.35), linewidth = 0.5)          # element mesh of the deformed shell
    hidedecorations!(ax); hidespines!(ax)

    # welds: rigid-tie nodes (C1 lip/flange corners and C2 web corners) within ± weld_length/2 of each weld center
    d = readdlm(path, ','; skipstart = 1)
    sid = Int.(d[:, 1]); xa = Float64.(d[:, 4]); ya = Float64.(d[:, 5]); za = Float64.(d[:, 6])
    xp = xa .- (sid .- 1) .* B                                   # part coordinates (C2 shifted back by B)
    nw = floor(Int, (L - weld_length) / weld_spacing) + 1
    centers = range(L / 2 - (nw - 1) * weld_spacing / 2, L / 2 + (nw - 1) * weld_spacing / 2, nw)
    near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
    corner = ((sid .== 1) .& near.(xp, B, R) .& (near.(ya, D, R) .| near.(ya, 0.0, R))) .|
             ((sid .== 2) .& near.(xp, 0.0, R) .& (near.(ya, D, R) .| near.(ya, 0.0, R)))
    inweld = [any(near(z, c, weld_length / 2) for c in centers) for z in za]
    w = findall(corner .& inweld .& (za .> 0) .& (za .< L))
    scatter!(ax, [Point3f(def[k, 1], def[k, 2], def[k, 3]) for k in w]; color = "#e34948", markersize = 7, label = "welds")

    # applied load: uniform axial compression at both ends (arrows on every third end node)
    ends = findall((za .≈ 0.0) .| (za .≈ L))
    ends = ends[1:4:end]
    pts = [Point3f(xyz[k, 1], xyz[k, 2], xyz[k, 3] + (xyz[k, 3] > L / 2 ? 2.5 : -2.5)) for k in ends]
    dirs = [Point3f(0, 0, xyz[k, 3] > L / 2 ? -2.2 : 2.2) for k in ends]
    if isdefined(WGLMakie.Makie, :arrows3d!)
        WGLMakie.Makie.arrows3d!(ax, pts, dirs; color = "#eb6834", shaftradius = 0.03, tipradius = 0.09, tiplength = 0.35)
    else
        arrows!(ax, pts, dirs; color = "#eb6834", linewidth = 0.06, tipradius = 0.09, tiplength = 0.35)
    end
    Colorbar(fig[1, 2]; colormap = :viridis, limits = (0, 1), label = "normalized displacement magnitude", height = Relative(0.5))

    # mid-length cross section: undeformed and deformed (same amplification)
    jm = (nz + 1) ÷ 2
    d = readdlm(path, ','; skipstart = 1)
    sel = Int.(d[:, 3]) .== jm
    ax2 = Axis(fig[1, 3]; aspect = DataAspect(), xlabel = "X (in)", ylabel = "Y (in)")
    for ss in 1:2
        rows = findall(sel .& (Int.(d[:, 1]) .== ss))
        rows = rows[sortperm(Int.(d[rows, 2]))]
        x0 = Float64.(d[rows, 4]); y0 = Float64.(d[rows, 5]); ux = Float64.(d[rows, 7]); uy = Float64.(d[rows, 8])
        lines!(ax2, x0, y0; color = :gray55, linewidth = 3)
        lines!(ax2, x0 .+ sc .* ux, y0 .+ sc .* uy; color = "#2a78d6", linewidth = 3)
    end
    colsize!(fig.layout, 1, Relative(0.45)); colsize!(fig.layout, 3, Relative(0.45))
    return fig
end

