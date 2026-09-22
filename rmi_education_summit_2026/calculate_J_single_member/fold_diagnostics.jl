# fold_diagnostics.jl — why the drilling treatment matters for torsion of folded / curved shells.
# Twists simple polyline sections (flat strip, two strips at a fold angle, slit tubes, sharp and rounded
# lipped C) with the three drilling options and reports J relative to the thin-walled Σbt³/3.
# Run:  julia --project=. fold_diagnostics.jl        (writes fold_diagnostics_log.txt)
include(joinpath(@__DIR__, "calculate_J_single_member.jl"))

function torsion_J_XY(X, Y; L = 28.0, n_z = 56, Cs = 0.1, drilling = :hughes_brezzi)
    t = shape.t
    Z = collect(range(0.0, L, n_z + 1)); nn = length(X); nz = length(Z)
    Jt = sum(hypot(X[i+1] - X[i], Y[i+1] - Y[i]) for i in 1:nn-1) * t^3 / 3
    sp = SectionProperties.open_thin_walled(X, Y, fill(t, nn - 1)); Xc, Yc = sp.xs, sp.ys
    grid, id = extruded_grid(X, Y, Z); dh, nd = shell_dofhandler(grid)
    end0 = Set(id(i, 1) for i in 1:nn); endL = Set(id(i, nz) for i in 1:nn)
    K = allocate_matrix(dh)
    K = Q.assemble_global_Ke!(K, dh, QuadratureRule{RefQuadrilateral}(2), QuadratureRule{RefQuadrilateral}(3), IP4(), IP6(), E, ν, t; Cs, drilling)
    K0 = copy(K)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, end0, (x, t_) -> [0.0, 0.0], [1, 2]))
    add!(ch, Dirichlet(:u, Set([id(1, 1)]), (x, t_) -> [0.0], [3]))
    add!(ch, Dirichlet(:u, endL, (x, t_) -> [-βo * (x[2] - Yc), βo * (x[1] - Xc)], [1, 2]))
    close!(ch); update!(ch, 0.0)
    f = zeros(ndofs(dh)); apply!(K, f, ch); u = K \ f; apply!(u, ch)
    R = K0 * u
    T = sum((grid.nodes[n].x[1] - Xc) * R[nd[n][2]] - (grid.nodes[n].x[2] - Yc) * R[nd[n][1]] for n in endL)
    F = hypot(sum(R[nd[n][1]] for n in endL), sum(R[nd[n][2]] for n in endL))
    return T * L / (G * βo) / Jt, F / T
end
seg(x0, y0, x1, y1, n) = ([x0 + (x1 - x0) * k / n for k in 0:n], [y0 + (y1 - y0) * k / n for k in 0:n])
function sections()
    s = Pair{String, Tuple{Vector{Float64}, Vector{Float64}}}[]
    push!(s, "flat 3 in strip, 8 elements across" => seg(0.0, 0.0, 3.0, 0.0, 8))
    for α in (5.0, 30.0, 90.0)
        a = deg2rad(α); X1, Y1 = seg(0.0, 0.0, 1.5, 0.0, 4); X2, Y2 = seg(1.5, 0.0, 1.5 + 1.5cos(a), 1.5sin(a), 4)
        push!(s, "two 1.5 in strips folded $(Int(α)) deg" => ([X1; X2[2:end]], [Y1; Y2[2:end]]))
    end
    X1, Y1 = seg(0.0, 3.0, 0.0, 0.0, 8); X2, Y2 = seg(0.0, 0.0, 3.0, 0.0, 8)
    push!(s, "sharp 3x3 angle, 8 per leg" => ([X1; X2[2:end]], [Y1; Y2[2:end]]))
    Lc = [shape.lip, shape.B, shape.D, shape.B, shape.lip]; θc = [π / 2, π, -π / 2, 0.0, π / 2]
    cs = CrossSectionGeometry.create_thin_walled_cross_section_geometry(Lc, θc, fill(8, 5), shape.t; centerline = "to left", offset = (0.0, 0.0))
    push!(s, "sharp-corner lipped C, 8 per flat" => ([p[1] for p in cs.centerline_node_XY], [p[2] for p in cs.centerline_node_XY]))
    push!(s, "rounded lipped C, 4/flat 5/corner (Abaqus mesh)" => centerline(shape; n_flat = 4, n_corner = 5))
    push!(s, "rounded lipped C, 8/flat 10/corner" => centerline(shape; n_flat = 8, n_corner = 10))
    push!(s, "rounded lipped C, 8/flat 20/corner" => centerline(shape; n_flat = 8, n_corner = 20))
    for (Rr, nseg) in ((shape.R, 30), (1.5, 60))
        φ = range(0, 1.5π, nseg + 1)
        push!(s, "270 deg slit tube R = $(Rr), $(nseg) segments" => (collect(Rr .* cos.(φ)), collect(Rr .* sin.(φ))))
    end
    return s
end
open(joinpath(@__DIR__, "fold_diagnostics_log.txt"), "w") do io
    for out in (stdout, io)
        println(out, "J / (Σbt³/3) for L = 28 in, 56 elements along, Cs = 0.1.  |ΣF|/T is the net in-plane end reaction (0 for pure torsion).")
        @printf(out, "%-50s %22s %22s\n", "section", "absolute penalty 1/100", "Hughes–Brezzi γ = Gt")
    end
    for (name, (X, Y)) in sections()
        vals = String[]
        for d in (:penalty, :hughes_brezzi)
            r, F = torsion_J_XY(X, Y; drilling = d)
            push!(vals, @sprintf("%7.4f  (|ΣF|/T %.0e)", r, F))
        end
        for out in (stdout, io)
            @printf(out, "%-50s %22s %22s\n", name, vals...)
        end
    end
end
