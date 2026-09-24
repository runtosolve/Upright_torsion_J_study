# perforated_J_study.jl
#
# Static twist J study (Moen 2008 §4.2.7.3.2.3) with the upright's perforation pattern, for one C and for the r5
# welded pair (7 × 3 in welds at 18 in, L = 111 in). Holes are modelled by removing shell elements whose centre
# lies inside the hole outline, on a refined extruded mesh.
#
# Perforation pattern (strip layout drawing, 14 ga):
#   web     teardrop holes, two rows 0.797 in from each outside face of the web (1.406 in apart), every 2.0 in
#           along the member; teardrop = Ø0.719 circle + Ø0.375 lobe pointing along the member (≈ 0.90 in long)
#   flanges square holes 0.562 × 0.562, centred 0.856 in from the web outside face, every 2.0 in, on both flanges,
#           placed 1.0 in (half a pitch) along the member from the teardrops (assumed stagger)
#
# Run:  julia --project=. perforated_J_study.jl

include(joinpath(@__DIR__, "calculate_J_built_up_member.jl"))

const hole = (pitch = 2.0, z0 = 1.0,                                   # teardrop centres at z0 + k·pitch
              td_y_from_face = 0.797, td_r_big = 0.719 / 2, td_r_small = 0.375 / 2, td_lobe = 0.546 - 0.375 / 2,
              sq = 0.562, sq_x_from_face = 0.856, sq_z_offset = 1.0)

"""Is the point (Xp, Yp) of the part centerline (min-shifted coords: web at X = 0, flanges at Y = 0 and Y = Ymax),
at member coordinate z, inside a hole?"""
function in_hole(Xp, Yp, z, Ymax, L)
    t = shape.t
    z < 0.6 || z > L - 0.6 && return false                                # keep the end sections intact
    zt = hole.z0 + hole.pitch * round((z - hole.z0) / hole.pitch)        # nearest teardrop centre
    if abs(Xp) < 2e-3                                                     # web flat (generator noise ~2e-4)
        for yc in (hole.td_y_from_face - t / 2, Ymax - (hole.td_y_from_face - t / 2))
            dy = Yp - yc
            hypot(dy, z - zt) <= hole.td_r_big && return true
            hypot(dy, z - (zt + hole.td_lobe)) <= hole.td_r_small && return true
            if zt <= z <= zt + hole.td_lobe                               # tangent hull between the two circles
                r = hole.td_r_big + (hole.td_r_small - hole.td_r_big) * (z - zt) / hole.td_lobe
                abs(dy) <= r && return true
            end
        end
    elseif abs(Yp) < 2e-3 || abs(Yp - Ymax) < 2e-3                        # flange flats
        xc = hole.sq_x_from_face - t / 2
        zs = hole.z0 + hole.sq_z_offset + hole.pitch * round((z - hole.z0 - hole.sq_z_offset) / hole.pitch)
        abs(Xp - xc) <= hole.sq / 2 && abs(z - zs) <= hole.sq / 2 && return true
    end
    return false
end

"""Extruded mesh of `nshapes` C's (shape 2 offset by B) with holes removed; returns grid, id, node presence."""
function perforated_grid(X, Y, Z, nshapes; perforated = true)
    nn = length(X); nz = length(Z); Ymax = maximum(Y)
    nodes = [Node(Vec((X[i] + (s - 1) * shape.B, Y[i], Z[j]))) for j in 1:nz for s in 1:nshapes for i in 1:nn]
    id(s, i, j) = (j - 1) * nshapes * nn + (s - 1) * nn + i
    cells = Quadrilateral[]
    nremoved = 0
    for j in 1:nz-1, s in 1:nshapes, i in 1:nn-1
        xc = (X[i] + X[i+1]) / 2; yc = (Y[i] + Y[i+1]) / 2; zc = (Z[j] + Z[j+1]) / 2
        if perforated && in_hole(xc, yc, zc, Ymax, Z[end])
            nremoved += 1; continue
        end
        push!(cells, Quadrilateral((id(s, i, j), id(s, i, j + 1), id(s, i + 1, j + 1), id(s, i + 1, j))))
    end
    return Grid(cells, nodes), id, nremoved
end

"""Static twist (Fig. 4.41) for 1 or 2 C's, welds as rigid ties (2 C's), J = T L / (G βo)."""
function twist_J(; nshapes = 1, perforated = true, L = shape.length, n_flat = 12, n_corner = 4, dz = 0.2,
                   weld_length = 3.0, weld_locations = collect(range(1.5, L - 1.5, 7)), center = :shear_center)
    X, Y = centerline(shape; n_flat, n_corner); nn = length(X)
    Z = collect(range(0.0, L, Int(round(L / dz)) + 1)); nz = length(Z)
    grid, id, nremoved = perforated_grid(X, Y, Z, nshapes; perforated)
    dh, nd = shell_dofhandler(grid)
    present(n) = haskey(nd, n)
    sp = section_properties(X, Y, shape.t)
    Xc, Yc = nshapes == 1 ? (center == :shear_center ? (sp.xs, sp.ys) : (sp.xc, sp.yc)) : (shape.B / 2, shape.D / 2)
    R = shape.R; B = shape.B; D = shape.D
    near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
    i_webmid = argmin(X .^ 2 .+ (Y .- D / 2) .^ 2)

    ch = ConstraintHandler(dh)
    end0 = Set(id(s, i, 1) for s in 1:nshapes for i in 1:nn); endL = Set(id(s, i, nz) for s in 1:nshapes for i in 1:nn)
    add!(ch, Dirichlet(:u, end0, (x, t) -> [0.0, 0.0], [1, 2]))
    add!(ch, Dirichlet(:u, Set([id(nshapes, i_webmid, 1)]), (x, t) -> [0.0], [3]))
    add!(ch, Dirichlet(:u, endL, (x, t) -> [-βo * (x[2] - Yc), βo * (x[1] - Xc)], [1, 2]))
    n_ties = 0
    if nshapes == 2
        c1_top = findall(i -> near(X[i], B, R) && near(Y[i], D, R), 1:nn); c1_bot = findall(i -> near(X[i], B, R) && near(Y[i], 0.0, R), 1:nn)
        c2_top = findall(i -> near(X[i], 0.0, R) && near(Y[i], D, R), 1:nn); c2_bot = findall(i -> near(X[i], 0.0, R) && near(Y[i], 0.0, R), 1:nn)
        claimed = falses(nz)
        for zw in weld_locations
            st = [j for j in 2:nz-1 if abs(Z[j] - zw) <= weld_length / 2 + 1e-9 && !claimed[j]]; claimed[st] .= true
            for (c1, c2) in ((c1_top, c2_top), (c1_bot, c2_bot))
                slaves = filter(present, vcat([id(1, i, j) for j in st for i in c1], [id(2, i, j) for j in st for i in c2]))
                add_rigid_tie!(ch, nd, grid, slaves[1], slaves); n_ties += length(slaves) - 1
            end
        end
    end
    close!(ch); update!(ch, 0.0)
    K = allocate_matrix(dh, ch)
    K = Q.assemble_global_Ke!(K, dh, QuadratureRule{RefQuadrilateral}(2), QuadratureRule{RefQuadrilateral}(3), IP4(), IP6(), E, ν, shape.t)
    K0 = copy(K); f = zeros(ndofs(dh)); apply!(K, f, ch); u = K \ f; apply!(u, ch)
    Rv = K0 * u
    T = 0.0; Fx = 0.0; Fy = 0.0
    for n in endL
        x, y = grid.nodes[n].x[1], grid.nodes[n].x[2]
        T += (x - Xc) * Rv[nd[n][2]] - (y - Yc) * Rv[nd[n][1]]; Fx += Rv[nd[n][1]]; Fy += Rv[nd[n][2]]
    end
    J = T * L / (G * βo)
    # twist profile (least-squares rigid rotation of the present nodes, shape 1)
    β = zeros(nz)
    for j in 1:nz
        ns = filter(present, [id(1, i, j) for i in 1:nn]); xs = [grid.nodes[n].x[1] for n in ns]; ys = [grid.nodes[n].x[2] for n in ns]
        β[j] = ls_rotation(xs, ys, [u[nd[n][1]] for n in ns], [u[nd[n][2]] for n in ns])
    end
    _, slope, R2, dev = linear_fit(Z, β)
    Aremoved = nremoved / (nshapes * (nn - 1) * (nz - 1))
    @printf("  %s, %s: %d cells (%d removed = %.1f%% of shell area), %d dofs, %d weld ties; T = %.4e lbf-in, |ΣF|/T = %.1e; J = %.4e in⁴; β(z) fit R² = %.6f\n",
            nshapes == 1 ? "single C" : "welded pair", perforated ? "perforated" : "no holes", getncells(grid), nremoved, 100Aremoved, ndofs(dh), n_ties, T, hypot(Fx, Fy) / abs(T), J, R2)
    return (; J, T, Z, β, nremoved, ncells = getncells(grid), ndofs = ndofs(dh))
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Perforation pattern: teardrops (Ø$(2hole.td_r_big) + Ø$(2hole.td_r_small) lobe) in the web at $(hole.td_y_from_face) in from each face every $(hole.pitch) in; $(hole.sq) in squares in both flanges at $(hole.sq_x_from_face) in from the web face every $(hole.pitch) in, offset $(hole.sq_z_offset) in\n")
    results = Dict{String,Any}()
    for (label, kw) in (("single C, no holes",        (; nshapes = 1, perforated = false)),
                        ("single C, perforated",      (; nshapes = 1, perforated = true)),
                        ("welded pair r5, no holes",  (; nshapes = 2, perforated = false)),
                        ("welded pair r5, perforated",(; nshapes = 2, perforated = true)))
        println(label); t_run = @elapsed r = twist_J(; kw...); @printf("  (%.0f s)\n", t_run); results[label] = r
    end
    J1 = results["single C, no holes"].J; J1h = results["single C, perforated"].J
    J2 = results["welded pair r5, no holes"].J; J2h = results["welded pair r5, perforated"].J
    @printf("\nSingle C:     J = %.4e (no holes)  →  %.4e (perforated), ratio %.3f   [exact Saint-Venant gross J = 1.3349e-3]\n", J1, J1h, J1h / J1)
    @printf("Welded pair:  J_eff = %.4f (no holes)  →  %.4f (perforated), ratio %.3f\n", J2, J2h, J2h / J2)
    open(joinpath(@__DIR__, "perforated_J_results.csv"), "w") do io
        println(io, "case,J_in4,T_lbf_in,cells,cells_removed,ndofs")
        for (k, r) in results; println(io, "\"$k\",$(r.J),$(r.T),$(r.ncells),$(r.nremoved),$(r.ndofs)"); end
    end
    r = results["welded pair r5, perforated"]
    writedlm(joinpath(@__DIR__, "twist_profile_built_up_perforated.csv"), ["z_in" "beta_C1_over_betao"; r.Z r.β ./ βo], ',')
    println("Wrote perforated_J_results.csv, twist_profile_built_up_perforated.csv")
end
