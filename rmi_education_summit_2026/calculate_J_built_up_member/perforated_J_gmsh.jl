# perforated_J_gmsh.jl
#
# Static twist J study (Moen 2008 §4.2.7.3.2.3) of the perforated upright with a Gmsh mesh that follows the true
# hole outlines. The developed (unfolded) strip, width S = centerline length, is drawn in Gmsh (OCC) with the teardrop
# and square holes cut out, fragmented by the bend lines at every centerline vertex so the mesh conforms to the
# folds, meshed with recombination (quads where possible, triangles at the hole boundaries), then folded onto the C:
# (s, z) → (X(s), Y(s), z). Quads are assembled with QuadShellFiniteElement.jl, triangles with
# TriShellFiniteElement.jl (both with the Hughes–Brezzi drilling term) through Ferrite SubDofHandlers.
#
# Perforations (14 ga strip layout): web teardrops Ø0.719 + Ø0.375 lobe, rows 0.797 in from each web face, every
# 2.0 in; flange squares 0.562 × 0.562 (R0.02), 0.856 in from the web face on both flanges, every 2.0 in, offset
# 1.0 in from the teardrops (assumed stagger).
#
# Run:  julia --project=. perforated_J_gmsh.jl        (writes perforated_J_gmsh_results.csv)

include(joinpath(@__DIR__, "calculate_J_built_up_member.jl"))
using FerriteGmsh; import FerriteGmsh: gmsh
using TriShellFiniteElement
const TS = TriShellFiniteElement

const hole = (pitch = 2.0, z0 = 1.0, td_y_from_face = 0.797, td_r_big = 0.719 / 2, td_r_small = 0.375 / 2,
              td_lobe = 0.546 - 0.375 / 2, sq = 0.562, sq_r = 0.02, sq_x_from_face = 0.856, sq_z_offset = 1.0)

# ---- centerline arc-length parametrization (min-shifted part coordinates: web at X = 0)
function centerline_param(; n_flat = 4, n_corner = 4)
    X, Y = centerline(shape; n_flat, n_corner)
    s = cumsum([0.0; hypot.(diff(X), diff(Y))])
    fold(sv) = begin                                            # (X, Y) at arc length sv, piecewise linear
        k = clamp(searchsortedlast(s, sv), 1, length(s) - 1)
        w = (sv - s[k]) / (s[k+1] - s[k])
        (X[k] + w * (X[k+1] - X[k]), Y[k] + w * (Y[k+1] - Y[k]))
    end
    # arc length of the nearest centerline point to (xp, yp)
    s_of(xp, yp) = begin
        best = (Inf, 0.0)
        for k in 1:length(s)-1
            dx, dy = X[k+1] - X[k], Y[k+1] - Y[k]; l2 = dx^2 + dy^2
            w = clamp(((xp - X[k]) * dx + (yp - Y[k]) * dy) / l2, 0.0, 1.0)
            d = hypot(xp - (X[k] + w * dx), yp - (Y[k] + w * dy))
            d < best[1] && (best = (d, s[k] + w * sqrt(l2)))
        end
        best[2]
    end
    return X, Y, s, fold, s_of
end

# ---- Gmsh mesh of the developed strip [0, S] × [0, L] with holes and fold lines; returns a 2D Ferrite grid (s, z)
function strip_mesh(L; perforated = true, mesh_size = 0.2, n_flat = 4, n_corner = 4)
    X, Y, s, fold, s_of = centerline_param(; n_flat, n_corner)
    S = s[end]; t = shape.t; Ymax = maximum(Y)
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("strip"); occ = gmsh.model.occ
    rect = occ.addRectangle(0.0, 0.0, 0.0, S, L)
    holes = Int[]
    if perforated
        for zt in (hole.z0):(hole.pitch):(L - 0.6), sy in (s_of(0.0, hole.td_y_from_face - t / 2), s_of(0.0, Ymax - (hole.td_y_from_face - t / 2)))
            big = occ.addDisk(sy, zt, 0.0, hole.td_r_big, hole.td_r_big)
            small = occ.addDisk(sy, zt + hole.td_lobe, 0.0, hole.td_r_small, hole.td_r_small)
            p = [occ.addPoint(sy - hole.td_r_big, zt, 0.0), occ.addPoint(sy + hole.td_r_big, zt, 0.0),
                 occ.addPoint(sy + hole.td_r_small, zt + hole.td_lobe, 0.0), occ.addPoint(sy - hole.td_r_small, zt + hole.td_lobe, 0.0)]
            l = [occ.addLine(p[1], p[2]), occ.addLine(p[2], p[3]), occ.addLine(p[3], p[4]), occ.addLine(p[4], p[1])]
            hull = occ.addPlaneSurface([occ.addCurveLoop(l)])
            fused, _ = occ.fuse([(2, big)], [(2, small), (2, hull)])
            append!(holes, [tag for (d, tag) in fused])
        end
        xc = hole.sq_x_from_face - t / 2
        for zs in (hole.z0 + hole.sq_z_offset):(hole.pitch):(L - 0.6), sx in (s_of(xc, 0.0), s_of(xc, Ymax))
            push!(holes, occ.addRectangle(sx - hole.sq / 2, zs - hole.sq / 2, 0.0, hole.sq, hole.sq, -1, hole.sq_r))
        end
    end
    surf = [(2, rect)]
    if !isempty(holes)
        surf, _ = occ.cut([(2, rect)], [(2, h) for h in holes])
    end
    # fold lines at every interior centerline vertex
    lines = [(1, occ.addLine(occ.addPoint(s[k], 0.0, 0.0), occ.addPoint(s[k], L, 0.0))) for k in 2:length(s)-1]
    frag, _ = occ.fragment(surf, lines)
    occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size); gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size / 3)
    gmsh.option.setNumber("Mesh.Algorithm", 6)                 # Frontal-Delaunay
    gmsh.option.setNumber("Mesh.RecombineAll", 1); gmsh.option.setNumber("Mesh.RecombinationAlgorithm", 1)
    gmsh.model.mesh.generate(2)
    grid2d = togrid()
    gmsh.finalize()
    return grid2d, fold, S
end

# ---- fold into 3D and duplicate for the second C; mixed Triangle/Quadrilateral grid with cell sets
function folded_grid(grid2d, fold, nshapes)
    n2 = getnnodes(grid2d)
    nodes = Node{3,Float64}[]
    for s in 1:nshapes, nd in grid2d.nodes
        x, y = fold(nd.x[1]); push!(nodes, Node(Vec((x + (s - 1) * shape.B, y, nd.x[2]))))
    end
    cells = Union{Triangle,Quadrilateral}[]
    for s in 1:nshapes, c in grid2d.cells
        off = (s - 1) * n2
        push!(cells, c isa Triangle ? Triangle(map(i -> i + off, c.nodes)) : Quadrilateral(map(i -> i + off, c.nodes)))
    end
    grid = Grid(cells, nodes)
    addcellset!(grid, "tri", Set(i for (i, c) in enumerate(cells) if c isa Triangle))
    addcellset!(grid, "quad", Set(i for (i, c) in enumerate(cells) if c isa Quadrilateral))
    return grid
end

function mixed_dofhandler(grid)
    dh = DofHandler(grid)
    tri = getcellset(grid, "tri"); quad = getcellset(grid, "quad")
    if !isempty(quad)
        sdh = SubDofHandler(dh, quad); add!(sdh, :u, Lagrange{RefQuadrilateral,1}()^3); add!(sdh, :θ, Lagrange{RefQuadrilateral,1}()^3)
    end
    if !isempty(tri)
        sdh = SubDofHandler(dh, tri); add!(sdh, :u, Lagrange{RefTriangle,1}()^3); add!(sdh, :θ, Lagrange{RefTriangle,1}()^3)
    end
    close!(dh)
    nd = Dict{Int,Vector{Int}}()
    for sdh in dh.subdofhandlers, cell in CellIterator(sdh)
        cd = celldofs(cell); k = length(cell.nodes)
        for (i, n) in enumerate(cell.nodes)
            nd[n] = [cd[3(i-1)+1], cd[3(i-1)+2], cd[3(i-1)+3], cd[3k+3(i-1)+1], cd[3k+3(i-1)+2], cd[3k+3(i-1)+3]]
        end
    end
    return dh, nd
end

const IND18 = [1, 2, 3, 7, 8, 9, 13, 14, 15, 4, 5, 6, 10, 11, 12, 16, 17, 18]
function assemble_mixed_Ke!(K, dh, t)
    qr_m = QuadratureRule{RefQuadrilateral}(2); qr_b = QuadratureRule{RefQuadrilateral}(3)
    qr1 = QuadratureRule{RefTriangle}(1); qr3 = QuadratureRule{RefTriangle}(2)
    assembler = start_assemble(K)
    for sdh in dh.subdofhandlers, cell in CellIterator(sdh)
        x = getcoordinates(cell)
        if length(x) == 4
            T = Q.calculation_rotation_matrix(x); xl = Q.global_nodal_coords_to_planar_coords(x, T)
            ke = Q.local_elastic_stiffness_matrix!(qr_m, qr_b, IP4(), IP6(), E, ν, t, xl)
            Te = Q.rotation_matrix_for_element_stiffness_drilling(T)
            keg = (Te * ke * Te')[Q.FIELD_ORDER_24, Q.FIELD_ORDER_24]
        else
            T = TS.calculation_rotation_matrix(x); xl = TS.global_nodal_coords_to_planar_coords(x, T)
            ke = TS.local_elastic_stiffness_matrix!(qr1, qr3, TS.IP3(), TS.IP6(), E, ν, t, xl)
            Te = TS.rotation_matrix_for_element_stiffness_drilling(T)
            keg = (Te * ke * Te')[IND18, IND18]
        end
        assemble!(assembler, celldofs(cell), keg)
    end
    return K
end

"""Static twist (Fig. 4.41) on the Gmsh mesh for 1 or 2 C's; J = T L / (G βo)."""
function twist_J_gmsh(; nshapes = 1, perforated = true, L = shape.length, mesh_size = 0.2, n_flat = 4, n_corner = 4,
                        weld_length = 3.0, weld_locations = collect(range(1.5, L - 1.5, 7)))
    grid2d, fold, S = strip_mesh(L; perforated, mesh_size, n_flat, n_corner)
    grid = folded_grid(grid2d, fold, nshapes)
    dh, nd = mixed_dofhandler(grid)
    X, Y = centerline(shape; n_flat, n_corner); sp = section_properties(X, Y, shape.t)
    Xc, Yc = nshapes == 1 ? (sp.xs, sp.ys) : (shape.B / 2, shape.D / 2)
    R = shape.R; B = shape.B; D = shape.D
    near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
    P = [n.x for n in grid.nodes]
    end0 = Set(i for i in eachindex(P) if near(P[i][3], 0.0, 1e-6)); endL = Set(i for i in eachindex(P) if near(P[i][3], L, 1e-6))
    part(i) = (P[i][1] - (P[i][1] > B + 0.5 ? B : 0.0), P[i][2])            # part coordinates (shape 2 shifted back)
    shapeof(i) = P[i][1] > B + 0.5 ? 2 : 1
    datum = argmin([near(P[i][3], 0.0, 1e-6) && shapeof(i) == nshapes ? part(i)[1]^2 + (part(i)[2] - D / 2)^2 : Inf for i in eachindex(P)])
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, end0, (x, t) -> [0.0, 0.0], [1, 2]))
    add!(ch, Dirichlet(:u, Set([datum]), (x, t) -> [0.0], [3]))
    add!(ch, Dirichlet(:u, endL, (x, t) -> [-βo * (x[2] - Yc), βo * (x[1] - Xc)], [1, 2]))
    n_ties = 0
    if nshapes == 2
        claimed = falses(length(P))                                              # abutting welds: a node belongs to the first weld
        for zw in weld_locations, (yw) in (D, 0.0)
            slaves = [i for i in eachindex(P) if !claimed[i] && abs(P[i][3] - zw) <= weld_length / 2 + 1e-9 && P[i][3] > 1e-6 && P[i][3] < L - 1e-6 &&
                      near(part(i)[2], yw, R) && ((shapeof(i) == 1 && near(part(i)[1], B, R)) || (shapeof(i) == 2 && near(part(i)[1], 0.0, R)))]
            isempty(slaves) && continue
            claimed[slaves] .= true
            add_rigid_tie!(ch, nd, grid, slaves[1], slaves); n_ties += length(slaves) - 1
        end
    end
    close!(ch); update!(ch, 0.0)
    K = allocate_matrix(dh, ch); K = assemble_mixed_Ke!(K, dh, shape.t)
    K0 = copy(K); f = zeros(ndofs(dh)); apply!(K, f, ch); u = K \ f; apply!(u, ch)
    Rv = K0 * u
    T = 0.0; Fx = 0.0; Fy = 0.0
    for n in endL
        T += (P[n][1] - Xc) * Rv[nd[n][2]] - (P[n][2] - Yc) * Rv[nd[n][1]]; Fx += Rv[nd[n][1]]; Fy += Rv[nd[n][2]]
    end
    J = T * L / (G * βo)
    ntri = length(getcellset(grid, "tri")); nquad = length(getcellset(grid, "quad"))
    @printf("  %s, %s (Gmsh h = %.2f): %d quads + %d tris, %d dofs, %d weld ties; T = %.4e lbf-in, |ΣF|/T = %.1e; J = %.4e in⁴\n",
            nshapes == 1 ? "single C" : "welded pair", perforated ? "perforated" : "no holes", mesh_size, nquad, ntri, ndofs(dh), n_ties, T, hypot(Fx, Fy) / abs(T), J)
    return (; J, T, nquad, ntri, ndofs = ndofs(dh), grid, u, nd)
end

if abspath(PROGRAM_FILE) == @__FILE__
    results = Dict{String,Any}()
    for (label, kw) in (("single C, no holes",         (; nshapes = 1, perforated = false)),
                        ("single C, perforated",       (; nshapes = 1, perforated = true)),
                        ("welded pair r5, no holes",   (; nshapes = 2, perforated = false)),
                        ("welded pair r5, perforated", (; nshapes = 2, perforated = true)))
        println(label); t_run = @elapsed r = twist_J_gmsh(; kw...); @printf("  (%.0f s)\n", t_run); results[label] = r
    end
    J1 = results["single C, no holes"].J; J1h = results["single C, perforated"].J
    J2 = results["welded pair r5, no holes"].J; J2h = results["welded pair r5, perforated"].J
    @printf("\nSingle C:     J = %.4e (no holes) → %.4e (perforated), ratio %.3f   [structured quad mesh: 1.3271e-3 → 1.2146e-3, 0.915; exact gross 1.3349e-3]\n", J1, J1h, J1h / J1)
    @printf("Welded pair:  J_eff = %.4f (no holes) → %.4f (perforated), ratio %.3f   [structured quad mesh: 0.7531 → 0.5915, 0.785]\n", J2, J2h, J2h / J2)
    open(joinpath(@__DIR__, "perforated_J_gmsh_results.csv"), "w") do io
        println(io, "case,J_in4,T_lbf_in,quads,tris,ndofs")
        for (k, r) in results; println(io, "\"$k\",$(r.J),$(r.T),$(r.nquad),$(r.ntri),$(r.ndofs)"); end
    end
end
