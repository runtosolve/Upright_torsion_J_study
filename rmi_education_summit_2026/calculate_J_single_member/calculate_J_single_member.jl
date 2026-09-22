# calculate_J_single_member.jl
#
# Saint-Venant torsion constant J of a single OneRack upright (the SHAPE_1 lipped C of
# combo_1/J/combo_1_J_r5.jl, i.e. one of the two members in the Abaqus models in combo_1/J/model_runs)
# computed with a Ferrite.jl shell model built from QuadShellFiniteElement.jl, using the static twist
# method of Moen (2008), Section 4.2.7.3.2.3:
#
#     T = G J dβ/dz − E Cw d³β/dz³                                   (Eq. 4.12)
#
# With both ends free to warp, d³β/dz³ = 0, the twist β is linear along the member and
#
#     J = T L / (G βo)          or, equivalently,      J = T / (G dβ/dz)         (Eq. 4.13)
#
# Boundary conditions follow Figure 4.41 of the dissertation:
#   z = 0 : every cross-section node fixed in X and Y (twist and in-plane translation restrained,
#           warping free); the web mid-height node fixed in Z to remove the axial rigid body mode.
#   z = L : cross section rotated rigidly by βo about the shear center (the Abaqus kinematic coupling
#           of the section to the shear-center reference node in dof 1, 2), warping free.
# The torque T is the resultant moment of the reactions at z = L about the shear center.
#
# The twist β(z) is measured two ways at every cross section along the member:
#   (a) as the relative X displacement of the two web/flange corners divided by their distance
#       (the dissertation's "relative rotation of the flange-web corners"), and
#   (b) as the least-squares rigid rotation of all nodes in the cross section.
# J is reported from the end twist (Eq. 4.13) and from the fitted slope dβ/dz, and compared with
#   - the thin-walled formula J = Σ b t³ / 3 (SectionProperties.jl, CUFSM cutwp), and
#   - the exact Saint-Venant J of the real cross section with its rounded corners and thickness,
#     from the Prandtl stress function solved on a 2D mesh of the section (J_saint_venant_2D.jl).
#
# The drilling-dof treatment of the shell element is selectable through the `drilling` keyword of
# QuadShellFiniteElement (added 2026-09-21 as a result of this study): `:hughes_brezzi` (now the
# package default) or `:penalty` (the original MATLAB q42 absolute penalty, which is wrong for torsion
# of folded/curved shells). Both are run here for comparison.
#
# Run from this folder:   julia --project=. calculate_J_single_member.jl

using Ferrite, Tensors, LinearAlgebra, SparseArrays, Statistics, Printf, DelimitedFiles
using QuadShellFiniteElement, CrossSectionGeometry, SectionProperties

include(joinpath(@__DIR__, "J_saint_venant_2D.jl"))     # defines `shape`, strip_grid, saint_venant_J, upright_section

const Q = QuadShellFiniteElement

# ---------------------------------------------------------------------------------------------
# Material, identical to the Abaqus model (steel, E = 29500 ksi, ν = 0.3, G = E/2.6)
# ---------------------------------------------------------------------------------------------
const E = 29_500_000.0            # psi
const ν = 0.30
const G = E / (2 * (1 + ν))
const βo = 0.01                   # applied end twist, rad (linear analysis, magnitude arbitrary)

# ---------------------------------------------------------------------------------------------
# Cross-section centerline (same CrossSectionGeometry call as the Abaqus model generator)
# ---------------------------------------------------------------------------------------------
function centerline(shape; n_flat = 4, n_corner = 5)
    L = [shape.lip, shape.B, shape.D, shape.B, shape.lip]
    θ = [π / 2, π, -π / 2, 0.0, π / 2]
    cs = CrossSectionGeometry.create_thin_walled_cross_section_geometry(L, θ, fill(n_flat, 5), fill(shape.R, 4), fill(n_corner, 4), shape.t;
                                                                         centerline = "to left", offset = (0.0, 0.0))
    X = [p[1] for p in cs.centerline_node_XY]
    Y = [p[2] for p in cs.centerline_node_XY]
    X .-= minimum(X)                      # same origin shift as the Abaqus model generator
    Y .-= minimum(Y)
    return X, Y
end

# thin-walled section properties of the centerline polyline (CUFSM cutwp_prop2)
section_properties(X, Y, t) = SectionProperties.open_thin_walled(X, Y, fill(t, length(X) - 1))

# ---------------------------------------------------------------------------------------------
# Extruded quadrilateral shell mesh with 3D nodes; node and element ordering match
# MeshExtrusions.open_cross_section_with_shell_elements (station-major node numbering)
# ---------------------------------------------------------------------------------------------
function extruded_grid(X, Y, Z)
    nn = length(X)
    nz = length(Z)
    nodes = [Node(Vec((X[i], Y[i], Z[j]))) for j in 1:nz for i in 1:nn]
    id(i, j) = (j - 1) * nn + i
    cells = [Quadrilateral((id(i, j), id(i, j + 1), id(i + 1, j + 1), id(i + 1, j))) for j in 1:nz-1 for i in 1:nn-1]
    return Grid(cells, nodes), id
end

function shell_dofhandler(grid)
    ip = Lagrange{RefQuadrilateral, 1}()
    dh = DofHandler(grid)
    add!(dh, :u, ip^3)
    add!(dh, :θ, ip^3)
    close!(dh)
    node_to_dofs = Dict{Int, Vector{Int}}()
    for cell in CellIterator(dh)
        cd = celldofs(cell)
        for (i, n) in enumerate(cell.nodes)
            node_to_dofs[n] = [cd[3(i - 1) + 1], cd[3(i - 1) + 2], cd[3(i - 1) + 3],
                               cd[12 + 3(i - 1) + 1], cd[12 + 3(i - 1) + 2], cd[12 + 3(i - 1) + 3]]
        end
    end
    return dh, node_to_dofs
end

# least-squares rigid rotation (about the section's mean point) of the in-plane displacements
function ls_rotation(X, Y, ux, uy)
    xm = mean(X); ym = mean(Y)
    return sum((X .- xm) .* uy .- (Y .- ym) .* ux) / sum((X .- xm) .^ 2 .+ (Y .- ym) .^ 2)
end

# straight-line least-squares fit y = a + b z, with R² and max deviation
function linear_fit(z, y)
    A = [ones(length(z)) z]
    c = A \ y
    yhat = A * c
    return c[1], c[2], 1 - sum((y .- yhat) .^ 2) / sum((y .- mean(y)) .^ 2), maximum(abs, y .- yhat)
end

# ---------------------------------------------------------------------------------------------
# One torsion analysis
# ---------------------------------------------------------------------------------------------
"""
    torsion_J(; n_flat, n_corner, n_z, Cs, rotate_about, L, drilling, drill_factor, hb_gamma)

Warping-free static twist analysis of the upright (Moen 2008, Fig. 4.41). Returns a NamedTuple with
the torque, the J estimates from Eq. 4.13 and from the slope of the twist, the linearity check of
β(z), the thin-walled J and the twist profile.
"""
function torsion_J(; n_flat = 4, n_corner = 5, n_z = Int(2 * shape.length), Cs = Q.DEFAULT_SHEAR_RELAXATION,
                     rotate_about = :shear_center, L = shape.length,
                     drilling = :hughes_brezzi, drilling_gamma = 1.0, verbose = true)

    X, Y = centerline(shape; n_flat, n_corner)
    Z = collect(range(0.0, L, n_z + 1))
    nn = length(X)
    nz = length(Z)

    sp = section_properties(X, Y, shape.t)
    J_thin = sp.J                                        # Σ b t³/3 on the discretized centerline
    Xc, Yc = rotate_about == :shear_center ? (sp.xs, sp.ys) : (sp.xc, sp.yc)

    grid, id = extruded_grid(X, Y, Z)
    dh, node_to_dofs = shell_dofhandler(grid)

    # node sets
    end0 = Set(id(i, 1) for i in 1:nn)
    endL = Set(id(i, nz) for i in 1:nn)
    i_webmid = argmin(X .^ 2 .+ (Y .- shape.D / 2) .^ 2)                # web centerline is at X = 0
    web_nodes = findall(x -> abs(x) < 1e-6, X)
    i_web_top = web_nodes[argmax(Y[web_nodes])]                        # web/flange corner, top
    i_web_bot = web_nodes[argmin(Y[web_nodes])]                        # web/flange corner, bottom

    # stiffness
    ip4 = IP4(); ip6 = IP6()
    qr_m = QuadratureRule{RefQuadrilateral}(2)
    qr_b = QuadratureRule{RefQuadrilateral}(3)
    K = allocate_matrix(dh)
    K = Q.assemble_global_Ke!(K, dh, qr_m, qr_b, ip4, ip6, E, ν, shape.t; Cs, drilling, drilling_gamma)
    K0 = copy(K)                                                      # unconstrained copy for reactions

    # boundary conditions, Figure 4.41
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, end0, (x, t) -> [0.0, 0.0], [1, 2]))                            # twist fixed, warping free
    add!(ch, Dirichlet(:u, Set([id(i_webmid, 1)]), (x, t) -> [0.0], [3]))                   # axial rigid body mode
    add!(ch, Dirichlet(:u, endL, (x, t) -> [-βo * (x[2] - Yc), βo * (x[1] - Xc)], [1, 2]))  # rigid twist βo about (Xc, Yc)
    close!(ch)
    update!(ch, 0.0)

    f = zeros(ndofs(dh))
    apply!(K, f, ch)
    u = K \ f
    apply!(u, ch)

    # reactions and torque at the twisted end (and at the fixed end, for equilibrium)
    R = K0 * u
    T = 0.0; Fx = 0.0; Fy = 0.0; T0 = 0.0
    for n in endL
        x, y = grid.nodes[n].x[1], grid.nodes[n].x[2]
        Rx, Ry = R[node_to_dofs[n][1]], R[node_to_dofs[n][2]]
        T += (x - Xc) * Ry - (y - Yc) * Rx
        Fx += Rx; Fy += Ry
    end
    for n in end0
        x, y = grid.nodes[n].x[1], grid.nodes[n].x[2]
        T0 += (x - Xc) * R[node_to_dofs[n][2]] - (y - Yc) * R[node_to_dofs[n][1]]
    end

    # twist along the member
    β_ls = zeros(nz); β_web = zeros(nz)
    for j in 1:nz
        ux = [u[node_to_dofs[id(i, j)][1]] for i in 1:nn]
        uy = [u[node_to_dofs[id(i, j)][2]] for i in 1:nn]
        β_ls[j] = ls_rotation(X, Y, ux, uy)
        β_web[j] = -(ux[i_web_top] - ux[i_web_bot]) / (Y[i_web_top] - Y[i_web_bot])
    end
    _, slope_ls, R2_ls, dev_ls = linear_fit(Z, β_ls)
    _, slope_web, R2_web, dev_web = linear_fit(Z, β_web)

    J_end = T * L / (G * βo)                 # Eq. 4.13
    J_slope_ls = T / (G * slope_ls)          # T = G J dβ/dz with the fitted slope
    J_slope_web = T / (G * slope_web)

    if verbose
        @printf("  mesh: %d nodes/section, %d stations, %d elements, %d dofs; Cs = %.2f; drilling = %s; rotation about %s (%.4f, %.4f)\n",
                nn, nz, getncells(grid), ndofs(dh), Cs, String(drilling), String(rotate_about), Xc, Yc)
        @printf("  T = %.6e lbf-in at z = L,  T at z = 0 = %.6e,  ΣFx = %.2e, ΣFy = %.2e lbf\n", T, T0, Fx, Fy)
        @printf("  β(L)/βo = %.6f (ls) %.6f (web corners);  linear fit R² = %.8f (ls) %.8f (web);  max |β − fit|/βo = %.1e (ls) %.1e (web)\n",
                β_ls[end] / βo, β_web[end] / βo, R2_ls, R2_web, dev_ls / βo, dev_web / βo)
        @printf("  J = T L / (G βo)             = %.6e in⁴   (%.4f × Σbt³/3)\n", J_end, J_end / J_thin)
        @printf("  J = T / (G dβ/dz), ls twist  = %.6e in⁴   (%.4f × Σbt³/3)\n", J_slope_ls, J_slope_ls / J_thin)
        @printf("  J = T / (G dβ/dz), web twist = %.6e in⁴   (%.4f × Σbt³/3)\n", J_slope_web, J_slope_web / J_thin)
    end

    return (; n_flat, n_corner, n_z, Cs, rotate_about, L, drilling, drilling_gamma, nn, nz, ncells = getncells(grid), ndofs = ndofs(dh),
              T, T0, Fx, Fy, J_end, J_slope_ls, J_slope_web, J_thin, slope_ls, slope_web,
              R2_ls, R2_web, dev_ls, dev_web, Z, β_ls, β_web, Xc, Yc, sp)
end

# ---------------------------------------------------------------------------------------------
# Run the study
# ---------------------------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__

    println("Single upright: t = $(shape.t), lip = $(shape.lip), D = $(shape.D), B = $(shape.B), R = $(shape.R) (centerline), L = $(shape.length) in")
    X, Y = centerline(shape)
    sp = section_properties(X, Y, shape.t)
    @printf("Thin-walled (cutwp): A = %.4f in², J = Σbt³/3 = %.6e in⁴, Cw = %.5f in⁶, centroid (%.4f, %.4f), shear center (%.4f, %.4f)\n",
            sp.A, sp.J, sp.Cw, sp.xc, sp.yc, sp.xs, sp.ys)

    # exact Saint-Venant J of the real section (2D Prandtl stress function, quadratic elements)
    left, right, bcl = upright_section(shape; n_flat = 160, n_corner = 48)
    g2, bd2 = strip_grid(left, right; n_t = 8)
    J_SV, A_SV, nd_SV = saint_venant_J(g2, bd2; order = 2)
    @printf("Saint-Venant 2D reference:   A = %.4f in², J_SV = %.6e in⁴  (%.4f × Σbt³/3; %d dofs)\n\n", A_SV, J_SV, J_SV / (bcl * shape.t^3 / 3), nd_SV)

    cases = [
        (label = "Abaqus mesh (4/flat, 5/corner, 222 along), Hughes–Brezzi drilling", kw = (;)),
        (label = "  same, Cs = 0 (unrelaxed transverse shear)",                        kw = (; Cs = 0.0)),
        (label = "  same, rotation applied about the centroid",                        kw = (; rotate_about = :centroid)),
        (label = "  same, half length L = 55.5 in",                                    kw = (; L = shape.length / 2, n_z = Int(shape.length))),
        (label = "Coarse mesh (2/flat, 2/corner, 111 along), Hughes–Brezzi",           kw = (; n_flat = 2, n_corner = 2, n_z = Int(shape.length))),
        (label = "Fine mesh (8/flat, 10/corner, 444 along), Hughes–Brezzi",            kw = (; n_flat = 8, n_corner = 10, n_z = Int(4 * shape.length))),
        (label = "Abaqus mesh, original absolute drilling penalty (1/100)",            kw = (; drilling = :penalty)),
        (label = "Fine mesh, original absolute drilling penalty (1/100)",              kw = (; n_flat = 8, n_corner = 10, n_z = Int(4 * shape.length), drilling = :penalty)),
    ]

    results = []
    for c in cases
        println(c.label)
        t_run = @elapsed r = torsion_J(; c.kw..., verbose = true)
        @printf("  (%.1f s)\n\n", t_run)
        push!(results, (label = c.label, r = r))
    end

    println("Summary   (J_SV = $(@sprintf("%.6e", J_SV)) in⁴ exact Saint-Venant, Σbt³/3 = $(@sprintf("%.6e", sp.J)) in⁴)")
    @printf("%-72s %8s %13s %13s %9s %9s %11s\n", "case", "dofs", "J_end", "J_slope", "J/J_SV", "J/Jthin", "R² β(z)")
    for (label, r) in results
        @printf("%-72s %8d %13.6e %13.6e %9.4f %9.4f %11.8f\n", label, r.ndofs, r.J_end, r.J_slope_ls, r.J_end / J_SV, r.J_end / r.J_thin, r.R2_ls)
    end

    # write results
    open(joinpath(@__DIR__, "J_results.csv"), "w") do io
        println(io, "case,drilling,n_flat,n_corner,n_z,Cs,rotate_about,L,ndofs,T_lbf_in,J_end_in4,J_slope_ls_in4,J_slope_web_in4,J_thin_walled_in4,J_saint_venant_2D_in4,J_end_over_J_SV,J_end_over_J_thin,R2_ls,R2_web,max_dev_ls_over_betao,max_dev_web_over_betao,sumFx,sumFy")
        for (label, r) in results
            println(io, join(["\"" * strip(label) * "\"", r.drilling, r.n_flat, r.n_corner, r.n_z, r.Cs, r.rotate_about, r.L, r.ndofs, r.T, r.J_end, r.J_slope_ls, r.J_slope_web,
                              r.J_thin, J_SV, r.J_end / J_SV, r.J_end / r.J_thin, r.R2_ls, r.R2_web, r.dev_ls / βo, r.dev_web / βo, r.Fx, r.Fy], ","))
        end
    end
    r1 = results[1].r
    writedlm(joinpath(@__DIR__, "twist_profile_abaqus_mesh.csv"),
             ["z_in" "beta_ls_over_betao" "beta_web_over_betao"; r1.Z r1.β_ls ./ βo r1.β_web ./ βo], ',')
    println("\nWrote J_results.csv and twist_profile_abaqus_mesh.csv")
end
