# calculate_J_built_up_member.jl
#
# Effective Saint-Venant torsion constant of the OneRack built-up column of combo_1/J/combo_1_J_r5.jl
# (Abaqus model combo_1_J_r5.inp) computed with a Ferrite.jl shell model (QuadShellFiniteElement.jl,
# Hughes–Brezzi drilling) by the static twist method of Moen (2008), Section 4.2.7.3.2.3:
#
#     J_eff = T L / (G βo)                                                   (Eq. 4.13)
#
# Layout (identical to the Abaqus assembly): two 3 × 3 × 0.074 in lipped C uprights, SHAPE_2 = SHAPE_1
# translated by B = 3.0 in in X, so the lip/flange corners of C1 bear on the web of C2. Seven 3 in long
# welds at z = 1.5, 19.5, …, 109.5 in tie the C1 top and bottom lip-corner nodes to the C2 web-corner
# nodes (nodes within R of the corner, |Δz| ≤ 1.5 in) as rigid bodies — the Abaqus *Rigid Body, tie
# definitions — implemented here as Ferrite affine constraints to one master node per weld.
#
# Boundary conditions (Fig. 4.41 / the Abaqus step): z = 0 all nodes of both shapes fixed in X and Y,
# C2 web mid-height node fixed in Z; z = L all nodes rotated rigidly by βo about the Abaqus reference
# point (B/2, D/2) = (1.5, 1.5) (the kinematic coupling to node 10000 in dofs 1, 2); warping free.
#
# Because the welds are intermittent the twist is not uniform: between welds the two open C sections
# twist as short warping-restrained members, at the welds the C1 web + flanges + C2 web form a closed
# cell. J_eff is therefore an equivalent member property (it contains warping stiffness), which is what
# the dissertation method and calculate_J_r5.jl report for Abaqus.
#
# Cases: (1) welds as in Abaqus; (2) fully welded (back-to-back 3 in welds) → closed-cell check against Bredt;
# (3) no welds → 2 × single C check; (4) welds + bilateral normal tie of the C1 lips to the C2 web
# between welds (a linear stand-in for the frictionless contact in the Abaqus model).
#
# Run from this folder:   julia --project=. calculate_J_built_up_member.jl

include(joinpath(@__DIR__, "..", "calculate_J_single_member", "calculate_J_single_member.jl"))   # shape, E, ν, G, βo, centerline, shell_dofhandler, ls_rotation, linear_fit, saint_venant_J …

const weld = (length = 3.0, locations = collect(range(0.0, shape.length - 3.0, 7) .+ 3.0 / 2))   # as in combo_1_J_r5.jl
const twist_center = (shape.B / 2, shape.D / 2)                                                  # Abaqus node 10000

# ---------------------------------------------------------------------------------------------
# Two-shape extruded mesh. Node id(s, i, j): shape s ∈ (1, 2), section node i, station j
# ---------------------------------------------------------------------------------------------
function built_up_grid(X1, Y1, ΔX, Z)
    nn = length(X1); nz = length(Z)
    Xs = (X1, X1 .+ ΔX); Ys = (Y1, Y1)
    nodes = [Node(Vec((Xs[s][i], Ys[s][i], Z[j]))) for j in 1:nz for s in 1:2 for i in 1:nn]
    id(s, i, j) = (j - 1) * 2nn + (s - 1) * nn + i
    cells = [Quadrilateral((id(s, i, j), id(s, i, j + 1), id(s, i + 1, j + 1), id(s, i + 1, j))) for j in 1:nz-1 for s in 1:2 for i in 1:nn-1]
    return Grid(cells, nodes), id, Xs, Ys
end

# rigid body tie of `slaves` to `master` (6 dofs each): u_s = u_m + θ_m × (x_s − x_m), θ_s = θ_m
function add_rigid_tie!(ch, node_to_dofs, grid, master, slaves)
    xm = grid.nodes[master].x
    dm = node_to_dofs[master]
    for s in slaves
        s == master && continue
        r = grid.nodes[s].x - xm
        ds = node_to_dofs[s]
        add!(ch, AffineConstraint(ds[1], [dm[1] => 1.0, dm[5] => r[3], dm[6] => -r[2]], 0.0))
        add!(ch, AffineConstraint(ds[2], [dm[2] => 1.0, dm[6] => r[1], dm[4] => -r[3]], 0.0))
        add!(ch, AffineConstraint(ds[3], [dm[3] => 1.0, dm[4] => r[2], dm[5] => -r[1]], 0.0))
        for k in 4:6
            add!(ch, AffineConstraint(ds[k], [dm[k] => 1.0], 0.0))
        end
    end
end

"""
    built_up_J(; welds = :abaqus | :full | :none, contact = false, n_flat, n_corner, n_z, Cs, center)

Static twist of the two-C built-up member. Returns J_eff = T L / (G βo), the end torque, net end
force, and the twist profiles β₁(z), β₂(z) of the two shapes (least-squares rigid rotation).
"""
function built_up_J(; welds = :abaqus, contact = false, n_flat = 4, n_corner = 5, L = shape.length, n_z = Int(round(2L)),
                      weld_length = weld.length, weld_locations = weld.locations,
                      Cs = Q.DEFAULT_SHEAR_RELAXATION, center = twist_center, verbose = true)

    X1, Y1 = centerline(shape; n_flat, n_corner)
    Z = collect(range(0.0, L, n_z + 1))
    nn = length(X1); nz = length(Z)
    grid, id, Xs, Ys = built_up_grid(X1, Y1, shape.B, Z)
    dh, nd = shell_dofhandler(grid)
    R = shape.R; B = shape.B; D = shape.D
    Xc, Yc = center

    # ---- node sets in part coordinates (same criteria as LinesCurvesNodes.find_nodes in the Abaqus script)
    near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
    c1_top = findall(i -> near(X1[i], B, R) && near(Y1[i], D, R), 1:nn)      # C1 top lip/flange corner
    c1_bot = findall(i -> near(X1[i], B, R) && near(Y1[i], 0.0, R), 1:nn)    # C1 bottom lip/flange corner
    c2_top = findall(i -> near(X1[i], 0.0, R) && near(Y1[i], D, R), 1:nn)    # C2 top web/flange corner
    c2_bot = findall(i -> near(X1[i], 0.0, R) && near(Y1[i], 0.0, R), 1:nn)  # C2 bottom web/flange corner
    i_webmid = argmin(X1 .^ 2 .+ (Y1 .- D / 2) .^ 2)

    ch = ConstraintHandler(dh)
    end0 = Set(id(s, i, 1) for s in 1:2 for i in 1:nn)
    endL = Set(id(s, i, nz) for s in 1:2 for i in 1:nn)
    add!(ch, Dirichlet(:u, end0, (x, t) -> [0.0, 0.0], [1, 2]))
    add!(ch, Dirichlet(:u, Set([id(2, i_webmid, 1)]), (x, t) -> [0.0], [3]))
    add!(ch, Dirichlet(:u, endL, (x, t) -> [-βo * (x[2] - Yc), βo * (x[1] - Xc)], [1, 2]))

    # ---- welds: stations with |z − z_w| ≤ 1.5 (end stations carry Dirichlet conditions and are left out)
    interior = 2:nz-1
    weld_stations = if welds == :abaqus
        # stations within ± weld_length/2 of each weld center; a station shared by abutting welds goes to the first
        claimed = falses(nz)
        [begin
             st = [j for j in interior if abs(Z[j] - zw) <= weld_length / 2 + 1e-9 && !claimed[j]]
             claimed[st] .= true
             st
         end for zw in weld_locations]
    elseif welds == :full                       # back-to-back welds over the whole length, each its own rigid body
        [[j for j in interior if abs(Z[j] - zw) <= weld_length / 2 - 1e-9 || Z[j] == zw - weld_length / 2] for zw in (weld_length / 2):weld_length:(L - weld_length / 2)]
    else
        Vector{Int}[]
    end
    n_ties = 0
    welded_station = falses(nz)
    for st in weld_stations, (c1, c2) in ((c1_top, c2_top), (c1_bot, c2_bot))
        slaves = vcat([id(1, i, j) for j in st for i in c1], [id(2, i, j) for j in st for i in c2])
        add_rigid_tie!(ch, nd, grid, slaves[1], slaves)
        n_ties += length(slaves) - 1
        welded_station[st] .= true
    end

    # ---- optional bilateral normal tie of the C1 lips to the C2 web between welds (contact stand-in)
    n_contact = 0
    if contact
        lip = findall(i -> near(X1[i], maximum(X1), 1e-6), 1:nn)                      # C1 lip flat nodes
        web = sort(findall(i -> near(X1[i], 0.0, 1e-6), 1:nn); by = i -> Y1[i])          # C2 web flat nodes
        for j in interior
            welded_station[j] && continue
            for i in lip
                k = findfirst(m -> Y1[web[m]] <= Y1[i] <= Y1[web[m+1]], 1:length(web)-1)
                k === nothing && continue
                a, b = web[k], web[k+1]
                s = (Y1[i] - Y1[a]) / (Y1[b] - Y1[a])
                add!(ch, AffineConstraint(nd[id(1, i, j)][1], [nd[id(2, a, j)][1] => 1 - s, nd[id(2, b, j)][1] => s], 0.0))
                n_contact += 1
            end
        end
    end
    close!(ch); update!(ch, 0.0)

    # sparsity pattern must include the affine-constraint couplings
    K = allocate_matrix(dh, ch)
    K = Q.assemble_global_Ke!(K, dh, QuadratureRule{RefQuadrilateral}(2), QuadratureRule{RefQuadrilateral}(3), IP4(), IP6(), E, ν, shape.t; Cs)
    K0 = copy(K)

    f = zeros(ndofs(dh))
    apply!(K, f, ch)
    u = K \ f
    apply!(u, ch)

    # ---- end torque about the twist center and net force
    Rv = K0 * u
    T = 0.0; Fx = 0.0; Fy = 0.0; T0 = 0.0
    for n in endL
        x, y = grid.nodes[n].x[1], grid.nodes[n].x[2]
        Rx, Ry = Rv[nd[n][1]], Rv[nd[n][2]]
        T += (x - Xc) * Ry - (y - Yc) * Rx; Fx += Rx; Fy += Ry
    end
    for n in end0
        x, y = grid.nodes[n].x[1], grid.nodes[n].x[2]
        T0 += (x - Xc) * Rv[nd[n][2]] - (y - Yc) * Rv[nd[n][1]]
    end
    J_eff = T * L / (G * βo)

    # ---- twist of each shape along the length
    β = zeros(nz, 2)
    for s in 1:2, j in 1:nz
        ux = [u[nd[id(s, i, j)][1]] for i in 1:nn]; uy = [u[nd[id(s, i, j)][2]] for i in 1:nn]
        β[j, s] = ls_rotation(Xs[s], Ys[s], ux, uy)
    end
    _, slope_avg, R2, _ = linear_fit(Z, β[:, 1])

    if verbose
        @printf("  welds = %s (%d × %.1f in), L = %.1f, contact = %s: %d dofs, %d rigid-tie constraints, %d contact constraints\n", welds, length(weld_stations), weld_length, L, contact, ndofs(dh), n_ties, n_contact)
        @printf("  T = %.6e lbf-in (z = 0: %.6e),  ΣFx = %.2e, ΣFy = %.2e lbf  (|ΣF|/T = %.1e)\n", T, T0, Fx, Fy, hypot(Fx, Fy) / abs(T))
        @printf("  β₁(L)/βo = %.4f, β₂(L)/βo = %.4f;  linear-fit R² of β₁(z) = %.5f\n", β[end, 1] / βo, β[end, 2] / βo, R2)
        @printf("  J_eff = T L / (G βo) = %.5e in⁴\n", J_eff)
    end
    return (; welds, contact, L, weld_length, n_welds = length(weld_stations), J_eff, T, T0, Fx, Fy, Z, β, welded_station, ndofs = ndofs(dh), n_ties, n_contact,
              u, grid, nd, id, nn, nz, Xs, Ys, Xc, Yc)
end

# ---------------------------------------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    t = shape.t
    X1, Y1 = centerline(shape)
    sp = section_properties(X1, Y1, t)
    J_C = 1.334851e-3                                   # exact Saint-Venant J of one C (calculate_J_single_member.jl)
    a = shape.B - t; A_m = a^2; s_m = 4a                # closed cell: C1 web, C1 flanges, C2 web (centerline 2.926 in square)
    J_box = 4 * A_m^2 * t / s_m
    L_w = length(weld.locations) * weld.length
    println("Built-up member: two 3 × 3 × 0.074 lipped C, C2 offset by B = $(shape.B) in, $(length(weld.locations)) welds × $(weld.length) in at z = $(weld.locations)")
    @printf("References: 2 J_C = %.4e in⁴;  Bredt closed cell J_box = %.4f in⁴;  arithmetic weld-weighted (calculate_J_r5.jl) = %.4f in⁴;  harmonic (series) = %.4e in⁴\n",
            2J_C, J_box, L_w / shape.length * J_box + (1 - L_w / shape.length) * 2J_C, shape.length / (L_w / J_box + (shape.length - L_w) / (2J_C)))
    @printf("Abaqus combo_1_J_r5 (calculate_J_r5.jl, S4R, welds + friction contact): J = %.4f in⁴\n\n", 0.001 * (2.9688e-3 / 2.5) * 111.0 / (E / 2.6 * 3.9660216e-11))

    cases = [(label = "Abaqus weld layout (7 × 3 in), welds only", kw = (;)),
             (label = "Abaqus weld layout + bilateral normal lip/web tie between welds", kw = (; contact = true)),
             (label = "Fully welded (37 back-to-back 3 in welds)", kw = (; welds = :full)),
             (label = "No welds (two independent C)", kw = (; welds = :none)),
             (label = "Abaqus weld layout, twist applied about C1 shear center", kw = (; center = (sp.xs, sp.ys))),
             (label = "Abaqus weld layout, Cs = 0", kw = (; Cs = 0.0)),
             (label = "Abaqus weld layout, fine mesh (8/flat, 10/corner, 444 along)", kw = (; n_flat = 8, n_corner = 10, n_z = Int(4 * shape.length)))]
    results = []
    for c in cases
        println(c.label)
        t_run = @elapsed r = built_up_J(; c.kw...)
        @printf("  (%.1f s)\n\n", t_run)
        push!(results, (label = c.label, r = r))
    end

    println("Summary")
    @printf("%-72s %8s %14s %12s %12s\n", "case", "dofs", "J_eff (in⁴)", "J/J_box", "|ΣF|/T")
    for (label, r) in results
        @printf("%-72s %8d %14.5e %12.4f %12.1e\n", label, r.ndofs, r.J_eff, r.J_eff / J_box, hypot(r.Fx, r.Fy) / abs(r.T))
    end

    open(joinpath(@__DIR__, "J_built_up_results.csv"), "w") do io
        println(io, "case,welds,contact,ndofs,T_lbf_in,J_eff_in4,J_eff_over_J_box,J_box_in4,two_J_C_in4,J_abaqus_r5_in4,sumFx_over_T")
        for (label, r) in results
            println(io, join(["\"" * label * "\"", r.welds, r.contact, r.ndofs, r.T, r.J_eff, r.J_eff / J_box, J_box, 2J_C, 0.2929, hypot(r.Fx, r.Fy) / abs(r.T)], ","))
        end
    end
    r1 = results[1].r
    writedlm(joinpath(@__DIR__, "twist_profile_built_up.csv"),
             ["z_in" "beta_C1_over_betao" "beta_C2_over_betao" "welded_station"; r1.Z r1.β[:, 1] ./ βo r1.β[:, 2] ./ βo Int.(r1.welded_station)], ',')
    println("\nWrote J_built_up_results.csv and twist_profile_built_up.csv")
end
