# buckling_built_up.jl
#
# Elastic eigenbuckling of the two-C welded upright (Ferrite.jl + QuadShellFiniteElement.jl) under uniform
# axial compression, pinned and warping-free at both ends.
#
# Geometry and welds: the same two 3 × 3 × 0.074 in lipped C's as calculate_J_built_up_member.jl, member
# length L (default 44 in), 3 in welds at 18 in spacing placed symmetrically about mid-length (for L = 44:
# weld centers at z = 4, 22, 40 in), each weld a rigid body (Ferrite affine constraints) as in the Abaqus
# model.
#
# Boundary conditions ("pinned, warping free"): at z = 0 and z = L every node of both C's is fixed in X
# and Y (translation and twist restrained, rotations about X and Y free, warping free). Axial rigid body
# motion is removed at one mid-length node. The reference load is a self-equilibrated uniform compressive
# stress: tributary nodal forces −P/A at z = L and +P/A at z = 0, P_ref = 1 kip.
#
# Eigenproblem: K φ = P_cr (−K_g) φ on the constrained dofs, K_c = Cᵀ K C, K_g,c = Cᵀ K_g C with C the
# Ferrite constraint matrix (Dirichlet + affine), solved by shift-invert Arnoldi (ArnoldiMethod.jl) with a
# sparse LU of K_c.
#
# Two analyses:
#   :all     — the shell model as is: the lowest modes are local plate buckling of the 2.9 in flanges/webs.
#   :global  — the cross section of the whole built-up member is constrained to move rigidly in-plane at
#              every station (u, v, θz shared by all nodes of both C's, warping u_z left free), which
#              suppresses local and distortional buckling and leaves the global flexural and
#              flexural-torsional modes. The weld ties keep only their u_z, θx, θy parts (their in-plane
#              parts are implied by the rigid-section constraint). This is the shell-model counterpart of a
#              beam flexural-torsional buckling analysis, with warping and shear deformation from the FE.
#
# Run:  julia --project=. buckling_built_up.jl          (writes buckling_results.csv, mode_*.csv)

include(joinpath(@__DIR__, "calculate_J_built_up_member.jl"))
using ArnoldiMethod, LinearMaps, SparseArrays

const P_ref = 1000.0     # lbf (1 kip) total reference compression

# rigid tie of the axial / out-of-plane dofs only: u_z, θx, θy (in-plane handled by the section constraint)
function add_rigid_tie_axial!(ch, node_to_dofs, grid, master, slaves)
    xm = grid.nodes[master].x
    dm = node_to_dofs[master]
    for s in slaves
        s == master && continue
        r = grid.nodes[s].x - xm
        ds = node_to_dofs[s]
        add!(ch, AffineConstraint(ds[3], [dm[3] => 1.0, dm[4] => r[2], dm[5] => -r[1]], 0.0))
        add!(ch, AffineConstraint(ds[4], [dm[4] => 1.0], 0.0))
        add!(ch, AffineConstraint(ds[5], [dm[5] => 1.0], 0.0))
    end
end

# symmetric weld layout: welds of length wl at spacing ws, centered about L/2, kept ≥ wl/2 from the ends
function symmetric_weld_locations(L, wl, ws)
    n = floor(Int, (L - wl) / ws) + 1
    return collect(range(L / 2 - (n - 1) * ws / 2, L / 2 + (n - 1) * ws / 2, n))
end

"""
    built_up_buckling(; L, weld_length, weld_spacing, mode = :all | :global, nev, n_flat, n_corner, n_z)

Returns critical loads P_cr (lbf, ascending), the corresponding full-dof mode vectors, the grid, node dof
map, and per-mode classification data.
"""
function built_up_buckling(; L = 44.0, weld_length = 3.0, weld_spacing = 18.0, mode = :all, nev = 6,
                             n_flat = 4, n_corner = 5, n_z = Int(round(2L)), Cs = Q.DEFAULT_SHEAR_RELAXATION, verbose = true)

    X1, Y1 = centerline(shape; n_flat, n_corner)
    Z = collect(range(0.0, L, n_z + 1))
    nn = length(X1); nz = length(Z)
    grid, id, Xs, Ys = built_up_grid(X1, Y1, shape.B, Z)
    dh, nd = shell_dofhandler(grid)
    R = shape.R; B = shape.B; D = shape.D; t = shape.t
    near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
    c1_top = findall(i -> near(X1[i], B, R) && near(Y1[i], D, R), 1:nn)
    c1_bot = findall(i -> near(X1[i], B, R) && near(Y1[i], 0.0, R), 1:nn)
    c2_top = findall(i -> near(X1[i], 0.0, R) && near(Y1[i], D, R), 1:nn)
    c2_bot = findall(i -> near(X1[i], 0.0, R) && near(Y1[i], 0.0, R), 1:nn)
    i_webmid = argmin(X1 .^ 2 .+ (Y1 .- D / 2) .^ 2)
    j_mid = argmin(abs.(Z .- L / 2))

    # ---- constraints
    ch = ConstraintHandler(dh)
    end0 = Set(id(s, i, 1) for s in 1:2 for i in 1:nn)
    endL = Set(id(s, i, nz) for s in 1:2 for i in 1:nn)
    add!(ch, Dirichlet(:u, union(end0, endL), (x, t) -> [0.0, 0.0], [1, 2]))          # pinned, warping free
    add!(ch, Dirichlet(:u, Set([id(1, i_webmid, j_mid)]), (x, t) -> [0.0], [3]))         # axial datum (self-equilibrated load)

    weld_locations = symmetric_weld_locations(L, weld_length, weld_spacing)
    interior = 2:nz-1
    claimed = falses(nz)
    weld_stations = [begin
                         st = [j for j in interior if abs(Z[j] - zw) <= weld_length / 2 + 1e-9 && !claimed[j]]
                         claimed[st] .= true
                         st
                     end for zw in weld_locations]
    for st in weld_stations, (c1, c2) in ((c1_top, c2_top), (c1_bot, c2_bot))
        slaves = vcat([id(1, i, j) for j in st for i in c1], [id(2, i, j) for j in st for i in c2])
        if mode == :global
            add_rigid_tie_axial!(ch, nd, grid, slaves[1], slaves)
        else
            add_rigid_tie!(ch, nd, grid, slaves[1], slaves)
        end
    end
    if mode == :global
        # rigid in-plane cross section per interior station: u_x,i = u_x,m − θz,m (y_i − y_m), u_y,i = u_y,m + θz,m (x_i − x_m)
        for j in interior
            m = id(1, i_webmid, j); xm = grid.nodes[m].x; dm = nd[m]
            for s in 1:2, i in 1:nn
                n = id(s, i, j); n == m && continue
                r = grid.nodes[n].x - xm; dn = nd[n]
                add!(ch, AffineConstraint(dn[1], [dm[1] => 1.0, dm[6] => -r[2]], 0.0))
                add!(ch, AffineConstraint(dn[2], [dm[2] => 1.0, dm[6] => r[1]], 0.0))
            end
        end
    end
    close!(ch); update!(ch, 0.0)

    # ---- stiffness and reference load (uniform compression, tributary nodal forces at both ends)
    K = allocate_matrix(dh, ch)
    K = Q.assemble_global_Ke!(K, dh, QuadratureRule{RefQuadrilateral}(2), QuadratureRule{RefQuadrilateral}(3), IP4(), IP6(), E, ν, t; Cs)
    seg = [hypot(X1[i+1] - X1[i], Y1[i+1] - Y1[i]) for i in 1:nn-1]
    trib = [(i > 1 ? seg[i-1] / 2 : 0.0) + (i < nn ? seg[i] / 2 : 0.0) for i in 1:nn]     # tributary centerline length per node
    A_total = 2 * sum(seg) * t
    f = zeros(ndofs(dh))
    for s in 1:2, i in 1:nn
        f[nd[id(s, i, nz)][3]] -= P_ref / A_total * trib[i] * t
        f[nd[id(s, i, 1)][3]]  += P_ref / A_total * trib[i] * t
    end
    K0 = copy(K)
    apply!(K, f, ch)
    u = K \ f
    apply!(u, ch)

    # ---- geometric stiffness from the prebuckling membrane stresses
    qr_g = QuadratureRule{RefQuadrilateral}(2)
    σXX, σYY, τXY = Q.element_membrane_stresses(dh, u, IP4(), E, ν, t; qr = qr_g)
    σz_mean = mean(mean.(σXX))                     # element local x runs along the extrusion direction z
    Kg = allocate_matrix(dh)
    Kg = Q.assemble_global_Kg!(Kg, dh, qr_g, IP4(), σXX .* t, σYY .* t, τXY .* t)

    # ---- condensed eigenproblem  K_c φ = P (−K_g,c) φ  via shift-invert Arnoldi on  K_c⁻¹ (−K_g,c)
    C, _ = Ferrite.create_constraint_matrix(ch)
    Kc = Symmetric(sparse(C' * K0 * C))
    Kgc = sparse(C' * (-Kg) * C)
    F = lu(sparse(Kc))
    op = LinearMap{Float64}(x -> F \ (Kgc * x), size(Kc, 1); ismutating = false)
    decomp, hist = partialschur(op; nev = nev + 4, tol = 1e-8, which = isdefined(ArnoldiMethod, :LR) ? ArnoldiMethod.LR() : :LR)
    μ, Ψ = partialeigen(decomp)
    keep = findall(m -> real(m) > 1e-12, μ)
    order = sortperm(real.(μ[keep]); rev = true)[1:min(nev, length(keep))]
    P_cr = [P_ref / real(μ[keep[k]]) for k in order]
    modes = [C * real.(Ψ[:, keep[k]]) for k in order]

    # ---- classification: rigid-section participation of the in-plane displacements, and mid-length section motion
    info = map(modes) do a
        rig = 0.0; tot = 0.0
        for j in 2:nz-1
            xs = vcat([Xs[s][i] for s in 1:2 for i in 1:nn]); ys = vcat([Ys[s][i] for s in 1:2 for i in 1:nn])
            ux = [a[nd[id(s, i, j)][1]] for s in 1:2 for i in 1:nn]; uy = [a[nd[id(s, i, j)][2]] for s in 1:2 for i in 1:nn]
            xm = mean(xs); ym = mean(ys)
            U = mean(ux); V = mean(uy)
            θ = ls_rotation(xs, ys, ux, uy)
            uxr = U .- θ .* (ys .- ym); uyr = V .+ θ .* (xs .- xm)
            rig += sum(uxr .^ 2 .+ uyr .^ 2); tot += sum(ux .^ 2 .+ uy .^ 2)
        end
        xs = vcat([Xs[s][i] for s in 1:2 for i in 1:nn]); ys = vcat([Ys[s][i] for s in 1:2 for i in 1:nn])
        ux = [a[nd[id(s, i, j_mid)][1]] for s in 1:2 for i in 1:nn]; uy = [a[nd[id(s, i, j_mid)][2]] for s in 1:2 for i in 1:nn]
        U = mean(ux); V = mean(uy); θ = ls_rotation(xs, ys, ux, uy)
        rmax = maximum(hypot.(xs .- mean(xs), ys .- mean(ys)))
        (; participation = rig / tot, U, V, θ, twist_to_translation = abs(θ) * rmax / max(hypot(U, V), 1e-30))
    end

    if verbose
        @printf("  %s: L = %.1f in, welds at %s, %d dofs, %d constrained; A = %.4f in², σ_ref = %.2f psi (uniform check: mean σ_z = %.2f)\n",
                mode, L, string(weld_locations), ndofs(dh), length(ch.prescribed_dofs), A_total, P_ref / A_total, -σz_mean)
        for (k, P) in enumerate(P_cr)
            c = info[k]
            @printf("  mode %d: P_cr = %10.1f lbf = %8.2f kips, σ_cr = %8.2f ksi;  rigid-section participation %.3f;  mid-length U = %+.3f V = %+.3f θ·r_max/|U,V| = %.2f\n",
                    k, P, P / 1000, P / A_total / 1000, c.participation, c.U / maximum(abs, [c.U, c.V, 1e-30]), c.V / maximum(abs, [c.U, c.V, 1e-30]), c.twist_to_translation)
        end
    end
    return (; P_cr, modes, info, grid, nd, id, nn, nz, Xs, Ys, Z, A_total, weld_locations, weld_stations, ndofs = ndofs(dh))
end

# write a mode shape (node coordinates + displacements) for plotting
function write_mode(path, r, k)
    a = r.modes[k]
    open(path, "w") do io
        println(io, "shape,i,j,x,y,z,ux,uy,uz")
        for j in 1:r.nz, s in 1:2, i in 1:r.nn
            n = r.id(s, i, j); p = r.grid.nodes[n].x; d = r.nd[n]
            println(io, join([s, i, j, p[1], p[2], p[3], a[d[1]], a[d[2]], a[d[3]]], ","))
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    L = 44.0
    println("Built-up two-C upright, L = $L in, pinned warping-free ends, 3 in welds at 18 in spacing\n")
    println("Global modes (rigid in-plane cross section, warping free):")
    g = built_up_buckling(; L, mode = :global, nev = 4)
    println("\nAll modes (unconstrained shell):")
    a = built_up_buckling(; L, mode = :all, nev = 4)

    # beam-theory context: Euler loads of the fully composite section (two C's rigidly joined)
    X1, Y1 = centerline(shape); sp = section_properties(X1, Y1, shape.t)
    A2 = 2sp.A; xc = sp.xc + shape.B / 2
    Ix2 = 2sp.Ixx; Iy2 = 2 * (sp.Iyy + sp.A * (shape.B / 2)^2)
    Pey = π^2 * E * Ix2 / L^2; Pex = π^2 * E * Iy2 / L^2
    @printf("\nComposite section (rigidly joined): A = %.4f in², Ix = %.4f in⁴ (about horizontal axis), Iy = %.4f in⁴;  Euler P_ey = %.1f kips (bending about x, deflection in y), P_ex = %.1f kips\n",
            A2, Ix2, Iy2, Pey / 1000, Pex / 1000)

    open(joinpath(@__DIR__, "buckling_results.csv"), "w") do io
        println(io, "analysis,mode,P_cr_lbf,P_cr_kips,sigma_cr_ksi,rigid_section_participation,U_mid,V_mid,theta_mid,twist_to_translation")
        for (name, r) in (("global", g), ("all", a)), (k, P) in enumerate(r.P_cr)
            c = r.info[k]
            println(io, join([name, k, P, P / 1000, P / r.A_total / 1000, c.participation, c.U, c.V, c.θ, c.twist_to_translation], ","))
        end
        println(io, "euler_composite,P_ey,$Pey,$(Pey/1000),$(Pey/A2/1000),,,,,")
        println(io, "euler_composite,P_ex,$Pex,$(Pex/1000),$(Pex/A2/1000),,,,,")
    end
    write_mode(joinpath(@__DIR__, "mode_global_1.csv"), g, 1)
    write_mode(joinpath(@__DIR__, "mode_global_2.csv"), g, 2)
    write_mode(joinpath(@__DIR__, "mode_all_1.csv"), a, 1)
    println("\nWrote buckling_results.csv, mode_global_1.csv, mode_global_2.csv, mode_all_1.csv")
end
