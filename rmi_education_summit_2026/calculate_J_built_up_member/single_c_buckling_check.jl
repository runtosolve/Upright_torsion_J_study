# single_c_buckling_check.jl — validation of the rigid-section torsional / flexural-torsional buckling analysis on a
# single lipped C, whose J, C_w, shear center and r_o are known exactly, against the classical formulas
#     P_t  = (G J + π² E C_w / L²) / (A r_o²)          P_ey = π² E I_x / L²
#     P_FT = [(P_ey + P_t) − √((P_ey + P_t)² − 4 β P_ey P_t)] / (2β),   β = 1 − (x_o / r_o)²
# Same modelling as buckling_built_up.jl: pinned warping-free ends, rigid in-plane section per station, K_g from
# the uniform axial reference stress.   Run:  julia --project=. single_c_buckling_check.jl
include(joinpath(@__DIR__, "buckling_built_up.jl"))

function single_c_buckling(; L, restrain_translations = false, restrain = Symbol[], nev = 3, n_flat = 4, n_corner = 5, nz_per_in = 2, Cs = Q.DEFAULT_SHEAR_RELAXATION)
    X, Y = centerline(shape; n_flat, n_corner)
    Z = collect(range(0.0, L, Int(round(nz_per_in * L)) + 1))
    nn = length(X); nz = length(Z)
    grid, id = extruded_grid(X, Y, Z)
    dh, nd = shell_dofhandler(grid)
    sp = section_properties(X, Y, shape.t)
    xs, ys = sp.xs, sp.ys
    i_webmid = argmin(X .^ 2 .+ (Y .- shape.D / 2) .^ 2); j_mid = argmin(abs.(Z .- L / 2))
    ch = ConstraintHandler(dh)
    ends = Set(vcat([id(i, 1) for i in 1:nn], [id(i, nz) for i in 1:nn]))
    add!(ch, Dirichlet(:u, ends, (x, t) -> [0.0, 0.0], [1, 2]))
    add!(ch, Dirichlet(:u, Set([id(i_webmid, j_mid)]), (x, t) -> [0.0], [3]))
    ru = :u in restrain; rv = :v in restrain; rθ = :θ in restrain
    interior_nodes = Set(id(i, j) for j in 2:nz-1 for i in 1:nn)
    ru && add!(ch, Dirichlet(:u, interior_nodes, (x, t) -> [0.0], [1]))
    rv && add!(ch, Dirichlet(:u, interior_nodes, (x, t) -> [0.0], [2]))
    for j in 2:nz-1
        m = id(i_webmid, j); xm = grid.nodes[m].x; dm = nd[m]
        rθ && add!(ch, Dirichlet(:θ, Set([m]), (x, t) -> [0.0], [3]))
        for i in 1:nn
            n = id(i, j); p = grid.nodes[n].x; dn = nd[n]
            if restrain_translations                        # pure torsion about the shear center
                add!(ch, AffineConstraint(dn[1], [dm[6] => -(p[2] - ys)], 0.0))
                add!(ch, AffineConstraint(dn[2], [dm[6] => (p[1] - xs)], 0.0))
            elseif n != m                                   # free rigid-section motion (u, v, θ) minus restrained parts
                if !ru
                    terms = Pair{Int,Float64}[dm[1] => 1.0]; rθ || push!(terms, dm[6] => -(p[2] - xm[2]))
                    add!(ch, AffineConstraint(dn[1], terms, 0.0))
                end
                if !rv
                    terms = Pair{Int,Float64}[dm[2] => 1.0]; rθ || push!(terms, dm[6] => (p[1] - xm[1]))
                    add!(ch, AffineConstraint(dn[2], terms, 0.0))
                end
            end
        end
    end
    close!(ch); update!(ch, 0.0)
    K = allocate_matrix(dh, ch)
    K = Q.assemble_global_Ke!(K, dh, QuadratureRule{RefQuadrilateral}(2), QuadratureRule{RefQuadrilateral}(3), IP4(), IP6(), E, ν, shape.t; Cs)
    seg = [hypot(X[i+1] - X[i], Y[i+1] - Y[i]) for i in 1:nn-1]; A = sum(seg) * shape.t
    qr_g = QuadratureRule{RefQuadrilateral}(2); ncq = getnquadpoints(CellValues(qr_g, IP4(), IP4()))
    σ0 = -P_ref / A
    σXX = [fill(σ0, ncq) for _ in 1:getncells(grid)]; σYY = [zeros(ncq) for _ in 1:getncells(grid)]; τXY = [zeros(ncq) for _ in 1:getncells(grid)]
    Kg = allocate_matrix(dh); Kg = Q.assemble_global_Kg!(Kg, dh, qr_g, IP4(), σXX .* shape.t, σYY .* shape.t, τXY .* shape.t)
    C, _ = Ferrite.create_constraint_matrix(ch)
    Kc = sparse(C' * K * C); Kgc = sparse(C' * (-Kg) * C)
    F = lu(Kc)
    op = LinearMap{Float64}(x -> F \ (Kgc * x), size(Kc, 1); ismutating = false)
    decomp, _ = partialschur(op; nev = nev + 3, tol = 1e-8, which = ArnoldiMethod.LR())
    μ, _ = partialeigen(decomp)
    keep = sort(real.(μ[real.(μ) .> 1e-12]); rev = true)
    return [P_ref / m for m in keep[1:min(nev, length(keep))]], sp, A
end

J_SV = 1.334851e-3
L = 44.0
X, Y = centerline(shape); sp = section_properties(X, Y, shape.t)
x_o = sp.xs - sp.xc; r_o2 = (sp.Ixx + sp.Iyy) / sp.A + x_o^2
Pt_th = (G * J_SV + π^2 * E * sp.Cw / L^2) / r_o2; Pey = π^2 * E * sp.Ixx / L^2; Pex = π^2 * E * sp.Iyy / L^2   # torsional buckling load (force)
@printf("theory at L = 44: P_t = %.2f, P_ey = %.2f, P_ex = %.2f kips (Cw = %.4f in⁶, r_o = %.3f in)\n", Pt_th/1e3, Pey/1e3, Pex/1e3, sp.Cw, sqrt(r_o2))
for (lab, kw) in (("pure torsion, 4/5 mesh, 0.5 in",            (; restrain_translations = true)),
                  ("pure torsion, 8/10 mesh, 0.5 in",           (; restrain_translations = true, n_flat = 8, n_corner = 10)),
                  ("pure torsion, 8/10 mesh, 0.25 in",          (; restrain_translations = true, n_flat = 8, n_corner = 10, nz_per_in = 4)),
                  ("pure torsion, 4/5 mesh, Cs = 0",            (; restrain_translations = true, Cs = 0.0)),
                  ("pure y-flexure (u = θ = 0), 4/5 mesh",      (; restrain = [:u, :θ])),
                  ("pure x-flexure (v = θ = 0), 4/5 mesh",      (; restrain = [:v, :θ])),
                  ("pure x-flexure (v = θ = 0), 8/10, 0.25 in", (; restrain = [:v, :θ], n_flat = 8, n_corner = 10, nz_per_in = 4)))
    P, _, _ = single_c_buckling(; L, nev = 1, kw...)
    ref = occursin("torsion", lab) ? Pt_th : occursin("y-flexure", lab) ? Pey : Pex
    @printf("  %-44s shell %8.2f kips   theory %8.2f   ratio %.3f\n", lab, P[1]/1e3, ref/1e3, P[1]/ref)
end
