include(joinpath(@__DIR__, "buckling_built_up.jl"))
J_SV = 1.334851e-3
function torsion_mode_energy(; L, n_flat = 4, n_corner = 5)
    X, Y = centerline(shape; n_flat, n_corner); Z = collect(range(0.0, L, Int(round(2L)) + 1))
    nn = length(X); nz = length(Z); grid, id = extruded_grid(X, Y, Z); dh, nd = shell_dofhandler(grid)
    sp = section_properties(X, Y, shape.t); xs, ys = sp.xs, sp.ys
    i_webmid = argmin(X .^ 2 .+ (Y .- shape.D / 2) .^ 2); j_mid = argmin(abs.(Z .- L / 2))
    ch = ConstraintHandler(dh)
    ends = Set(vcat([id(i, 1) for i in 1:nn], [id(i, nz) for i in 1:nn]))
    add!(ch, Dirichlet(:u, ends, (x, t) -> [0.0, 0.0], [1, 2]))
    add!(ch, Dirichlet(:u, Set([id(i_webmid, j_mid)]), (x, t) -> [0.0], [3]))
    for j in 2:nz-1
        m = id(i_webmid, j); dm = nd[m]
        for i in 1:nn
            n = id(i, j); p = grid.nodes[n].x; dn = nd[n]
            add!(ch, AffineConstraint(dn[1], [dm[6] => -(p[2] - ys)], 0.0)); add!(ch, AffineConstraint(dn[2], [dm[6] => (p[1] - xs)], 0.0))
        end
    end
    close!(ch); update!(ch, 0.0)
    K = allocate_matrix(dh, ch)
    K = Q.assemble_global_Ke!(K, dh, QuadratureRule{RefQuadrilateral}(2), QuadratureRule{RefQuadrilateral}(3), IP4(), IP6(), E, ν, shape.t)
    seg = [hypot(X[i+1] - X[i], Y[i+1] - Y[i]) for i in 1:nn-1]; A = sum(seg) * shape.t
    qr_g = QuadratureRule{RefQuadrilateral}(2); ncq = getnquadpoints(CellValues(qr_g, IP4(), IP4())); nc = getncells(grid)
    σ0 = -P_ref / A
    Kg = allocate_matrix(dh); Kg = Q.assemble_global_Kg!(Kg, dh, qr_g, IP4(), [fill(σ0 * shape.t, ncq) for _ in 1:nc], [zeros(ncq) for _ in 1:nc], [zeros(ncq) for _ in 1:nc])
    C, _ = Ferrite.create_constraint_matrix(ch)
    Kc = sparse(C' * K * C); Kgc = sparse(C' * (-Kg) * C); F = lu(Kc)
    op = LinearMap{Float64}(x -> F \ (Kgc * x), size(Kc, 1); ismutating = false)
    decomp, _ = partialschur(op; nev = 3, tol = 1e-9, which = ArnoldiMethod.LR()); μ, Ψ = partialeigen(decomp)
    k = argmax(real.(μ)); φ = real.(Ψ[:, k]); P = P_ref / real(μ[k])
    Uel = φ' * Kc * φ; Ug = φ' * Kgc * φ                      # P_cr = P_ref * Uel / Ug
    a = C * φ
    θ = [a[nd[id(i_webmid, j)][6]] for j in 1:nz]            # twist per station (end stations: u = v = 0 → treat θ = 0)
    θ[1] = 0.0; θ[end] = 0.0
    h = Z[2] - Z[1]
    θ1 = diff(θ) ./ h                                        # θ' at mid-segments
    θ2 = [(θ[j+1] - 2θ[j] + θ[j-1]) / h^2 for j in 2:nz-1]  # θ'' at interior stations
    Iθ1 = sum(θ1 .^ 2) * h; Iθ2 = sum(θ2 .^ 2) * h
    x_o = xs - sp.xc; r_o2 = (sp.Ixx + sp.Iyy) / sp.A + x_o^2
    Uel_th = G * J_SV * Iθ1 + E * sp.Cw * Iθ2                # 2 × strain energy for the same θ(z)
    Ug_th = P_ref * r_o2 * Iθ1                               # 2 × geometric term with the reference load
    # fit θ to sin: fraction of the sine component
    s = sin.(π .* Z ./ L); frac = (θ' * s) / (s' * s); resid = norm(θ .- frac .* s) / norm(θ)
    @printf("L = %.0f in: P_t shell = %.2f kips; mode θ(z) sine content residual %.3f\n", L, P / 1e3, resid)
    @printf("  elastic:   φᵀKφ = %.6e   theory GJ∫θ'² + ECw∫θ''² = %.6e  (GJ part %.2e, ECw part %.2e)  ratio shell/theory = %.3f\n", Uel, Uel_th, G * J_SV * Iθ1, E * sp.Cw * Iθ2, Uel / Uel_th)
    @printf("  geometric: φᵀ(−Kg)φ = %.6e   theory P r_o² ∫θ'² = %.6e   ratio shell/theory = %.3f   (r_o² = %.3f in²)\n", Ug, Ug_th, Ug / Ug_th, r_o2)
    @printf("  → P_cr from theory pieces with this θ(z): %.2f kips; classical P_t = %.2f kips\n", P_ref * Uel_th / Ug_th / 1e3, (G * J_SV + π^2 * E * sp.Cw / L^2) / (sp.A * r_o2) / 1e3)
end
torsion_mode_energy(; L = 44.0)
torsion_mode_energy(; L = 400.0)
