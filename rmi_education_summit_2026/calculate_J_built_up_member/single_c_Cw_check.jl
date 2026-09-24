# single_c_Cw_check.jl — (1) independent sectorial-coordinate C_w and shear center of the C centerline polyline;
# (2) static warping-fixed twist of the shell model (u_z = 0 on both end sections, rigid twist φ0 at z = L,
#     rigid in-plane sections) → effective E C_w from  T = 12 E C_w φ0/L³ + G J φ0/L  (Vlasov, small G J).
include(joinpath(@__DIR__, "buckling_built_up.jl"))

# ---- (1) sectorial properties of an open thin-walled polyline
function sectorial(X, Y, t)
    n = length(X); ds = [hypot(X[i+1]-X[i], Y[i+1]-Y[i]) for i in 1:n-1]
    A = sum(ds) * t
    xc = sum((X[i]+X[i+1])/2 * ds[i] for i in 1:n-1) * t / A; yc = sum((Y[i]+Y[i+1])/2 * ds[i] for i in 1:n-1) * t / A
    # second moments (linear segments, exact)
    Ix = 0.0; Iy = 0.0; Ixy = 0.0
    for i in 1:n-1
        x1, x2 = X[i]-xc, X[i+1]-xc; y1, y2 = Y[i]-yc, Y[i+1]-yc
        Ix += t*ds[i]*(y1^2 + y1*y2 + y2^2)/3; Iy += t*ds[i]*(x1^2 + x1*x2 + x2^2)/3; Ixy += t*ds[i]*(2x1*y1 + x1*y2 + x2*y1 + 2x2*y2)/6
    end
    # sectorial area about the centroid, ω(s) = ∫ (x dy − y dx)
    ω = zeros(n)
    for i in 1:n-1
        x1, x2 = X[i]-xc, X[i+1]-xc; y1, y2 = Y[i]-yc, Y[i+1]-yc
        ω[i+1] = ω[i] + (x1*y2 - x2*y1)
    end
    # sectorial products ∫ω x dA, ∫ω y dA, ∫ω dA (linear ω along segments)
    Iωx = 0.0; Iωy = 0.0; Sω = 0.0
    for i in 1:n-1
        x1, x2 = X[i]-xc, X[i+1]-xc; y1, y2 = Y[i]-yc, Y[i+1]-yc; w1, w2 = ω[i], ω[i+1]
        Iωx += t*ds[i]*(2w1*x1 + w1*x2 + w2*x1 + 2w2*x2)/6; Iωy += t*ds[i]*(2w1*y1 + w1*y2 + w2*y1 + 2w2*y2)/6; Sω += t*ds[i]*(w1+w2)/2
    end
    # shear center (pole) relative to the centroid
    D = Ix*Iy - Ixy^2
    xs = (Iy*Iωy - Ixy*Iωx) / D;  ys = -(Ix*Iωx - Ixy*Iωy) / D
    # normalized warping about the shear center: ω_s = ω − ω̄ + ys x − xs y  (then remove the mean)
    ωs = [ω[i] + ys*(X[i]-xc) - xs*(Y[i]-yc) for i in 1:n]
    Sωs = sum(t*ds[i]*(ωs[i]+ωs[i+1])/2 for i in 1:n-1); ωs .-= Sωs / A
    Cw = sum(t*ds[i]*(ωs[i]^2 + ωs[i]*ωs[i+1] + ωs[i+1]^2)/3 for i in 1:n-1)
    return (; A, xc, yc, Ix, Iy, Ixy, xs = xc + xs, ys = yc + ys, Cw)
end
X, Y = centerline(shape); sp = section_properties(X, Y, shape.t); sec = sectorial(X, Y, shape.t)
@printf("cutwp:     A = %.4f, Ix = %.4f, Iy = %.4f, shear center (%.4f, %.4f), Cw = %.4f in⁶\n", sp.A, sp.Ixx, sp.Iyy, sp.xs, sp.ys, sp.Cw)
@printf("sectorial: A = %.4f, Ix = %.4f, Iy = %.4f, shear center (%.4f, %.4f), Cw = %.4f in⁶\n", sec.A, sec.Ix, sec.Iy, sec.xs, sec.ys, sec.Cw)

# ---- (2) warping-fixed static twist on the shell, rigid sections, rotation about the shear center
function warping_fixed_twist(; L, rigid = true, n_flat = 4, n_corner = 5)
    X, Y = centerline(shape; n_flat, n_corner); Z = collect(range(0.0, L, Int(round(2L)) + 1))
    nn = length(X); nz = length(Z); grid, id = extruded_grid(X, Y, Z); dh, nd = shell_dofhandler(grid)
    xs, ys = sp.xs, sp.ys; i_webmid = argmin(X .^ 2 .+ (Y .- shape.D / 2) .^ 2)
    ch = ConstraintHandler(dh)
    e0 = Set(id(i, 1) for i in 1:nn); eL = Set(id(i, nz) for i in 1:nn)
    add!(ch, Dirichlet(:u, e0, (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3]))                                     # fixed, warping fixed
    add!(ch, Dirichlet(:u, eL, (x, t) -> [-βo * (x[2] - ys), βo * (x[1] - xs), 0.0], [1, 2, 3]))             # twist βo, warping fixed
    if rigid
        for j in 2:nz-1
            m = id(i_webmid, j); dm = nd[m]
            for i in 1:nn
                n = id(i, j); n == m && continue; p = grid.nodes[n].x; dn = nd[n]; pm = grid.nodes[m].x
                add!(ch, AffineConstraint(dn[1], [dm[1] => 1.0, dm[6] => -(p[2] - pm[2])], 0.0))
                add!(ch, AffineConstraint(dn[2], [dm[2] => 1.0, dm[6] => (p[1] - pm[1])], 0.0))
            end
        end
    end
    close!(ch); update!(ch, 0.0)
    K = allocate_matrix(dh, ch)
    K = Q.assemble_global_Ke!(K, dh, QuadratureRule{RefQuadrilateral}(2), QuadratureRule{RefQuadrilateral}(3), IP4(), IP6(), E, ν, shape.t)
    K0 = copy(K); f = zeros(ndofs(dh)); apply!(K, f, ch); u = K \ f; apply!(u, ch); R = K0 * u
    T = sum((grid.nodes[n].x[1] - xs) * R[nd[n][2]] - (grid.nodes[n].x[2] - ys) * R[nd[n][1]] for n in eL)
    return T
end
J_SV = 1.334851e-3
for L in (44.0, 88.0), rigid in (false, true)
    T = warping_fixed_twist(; L, rigid)
    # Vlasov fixed-fixed warping: T = G J φ0 k / (k L − 2 tanh(k L / 2)),  k² = G J / (E C_w);  invert for E C_w by bisection
    f(ECw) = (k = sqrt(G * J_SV / ECw); G * J_SV * βo * k / (k * L - 2 * tanh(k * L / 2)))
    lo, hi = 1e3, 1e7
    for _ in 1:200; mid = sqrt(lo * hi); (f(mid) > T ? (hi = mid) : (lo = mid)); end
    ECw = sqrt(lo * hi)
    @printf("warping-fixed twist, L = %3.0f in, %s sections: T/φ0 = %.3f kip-in/rad → E C_w = %.0f kip-in⁴ → C_w,shell = %.4f in⁶ (cutwp %.4f, ratio %.3f);  small-GJ estimate C_w = %.4f\n",
            L, rigid ? "rigid" : "free ", T / βo / 1e3, ECw / 1e3, ECw / E, sp.Cw, ECw / E / sp.Cw, (T / βo - G * J_SV / L) * L^3 / (12E))
end
