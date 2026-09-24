# ft_analytical_check.jl — compare the shell-model global FT buckling load with Timoshenko's flexural-torsional
# formula for a section symmetric about the horizontal (x) axis:
#     P_FT = [ (P_ey + P_t) − sqrt((P_ey + P_t)² − 4 β P_ey P_t) ] / (2 β),   β = 1 − (x_o / r_o)²
# with P_ey = π² E I_x / L² (bending about x, deflection v), P_t = (G J + π² E C_w / L²) / (A r_o²),
# r_o² = (I_x + I_y)/A + x_o².  P_ey, P_t are also obtained from the shell model by restraining the rigid-section dofs,
# and x_o is backed out of the FT mode shape ratio  (P_ey − P) v = −P x_o φ.
include(joinpath(@__DIR__, "buckling_built_up.jl"))
L = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 44.0
tag = L == 44.0 ? "" : "_L$(Int(round(L)))"
X1, Y1 = centerline(shape); sp = section_properties(X1, Y1, shape.t)
A2 = 2sp.A; Ix2 = 2sp.Ixx; Iy2 = 2 * (sp.Iyy + sp.A * (shape.B / 2)^2)
P_ey = π^2 * E * Ix2 / L^2; P_ex = π^2 * E * Iy2 / L^2
J_eff = 0.7678                                            # in⁴, r5 weld layout (3 in welds at 18 in)

println("shell: global FT (u, v, θ free)"); g  = built_up_buckling(; L, mode = :global, nev = 2, verbose = false)
println("shell: pure y-flexure (u = θ = 0)"); gv = built_up_buckling(; L, mode = :global, nev = 2, restrain = [:u, :θ], verbose = false)
println("shell: pure x-flexure (v = θ = 0)"); gu = built_up_buckling(; L, mode = :global, nev = 2, restrain = [:v, :θ], verbose = false)
P_FT = g.P_cr[1]; P_ex_sh = gu.P_cr[1]

# shear center by statics: transverse force F_y at mid-length applied at x_p on the rigid-section model (pinned ends);
# the section twist at mid-length is linear in x_p and vanishes when x_p is the shear center
function twist_for_load_at(g, x_p)
    m = g.id(1, g.i_webmid, g.j_mid); dm = g.nd[m]; xm = g.grid.nodes[m].x[1]
    f = zeros(size(g.K0, 1)); f[dm[2]] = 1000.0; f[dm[6]] = 1000.0 * (x_p - xm)
    K = copy(g.K0); apply!(K, f, g.ch); u = K \ f; apply!(u, g.ch)
    return u[dm[6]]
end
θa = twist_for_load_at(g, 0.0); θb = twist_for_load_at(g, 3.0)
x_s_static = -θa / ((θb - θa) / 3.0)
@printf("Shear center by statics (zero twist under a mid-length transverse force): x_s = %.3f in\n", x_s_static)

# shear-center offset from the FT mode: v/φ at mid-length (section rigid rotation θ and translation V)
P_ey_sh = gv.P_cr[1]
c = g.info[1]
xs = vcat([g.Xs[s][i] for s in 1:2 for i in 1:g.nn]); ys = vcat([g.Ys[s][i] for s in 1:2 for i in 1:g.nn])
a = g.modes[1]; jm = argmin(abs.(g.Z .- L / 2))
ux = [a[g.nd[g.id(s, i, jm)][1]] for s in 1:2 for i in 1:g.nn]; uy = [a[g.nd[g.id(s, i, jm)][2]] for s in 1:2 for i in 1:g.nn]
V = mean(uy); φ = ls_rotation(xs, ys, ux, uy)
x_o_mode = -(P_ey_sh - P_FT) / P_FT * V / φ               # Timoshenko: (P_ey − P) v = −P x_o φ
xc2 = mean(xs); yc2 = mean(ys)                            # composite centroid (equal areas at every node is close enough for xc: use area-weighted below)
xc2 = (sp.xc + (sp.xc + shape.B)) / 2
x_o = x_s_static - xc2                                    # use the static shear center
@printf("x_o from statics = %.3f in;  x_o implied by the FT mode ratio = %.3f in\n", x_o, x_o_mode)
r_o2 = (Ix2 + Iy2) / A2 + x_o^2; β = 1 - x_o^2 / r_o2
yc2 = sp.yc
println("shell: pure torsion about the shear center (x = $(round(xc2 + x_o, digits = 3)), y = $(round(yc2, digits = 3)))")
gt = built_up_buckling(; L, mode = :global, nev = 2, restrain = [:u, :v], twist_center = (xc2 + x_o, yc2), verbose = false)
P_t_sh = gt.P_cr[1]
timo(Pey, Pt) = ((Pey + Pt) - sqrt((Pey + Pt)^2 - 4β * Pey * Pt)) / (2β)

P_t_J = G * J_eff / r_o2                                  # torsional buckling load from G J_eff alone (C_w = 0): P_t = (G J + π² E C_w/L²)/r_o²
Cw_eff = (P_t_sh * r_o2 - G * J_eff) * L^2 / (π^2 * E)        # warping constant implied by the shell P_t (P_t r_o² = G J + π² E C_w/L²)

@printf("\nComposite section: A = %.4f in², Ix = %.4f, Iy = %.4f in⁴, centroid x = %.3f in; L = %.1f in\n", A2, Ix2, Iy2, xc2, L)
@printf("Shell eigen loads (kips): P_FT = %.1f,  pure torsion P_t = %.1f,  pure y-flexure P_ey = %.1f (Euler %.1f),  pure x-flexure P_ex = %.1f (Euler %.1f)\n",
        P_FT/1e3, P_t_sh/1e3, P_ey_sh/1e3, P_ey/1e3, P_ex_sh/1e3, P_ex/1e3)
@printf("FT mode at mid-length: v/φ = %.3f in;  using x_o = %.3f in (shear center at x = %.3f in), r_o = %.3f in, β = 1 − (x_o/r_o)² = %.3f\n", V/φ, x_o, xc2 + x_o, sqrt(r_o2), β)
@printf("\nTimoshenko P_FT with shell P_ey and shell P_t:            %.1f kips   (shell FT %.1f, ratio %.3f)\n", timo(P_ey_sh, P_t_sh)/1e3, P_FT/1e3, timo(P_ey_sh, P_t_sh)/P_FT)
@printf("Timoshenko P_FT with Euler P_ey and P_t = G J_eff/r_o² (C_w = 0):   %.1f kips   (P_t = %.1f kips)\n", timo(P_ey, P_t_J)/1e3, P_t_J/1e3)
@printf("Shell pure-torsion load ↔ torsional stiffness P_t r_o² = %.0f kip-in²  (G J_eff = %.0f kip-in²)\n", P_t_sh * r_o2 / 1e3, G * J_eff / 1e3)
@printf("Warping constant implied by the shell pure-torsion load:  C_w,eff = %.3f in⁶  (G J_eff = %.0f, π²E C_w/L² = %.0f kip-in²)\n", Cw_eff, G*J_eff/1e3, π^2*E*Cw_eff/L^2/1e3)
open(joinpath(@__DIR__, "ft_analytical_check$(tag).csv"), "w") do io
    println(io, "quantity,value,unit")
    for (k, v, u) in (("P_FT_shell", P_FT/1e3, "kips"), ("P_t_shell", P_t_sh/1e3, "kips"), ("P_ey_shell", P_ey_sh/1e3, "kips"), ("P_ey_Euler", P_ey/1e3, "kips"),
                      ("P_ex_shell", P_ex_sh/1e3, "kips"), ("P_ex_Euler", P_ex/1e3, "kips"), ("x_o", x_o, "in"), ("r_o", sqrt(r_o2), "in"), ("beta", β, ""),
                      ("P_FT_Timoshenko_shell_inputs", timo(P_ey_sh, P_t_sh)/1e3, "kips"), ("P_t_from_GJeff_Cw0", P_t_J/1e3, "kips"), ("torsional_stiffness_shell_Pt_ro2", P_t_sh * r_o2 / 1e3, "kip-in^2"),
                      ("P_FT_Timoshenko_GJeff_Cw0", timo(P_ey, P_t_J)/1e3, "kips"), ("Cw_eff_implied", Cw_eff, "in^6"), ("J_eff", J_eff, "in^4"))
        println(io, "$k,$v,$u")
    end
end
