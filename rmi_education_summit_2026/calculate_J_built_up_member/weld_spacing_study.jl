# weld_spacing_study.jl
#
# Ferrite.jl study of the effective torsion constant J_eff of the two-C welded OneRack column as a
# function of weld spacing, compared with the proposed prediction equation
#
#     J_combined = J_tube · (w_length / w_spacing) + J_back + J_front · (1 − w_length / w_spacing)
#
# with J_tube the Bredt closed-cell constant of the C1 web + C1 flanges + C2 web (4 A_m² t / s_m),
# J_front = J of the lipped C whose lips are welded to the other web (C1), J_back = J of the C whose web
# closes the cell (C2); both are the single-C Saint-Venant value.
#
# Welds are w_length long, equally spaced from end to end (first and last weld flush with the member
# ends, as in combo_1_J_r5.jl): n welds, w_spacing = (L − w_length) / (n − 1). Series:
#   A. w_length = 3 in, L = 111 in, w_spacing = 108 … 3 in (2 … 37 welds; 3 in = fully welded)
#   B. w_length = 6 in, L = 111 in, w_spacing = 105 … 6 in
#   C. w_length = 1.5 in, L = 111 in, w_spacing = 109.5 … 6 in
#   D. w_length = 3 in, w_spacing = 18 in, L = 57, 111, 219 in (member-length dependence)
#
# Run:  julia --project=. weld_spacing_study.jl        (writes weld_spacing_results.csv)

include(joinpath(@__DIR__, "calculate_J_built_up_member.jl"))

const J_C = 1.334851e-3                                        # exact Saint-Venant J of one C (in⁴)
const a_cell = shape.B - shape.t                               # closed cell centerline side, 2.926 in
const J_tube = 4 * a_cell^4 * shape.t / (4 * a_cell)           # Bredt: 4 A_m² t / s_m, square cell
J_equation(wl, ws) = J_tube * wl / ws + J_C + J_C * (1 - wl / ws)

function spacing_locations(L, wl, n)
    n == 1 && return [L / 2]
    return collect(range(wl / 2, L - wl / 2, n))
end

function run_case(L, wl, n; kw...)
    locs = spacing_locations(L, wl, n)
    ws = n > 1 ? locs[2] - locs[1] : Inf
    r = built_up_J(; L, weld_length = wl, weld_locations = locs, verbose = false, kw...)
    return (; L, wl, n, ws, ratio = wl / ws, J_eff = r.J_eff, J_eq = J_equation(wl, ws), Fnet = hypot(r.Fx, r.Fy) / abs(r.T), ndofs = r.ndofs)
end

if abspath(PROGRAM_FILE) == @__FILE__
    @printf("J_tube (Bredt, %.3f in square cell) = %.4f in⁴,  J_front = J_back = J_C = %.4e in⁴\n\n", a_cell, J_tube, J_C)
    series = [
        ("A", 111.0, 3.0, [2, 3, 4, 5, 7, 10, 13, 19, 37]),          # spacing 108, 54, 36, 27, 18, 12, 9, 6, 3
        ("B", 111.0, 6.0, [2, 3, 4, 6, 8, 11, 16]),                   # spacing 105, 52.5, 35, 21, 15, 10.5, 7
        ("C", 111.0, 1.5, [2, 3, 4, 5, 7, 10, 13, 19, 37]),           # spacing 109.5, 54.75, 36.5, 27.4, 18.25, 12.2, 9.1, 6.1, 3.04
        ("D", 57.0, 3.0, [4]), ("D", 219.0, 3.0, [13]),               # spacing 18 at L = 57 and 219 (111 is in series A)
    ]
    rows = []
    for (name, L, wl, ns) in series, n in ns
        t_run = @elapsed c = run_case(L, wl, n)
        @printf("series %s: L = %5.1f, w_length = %.1f, %2d welds, w_spacing = %7.2f (w/s = %.3f): J_eff = %.4e in⁴, equation = %.4e, Ferrite/eq = %6.3f, |ΣF|/T = %.0e (%.1f s)\n",
                name, L, wl, n, c.ws, c.ratio, c.J_eff, c.J_eq, c.J_eff / c.J_eq, c.Fnet, t_run)
        push!(rows, (series = name, c...))
    end
    open(joinpath(@__DIR__, "weld_spacing_results.csv"), "w") do io
        println(io, "series,L_in,w_length_in,n_welds,w_spacing_in,w_length_over_w_spacing,J_eff_ferrite_in4,J_equation_in4,ferrite_over_equation,sumF_over_T,ndofs,J_tube_in4,J_C_in4")
        for r in rows
            println(io, join([r.series, r.L, r.wl, r.n, r.ws, r.ratio, r.J_eff, r.J_eq, r.J_eff / r.J_eq, r.Fnet, r.ndofs, J_tube, J_C], ","))
        end
    end
    println("\nWrote weld_spacing_results.csv")
end
