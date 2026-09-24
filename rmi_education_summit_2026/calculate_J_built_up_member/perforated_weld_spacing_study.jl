# perforated_weld_spacing_study.jl — J_eff of the welded pair vs weld spacing (3 in welds, L = 111 in) on the Gmsh mesh,
# with and without the perforation pattern.   Run:  julia --project=. perforated_weld_spacing_study.jl
include(joinpath(@__DIR__, "perforated_J_gmsh.jl"))
L = 111.0; wl = 3.0
open(joinpath(@__DIR__, "perforated_weld_spacing_results.csv"), "w") do io
    println(io, "n_welds,w_spacing_in,J_eff_noholes_in4,J_eff_perforated_in4,ratio")
    for n in (2, 3, 4, 5, 7, 10, 13, 19, 37)
        locs = collect(range(wl / 2, L - wl / 2, n)); ws = n > 1 ? locs[2] - locs[1] : Inf
        @printf("n = %2d welds, spacing %.2f in\n", n, ws)
        r0 = twist_J_gmsh(; nshapes = 2, perforated = false, L, weld_length = wl, weld_locations = locs)
        r1 = twist_J_gmsh(; nshapes = 2, perforated = true, L, weld_length = wl, weld_locations = locs)
        println(io, "$n,$ws,$(r0.J),$(r1.J),$(r1.J / r0.J)"); flush(io)
    end
end
println("Wrote perforated_weld_spacing_results.csv")
