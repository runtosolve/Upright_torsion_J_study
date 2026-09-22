# torsion_deformed_shape.jl — run the r5-layout static twist (3 in welds at 18 in, L = 111 in) and write the
# node coordinates and displacements for the interactive page (make_torsion_html.jl).
include(joinpath(@__DIR__, "calculate_J_built_up_member.jl"))
r = built_up_J(; verbose = true)
open(joinpath(@__DIR__, "torsion_deformed_r5.csv"), "w") do io
    println(io, "shape,i,j,x,y,z,ux,uy,uz")
    for j in 1:r.nz, s in 1:2, i in 1:r.nn
        n = r.id(s, i, j); p = r.grid.nodes[n].x; d = r.nd[n]
        println(io, join([s, i, j, p[1], p[2], p[3], r.u[d[1]], r.u[d[2]], r.u[d[3]]], ","))
    end
end
open(joinpath(@__DIR__, "torsion_summary_r5.csv"), "w") do io
    println(io, "L_in,weld_length_in,n_welds,T_lbf_in,J_eff_in4,betao_rad,Xc,Yc")
    println(io, join([r.L, r.weld_length, r.n_welds, r.T, r.J_eff, βo, r.Xc, r.Yc], ","))
end
println("wrote torsion_deformed_r5.csv, torsion_summary_r5.csv")
