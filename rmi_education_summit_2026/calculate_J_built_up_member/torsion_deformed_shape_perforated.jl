# torsion_deformed_shape_perforated.jl — static twist of the perforated welded pair (r5 layout) on the Gmsh mesh;
# writes the mesh and displacements for make_torsion_html_mesh.jl.
include(joinpath(@__DIR__, "perforated_J_gmsh.jl"))
r = twist_J_gmsh(; nshapes = 2, perforated = true)
g = r.grid; nd = r.nd; u = r.u
open(joinpath(@__DIR__, "torsion_perforated_nodes.csv"), "w") do io
    println(io, "x,y,z,ux,uy,uz")
    for (n, node) in enumerate(g.nodes)
        p = node.x
        if haskey(nd, n); d = nd[n]; println(io, join([p[1], p[2], p[3], u[d[1]], u[d[2]], u[d[3]]], ","))
        else; println(io, join([p[1], p[2], p[3], 0.0, 0.0, 0.0], ",")); end
    end
end
open(joinpath(@__DIR__, "torsion_perforated_cells.csv"), "w") do io
    for c in g.cells; println(io, join(c.nodes, ",")); end
end
open(joinpath(@__DIR__, "torsion_perforated_summary.csv"), "w") do io
    println(io, "L_in,T_lbf_in,J_eff_in4,betao_rad,Xc,Yc,quads,tris")
    println(io, join([shape.length, r.T, r.J, βo, shape.B / 2, shape.D / 2, r.nquad, r.ntri], ","))
end
println("wrote torsion_perforated_{nodes,cells,summary}.csv")
