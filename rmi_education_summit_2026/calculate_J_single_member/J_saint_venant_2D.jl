# J_saint_venant_2D.jl
#
# Reference Saint-Venant torsion constant of the upright cross section (with its rounded corners and
# real thickness) from the Prandtl stress function on the 2D section with Ferrite.jl:
#
#     ∇²φ = −2 in A,   φ = 0 on the boundary,   J = 2 ∫_A φ dA
#
# This is the exact elasticity value the shell model should approach. The section mesh is built from
# the left/right surface node lines that CrossSectionGeometry generates around the centerline.
#
# Run:  julia --project=. J_saint_venant_2D.jl

using Ferrite, Tensors, LinearAlgebra, Printf
using CrossSectionGeometry

const shape = (t = 0.074, lip = 0.75, D = 3.0, B = 3.0, R = 0.125 + 0.074, length = 108.0 + 3.0)

# quadrilateral 2D mesh of a thin-walled strip between two surface polylines, n_t layers through thickness
function strip_grid(left, right; n_t = 4)
    nn = length(left)
    nodes = Node{2, Float64}[]
    for k in 0:n_t, i in 1:nn
        p = left[i] .+ (right[i] .- left[i]) .* (k / n_t)
        push!(nodes, Node(Vec((p[1], p[2]))))
    end
    id(i, k) = k * nn + i
    # orientation check on the first cell
    p1 = nodes[id(1, 0)].x; p2 = nodes[id(2, 0)].x; p4 = nodes[id(1, 1)].x
    ccw = (p2[1] - p1[1]) * (p4[2] - p1[2]) - (p2[2] - p1[2]) * (p4[1] - p1[1]) > 0
    cells = Quadrilateral[]
    for k in 0:n_t-1, i in 1:nn-1
        c = ccw ? (id(i, k), id(i + 1, k), id(i + 1, k + 1), id(i, k + 1)) :
                  (id(i, k), id(i, k + 1), id(i + 1, k + 1), id(i + 1, k))
        push!(cells, Quadrilateral(c))
    end
    boundary = Set{Int}()
    for k in 0:n_t, i in 1:nn
        (k == 0 || k == n_t || i == 1 || i == nn) && push!(boundary, id(i, k))
    end
    return Grid(cells, nodes), boundary
end

function saint_venant_J(grid, boundary; order = 2)
    ip = Lagrange{RefQuadrilateral, order}()
    qr = QuadratureRule{RefQuadrilateral}(order + 1)
    cv = CellValues(qr, ip, Lagrange{RefQuadrilateral, 1}())
    dh = DofHandler(grid); add!(dh, :φ, ip); close!(dh)
    # φ = 0 on the boundary: all dofs on boundary facets (works for any interpolation order)
    facets = Set{FacetIndex}()
    for (c, cell) in enumerate(grid.cells), (fi, facet) in enumerate(Ferrite.facets(cell))
        all(n -> n in boundary, facet) && push!(facets, FacetIndex(c, fi))
    end
    ch = ConstraintHandler(dh); add!(ch, Dirichlet(:φ, facets, (x, t) -> 0.0)); close!(ch)
    K = allocate_matrix(dh); f = zeros(ndofs(dh))
    assembler = start_assemble(K, f)
    n = getnbasefunctions(cv); ke = zeros(n, n); fe = zeros(n)
    A = 0.0
    for cell in CellIterator(dh)
        reinit!(cv, cell); fill!(ke, 0.0); fill!(fe, 0.0)
        for q in 1:getnquadpoints(cv)
            dΩ = getdetJdV(cv, q); A += dΩ
            for i in 1:n
                δφ = shape_value(cv, q, i); ∇δφ = shape_gradient(cv, q, i)
                fe[i] += 2 * δφ * dΩ
                for j in 1:n
                    ke[i, j] += ∇δφ ⋅ shape_gradient(cv, q, j) * dΩ
                end
            end
        end
        assemble!(assembler, celldofs(cell), ke, fe)
    end
    apply!(K, f, ch)
    φ = K \ f
    apply!(φ, ch)
    J = 0.0
    for cell in CellIterator(dh)
        reinit!(cv, cell)
        φe = φ[celldofs(cell)]
        for q in 1:getnquadpoints(cv)
            J += 2 * function_value(cv, q, φe) * getdetJdV(cv, q)
        end
    end
    return J, A, ndofs(dh)
end

function upright_section(shape; n_flat, n_corner)
    L = [shape.lip, shape.B, shape.D, shape.B, shape.lip]
    θ = [π / 2, π, -π / 2, 0.0, π / 2]
    cs = CrossSectionGeometry.create_thin_walled_cross_section_geometry(L, θ, fill(n_flat, 5), fill(shape.R, 4), fill(n_corner, 4), shape.t;
                                                                         centerline = "to left", offset = (0.0, 0.0))
    C = cs.centerline_node_XY
    b = sum(hypot(C[i+1][1] - C[i][1], C[i+1][2] - C[i][2]) for i in 1:length(C)-1)
    return cs.left_surface_node_XY, cs.right_surface_node_XY, b
end

if abspath(PROGRAM_FILE) == @__FILE__
    t = shape.t
    println("Validation: rectangle 3 in × $(t) in")
    b = 3.0
    J_exact = b * t^3 / 3 * (1 - 0.63 * t / b + 0.052 * (t / b)^5)
    for (nb, nt, order) in ((120, 8, 1), (240, 12, 1), (60, 2, 2), (120, 4, 2), (240, 8, 2))
        left = [[b * i / nb, -t / 2] for i in 0:nb]; right = [[b * i / nb, t / 2] for i in 0:nb]
        g, bd = strip_grid(left, right; n_t = nt)
        J, A, nd = saint_venant_J(g, bd; order)
        @printf("  mesh %3d × %2d, order %d (%6d dofs): J = %.6e   J/(bt³/3) = %.5f   J/J_exact = %.5f\n", nb, nt, order, nd, J, J / (b * t^3 / 3), J / J_exact)
    end
    println("\nUpright lipped C with rounded corners (R = $(shape.R) to the centerline)")
    for (nf, nc, nt, order) in ((80, 24, 10, 1), (40, 12, 4, 2), (80, 24, 6, 2), (160, 48, 8, 2), (240, 72, 12, 2))
        left, right, bcl = upright_section(shape; n_flat = nf, n_corner = nc)
        g, bd = strip_grid(left, right; n_t = nt)
        J, A, nd = saint_venant_J(g, bd; order)
        Jt = bcl * t^3 / 3
        @printf("  mesh %3d/flat %2d/corner %2d thru-t, order %d (%7d dofs): A = %.5f in², J = %.6e in⁴, Σbt³/3 = %.6e, J/(Σbt³/3) = %.5f\n", nf, nc, nt, order, nd, A, J, Jt, J / Jt)
    end
end
