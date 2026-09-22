# J of a single OneRack upright with Ferrite.jl — Moen (2008) §4.2.7.3.2.3 static twist method

Goal: verify that the uniform torsion differential equation, T = G J dβ/dz (Eq. 4.12 with warping-free
ends so d³β/dz³ = 0), can be used with a Ferrite.jl + QuadShellFiniteElement.jl shell model to reliably
calculate J for one upright of the built-up column (SHAPE_1 in `combo_1/J/combo_1_J_r5.jl`: 3 × 3 × 0.074 in
lipped C, 0.75 in lips, 0.199 in centerline corner radius, L = 111 in, E = 29500 ksi, ν = 0.3).

## Files

| file | purpose |
|---|---|
| `calculate_J_single_member.jl` | main study: shell model, Fig. 4.41 boundary conditions, twist profile, J by Eq. 4.13 and by the slope of β(z), mesh / length / rotation-center / drilling sensitivity. Writes `J_results.csv`, `twist_profile_abaqus_mesh.csv`, `run_log.txt` |
| `J_saint_venant_2D.jl` | exact reference: Prandtl stress function ∇²φ = −2 on a 2D mesh of the real section (rounded corners, true thickness), J = 2∫φ dA, quadratic elements. Also validates against the exact rectangle formula |
| `fold_diagnostics.jl` | why the drilling treatment matters: flat strip, folded strips, slit tubes, sharp and rounded C. Writes `fold_diagnostics_log.txt` |
| `make_figure.jl` | `J_single_member_summary.png` (twist profile, Fig. 4.42 analog, and J summary) |

Run from this folder with `julia --project=.` (Project.toml pins Ferrite 1.7, QuadShellFiniteElement,
CrossSectionGeometry, SectionProperties from `~/.julia/dev`, CairoMakie). The drilling treatment is
selected with the `drilling` keyword of `QuadShellFiniteElement.assemble_global_Ke!`
(`:hughes_brezzi`, now the package default, or `:penalty`), added to the package on 2026-09-21 as a
result of this study.

## Model (Figure 4.41 of the dissertation)

- Mesh: the same centerline discretization and 0.5 in element length as the Abaqus model (4 elements per
  flat, 5 per corner, 222 along; 41 nodes per section, 54,858 dofs).
- z = 0: all nodes fixed in X and Y (twist and translation restrained, warping free); web mid-height node
  fixed in Z.
- z = L: all nodes given the rigid rotation βo about the shear center (u = −βo (Y − Ys), v = βo (X − Xs)),
  warping free — the Dirichlet equivalent of the Abaqus kinematic coupling to the shear-center node.
- T = resultant moment about the shear center of the X, Y reactions at z = L. β(z) is measured at every
  station as the relative rotation of the web/flange corners and as a least-squares rigid rotation.

## Results

Exact Saint-Venant J of the section (2D, converged): **J_SV = 1.33485e-3 in⁴** = 0.9956 × Σbt³/3
(Σbt³/3 = 1.34015e-3 in⁴; the 0.44 % difference is the free-edge effect at the two lip tips).

| case | J (in⁴) | J / J_SV |
|---|---|---|
| Abaqus mesh, Hughes–Brezzi drilling | 1.32558e-3 | 0.993 |
| same, Cs = 0 | 1.32909e-3 | 0.996 |
| same, twist applied about the centroid | 1.32558e-3 | 0.993 |
| same, L = 55.5 in | 1.32462e-3 | 0.992 |
| coarse 2/flat 2/corner 111 along | 1.31912e-3 | 0.988 |
| fine 8/flat 10/corner 444 along | 1.32735e-3 | 0.994 |
| Abaqus mesh, original absolute drilling penalty 1/100 | 2.86613e-3 | 2.147 |
| fine mesh, original absolute drilling penalty 1/100 | 4.57967e-3 | 3.431 |
| Abaqus mesh, drilling penalty removed (1e-6, diagnostic) | 1.18125e-3 | 0.885 |

In every case: β(z) is linear (R² = 1.00000000, max deviation from the fit 2e-5 βo), the torque at z = 0
balances the torque at z = L to 1e-6, and J from Eq. 4.13 (T L / G βo) equals J from the fitted slope
(T / G dβ/dz) to 4 digits. The method itself is therefore verified. The result is independent of the
member length and of the point about which the end twist is applied.

## Finding: the element's drilling treatment controls J for a folded / curved section

The original package treatment (absolute penalty k = min rotational diagonal / 100 on each drilling rotation) gives
J 2.1 × too large on the Abaqus mesh and 3.4 × too large on the finer mesh, with a spurious net in-plane
end reaction. Cause: in Saint-Venant torsion every plate strip rotates in-plane at the rate θn = β′h
(h = distance from the shear center to the strip), and the penalty on θn itself resists this true motion.
The error grows with the number of fold lines (corner refinement) because each corner node's drilling
direction picks up the neighbor's bending rotation.

Removing the penalty (1e-6) makes J 12 % too small. Cause: without drilling stiffness a strip's twisting
moment cannot be handed across a fold line into the next strip, so every fold acts as a Mindlin free edge
and costs ≈ 0.63 t⁴/3 (`fold_diagnostics.jl`: a 3 in strip folded in two at any angle loses exactly one
strip's free-edge correction; a 60-segment slit tube drops to 0.60 Σbt³/3).

A Hughes–Brezzi drilling term, γ∫(θn − ½(∂v/∂x − ∂u/∂y))² dA with γ = G t (insensitive to γ between
0.1 G t and 10 G t), is energy-free for the true in-plane rotation and transfers the twisting moment across
folds. With it the rounded C gives 0.987 Σbt³/3 on every mesh tested, the flat strip result is unchanged,
and J_end / J_SV = 0.993 on the Abaqus mesh (the remaining 0.7 % is the bilinear Mindlin element's
discretization error, which shrinks with refinement).

Outcome: the Hughes–Brezzi drilling term was added to QuadShellFiniteElement.jl and
TriShellFiniteElement.jl on 2026-09-21 as the default (`drilling = :hughes_brezzi`; `:penalty` keeps the
original form, which the MATLAB/Octave reference tests use). Both packages now include a
"torsion of folded strips" test. Any earlier torsion result from either package on a folded or curved
mesh should be treated as unreliable.

## Comparison with the Abaqus built-up runs

`combo_1/J/model_runs/r5/calculate_J_r5.jl` used the same equation, J = T L / (G β), for the two-member
welded built-up section. For that model T was read from Abaqus as the applied torque and β from the
shear-center reference node rotation, exactly as done here with the shell model's reactions and the
section twist. A single S4R upright has not been run in Abaqus; if one is, its J should land near
1.33e-3 in⁴.
