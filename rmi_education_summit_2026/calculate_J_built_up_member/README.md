# Effective J of the two-C welded OneRack column with Ferrite.jl — Moen (2008) §4.2.7.3.2.3 method

Companion to `../calculate_J_single_member/` (which it includes for the geometry, material, mesh and
helper functions). Same Julia environment (`julia --project=.`).

## Model (mirrors `combo_1/J/combo_1_J_r5.jl` / `combo_1_J_r5.inp`)

- Two 3 × 3 × 0.074 in lipped C uprights (0.75 in lips, R = 0.199 in centerline corners, L = 111 in),
  SHAPE_2 = SHAPE_1 translated by B = 3.0 in in X, so the lip/flange corners of C1 bear on the web of C2.
  Mesh as in Abaqus: 4 elements per flat, 5 per corner, 0.5 in along; 109,716 dofs.
- Welds: seven 3 in welds at z = 1.5, 19.5, …, 109.5 in. Each weld ties the C1 corner nodes
  (|X − B| ≤ R, |Y − D| ≤ R for the top, |Y| ≤ R for the bottom) and the C2 web-corner nodes
  (|X| ≤ R, same Y window) within |Δz| ≤ 1.5 in into one rigid body — the Abaqus `*Rigid Body, tie`
  definitions — as Ferrite affine constraints to one master node per weld (926 constraints). The z = 0
  and z = L stations are left out of the weld sets because they carry the Dirichlet conditions.
- Ends (Fig. 4.41 / the Abaqus step): z = 0 all nodes of both shapes fixed in X, Y and the C2 web
  mid-height node fixed in Z; z = L all nodes rotated rigidly by βo about the Abaqus reference point
  (1.5, 1.5) (kinematic coupling to node 10000 in dofs 1, 2); warping free at both ends.
- J_eff = T L / (G βo) with T the resultant moment of the z = L reactions (Eq. 4.13). QuadShellFiniteElement
  with Hughes–Brezzi drilling (see the single-member study for why).

## Results

| case | J_eff (in⁴) | J_eff / J_box |
|---|---|---|
| Abaqus weld layout, 7 × 3 in welds | **0.7678** | 0.414 |
| same, fine mesh (8/flat, 10/corner, 444 along; 432,540 dofs) | 0.7576 | 0.409 |
| same + bilateral normal tie of the C1 lips to the C2 web between welds (contact stand-in) | 0.7678 | 0.414 |
| same, twist applied about the C1 shear center instead of (1.5, 1.5) | 0.7678 | 0.414 |
| same, Cs = 0 | 0.7683 | 0.414 |
| fully welded, 37 back-to-back 3 in welds | 2.0665 | 1.115 |
| no welds | 0.002651 | 0.0014 |

References: one C exact Saint-Venant J_C = 1.3349e-3 in⁴ (2 J_C = 2.670e-3); Bredt closed cell
(C1 web + flanges + C2 web, 2.926 in square centerline, t = 0.074) J_box = 1.854 in⁴; Abaqus
`combo_1_J_r5` as read in `calculate_J_r5.jl` (S4R, rigid-body welds, friction contact, NLGEOM with
stabilization): J = 0.293 in⁴.

Checks: no welds returns exactly twice the single-C value; every case has zero net in-plane end force
(|ΣF|/T ≈ 1e-10), so the result does not depend on the point the twist is applied about; the fully
welded member is a closed section (1.11 × Bredt — the 3 in rigid weld blocks along the corners and the
doubled wall where the C1 lips lie on the C2 web add stiffness beyond the thin-wall formula); mesh
refinement changes J_eff by 1.3 %.

## What J_eff means here

The twist is not uniform along the intermittently welded member (linear-fit R² = 0.9992 instead of
1.0000). C2 twists at almost the uniform rate everywhere (0.93 to 1.15 × βo/L); C1 twists at 1.17 × βo/L
between welds and 0.42 × in the 3 in welded zones, where the C1 web + flanges + C2 web form a closed
cell, so the relative twist of the two shapes oscillates by about ±0.012 βo from weld to weld. Because
the two shapes are forced to rotate about a common axis while their own shear centers sit 3 in apart,
each weld-to-weld segment also bends and warps, and that, not Saint-Venant shear, carries the torque:
a Saint-Venant-only series estimate, L / Σ(Lᵢ/Jᵢ) with J_box in the welded lengths and 2 J_C between
them, gives only 3.3e-3 in⁴, while the model gives 0.77 in⁴. More than 99 % of the member's torsional
stiffness therefore comes from the E Cw-type terms (warping restraint at the welds and coupled bending
of the two C's between them), not from G J. J_eff is an equivalent member property that depends on weld
spacing, weld length and member length, and it is what a beam element with only G J would need to
reproduce the end twist. The arithmetic weld-weighted estimate in `calculate_J_r5.jl`
(L_weld/L × J_box = 0.351 in⁴) has no mechanical basis but happens to land in the same decade.

## Comparison with the Abaqus reading

The Ferrite value (0.77 in⁴) is 2.6 × the 0.293 in⁴ derived from the Abaqus odb. The model
ingredients are the same (geometry, mesh density, rigid-body welds, end conditions); the Abaqus
difference is S4R, NLGEOM with automatic stabilization (`stabilize = 0.002`), general contact with
friction, and the fact that both T and β in `calculate_J_r5.jl` were read at the first increment of a
nonlinear stabilized step (T = 0.001 × 2.9688e-3/2.5 = 1.19e-6 lbf-in, β = 3.97e-11 rad). The linear
shell model here shows contact (as a normal tie) changes nothing, so the discrepancy most likely sits in
the Abaqus reading: stabilization damping and the load/rotation frame pairing at a tiny first increment.
Recommended check: rerun `combo_1_J_r5` as a linear `*STATIC` step without `stabilize` and NLGEOM (or a
`*STATIC, PERTURBATION`), and read RM3 of node 10000 and UR3 from the same frame; the effective J should
then be close to 0.76 in⁴.

## Weld spacing study vs. the proposed prediction equation (`weld_spacing_study.jl`)

Proposed equation: J_combined = J_tube · (w_length / w_spacing) + J_back + J_front · (1 − w_length / w_spacing),
with J_tube = 1.854 in⁴ (Bredt, C1 web + flanges + C2 web) and J_front = J_back = J_C = 1.335e-3 in⁴.
Welds equally spaced from end to end (first and last flush with the ends, as in the r5 model),
w_spacing = (L − w_length)/(n − 1). Results in `weld_spacing_results.csv`, figure `weld_spacing_summary.png`.

| w_length (in) | w_spacing (in) | n welds | J_eff Ferrite (in⁴) | equation (in⁴) | Ferrite / equation |
|---|---|---|---|---|---|
| 3 | 108 | 2 | 0.0309 | 0.0541 | 0.57 |
| 3 | 54 | 3 | 0.1108 | 0.1056 | 1.05 |
| 3 | 36 | 4 | 0.2350 | 0.1570 | 1.50 |
| 3 | 27 | 5 | 0.3941 | 0.2085 | 1.89 |
| 3 | 18 | 7 | 0.7678 | 0.3114 | 2.47 |
| 3 | 12 | 10 | 1.2788 | 0.4658 | 2.75 |
| 3 | 9 | 13 | 1.6163 | 0.6202 | 2.61 |
| 3 | 6 | 19 | 1.9230 | 0.9289 | 2.07 |
| 3 | 3 (continuous) | 37 | 2.0665 | 1.8551 | 1.11 |
| 6 | 105 … 7 | 2 … 16 | 0.0367 … 2.064 | 0.1085 … 1.590 | 0.34 → 1.89 → 1.30 |
| 1.5 | 109.5 … 3.04 | 2 … 37 | 0.0282 … 1.934 | 0.0280 … 0.916 | 1.0 → 4.15 → 2.11 |
| 3 | 18, L = 57 / 219 | 4 / 13 | 0.7814 / 0.7608 | 0.3114 | 2.51 / 2.44 |

Findings:

- The equation and the shell model do not share a governing variable. The equation scales J_eff with
  the welded fraction w_length/w_spacing; the model shows J_eff is set almost entirely by the weld
  **spacing** and only weakly by the weld length. At an 18 in spacing, 1.5 in, 3 in and 6 in welds give
  0.58, 0.77 and 0.86 in⁴ (the equation says 0.155, 0.311 and 0.53). Doubling the weld length raises
  J_eff by ~15 %, not by 2 ×.
- Consequently the ratio Ferrite/equation is not a constant: for 3 in welds it runs from 0.57 (two end
  welds) through a peak of 2.75 at 12 in spacing back to 1.11 for the continuous weld, and for 1.5 in
  welds it peaks at 4.15. The equation is unconservative only for very few, long welds (2 welds of 6 in:
  0.34) and conservative by 2 to 4 × in the practical range of 9 to 36 in spacings.
- J_eff is nearly independent of member length (0.78, 0.77, 0.76 in⁴ for L = 57, 111, 219 in at 18 in
  spacing), so it is a property of the weld pattern, not of the member.
- Mechanism: each weld is a discrete restraint forcing the two C's, whose shear centers are 3 in apart,
  to twist together. Between welds the segments resist the relative twist by bending and warping, whose
  stiffness scales roughly with 1/spacing², so J_eff/J_tube follows approximately 1/(1 + (w_spacing/s₀)^p)
  with s₀ ≈ 15 in and p ≈ 2.1 for 3 in welds on this section (an empirical fit to series A; the
  continuous-weld limit 1.11 J_tube exceeds Bredt because of the doubled wall and rigid weld blocks).
  A predictor built on (w_spacing / section depth) would track the model far better than one built on
  w_length / w_spacing.

## Perforated uprights (`perforated_J_study.jl`, `perforated_J_gmsh.jl`)

The same static twist study with the upright's perforation pattern (14 ga strip layout drawing): teardrop holes
(Ø0.719 in circle + Ø0.375 in lobe pointing along the member, ≈ 0.90 in long) in two rows 0.797 in from each
web face (1.406 in apart) every 2.0 in; 0.562 × 0.562 in square holes (R0.02) centred 0.856 in from the web face on
both flanges every 2.0 in, placed 1.0 in (half a pitch) from the teardrops along the member (assumed stagger; the
drawing's 0.9215 / 1.517 in dimensions fix the true phase). Holes remove 3.9 % of the shell area.

Two meshes:

- `perforated_J_study.jl`: structured extruded quad mesh (0.21 in around, 0.2 in along), holes as removed cells
  whose centres lie inside the outlines (pixelated edges).
- `perforated_J_gmsh.jl`: the developed strip is drawn in Gmsh (OCC) with the true hole outlines, fragmented by the
  bend lines at every centerline vertex, meshed with recombination (h = 0.2 in; quads plus triangles at the hole
  boundaries), folded onto the C, and assembled through Ferrite SubDofHandlers with QuadShellFiniteElement.jl for the
  quads and TriShellFiniteElement.jl for the triangles (both Hughes–Brezzi drilling). `gmsh_strip_mesh.png` shows the
  strip mesh. The unperforated Gmsh mesh reproduces the structured-mesh J (1.3287e-3 vs 1.3271e-3 in⁴).

| L = 111 in | no holes | perforated | ratio |
|---|---|---|---|
| single C, J, Gmsh mesh (44,947 quads + 3,172 tris) | 1.3287e-3 in⁴ | 1.2034e-3 in⁴ | 0.906 |
| single C, J, structured mesh | 1.3271e-3 | 1.2146e-3 | 0.915 |
| welded pair r5 (7 × 3 in welds at 18 in), J_eff, Gmsh mesh (89,894 quads + 6,344 tris) | 0.737 in⁴ | 0.572 in⁴ | 0.777 |
| welded pair r5, J_eff, structured mesh | 0.753 | 0.592 | 0.785 |

Weld spacing sweep with perforations (`perforated_weld_spacing_study.jl`, Gmsh mesh, 3 in welds, L = 111 in;
`weld_spacing_comparison_perforated.png`, `perforated_weld_spacing_results.csv`):

| w_spacing (in) | 108 | 54 | 36 | 27 | 18 | 12 | 9 | 6 | 3 (continuous) |
|---|---|---|---|---|---|---|---|---|---|
| J_eff gross (in⁴) | 0.031 | 0.109 | 0.229 | 0.382 | 0.737 | 1.226 | 1.557 | 1.876 | 2.066 |
| J_eff perforated (in⁴) | 0.027 | 0.095 | 0.194 | 0.315 | 0.572 | 0.882 | 1.067 | 1.231 | 1.342 |
| perforated / gross | 0.88 | 0.87 | 0.85 | 0.82 | 0.78 | 0.72 | 0.69 | 0.66 | 0.65 |

The hole penalty grows as the welds get closer: 12 % at 108 in spacing, 22 % at the r5 spacing of 18 in, and 35 % for
the continuous weld, where the perforated closed cell keeps only 0.65 of the gross J_tube-type stiffness because the
teardrops interrupt the web shear flow of the cell.

The two meshes agree within 1 % on the reduction. For one C the holes cut J by 9 to 10 %, about 2.4 times the
removed area fraction: each hole adds free edges and disturbs the Saint-Venant shear flow over a zone longer than
the hole (compare the free-edge loss of ≈ 0.63 t⁴/3 per edge from the single-member study). For the welded pair the
loss is 22 %, larger because the closed-cell shear flow through the welds runs through the perforated webs and the
flange squares sit in the load path between the welds and the C1 lips. Results in `perforated_J_results.csv` and
`perforated_J_gmsh_results.csv`.

## Global flexural-torsional buckling (`buckling_built_up.jl`, `ft_analytical_check.jl`, `make_mode_html.jl`)

Elastic eigenbuckling of the welded pair under uniform axial compression, both ends pinned and warping free (all
end nodes fixed in X and Y, rotations and u_z free), 3 in welds at 18 in spacing placed symmetrically about
mid-length. Reference load: self-equilibrated tributary end forces, 1 kip total. K φ = P_cr (−K_g) φ solved on the
Ferrite-condensed dofs (Cᵀ K C) by shift-invert Arnoldi (ArnoldiMethod.jl, sparse LU). Two models:

- **unconstrained shell** — the true lowest modes, local plate buckling included;
- **rigid-section (global) model** — every station's in-plane motion reduced to one (u, v, θ), warping free, weld
  ties kept for u_z/θx/θy. It isolates the global modes where local modes come first (short members). Its K_g is
  built from the uniform axial reference stress −P/A: the constrained prebuckling solve blocks Poisson expansion
  and would otherwise carry a spurious transverse compression ν σ_z that destabilizes twist (an early version of
  this study reported 160.6 kips for L = 44 in from that contaminated K_g; the corrected value is 298.9 kips).
  The same constraint also suppresses transverse membrane strain, so the plates act with E/(1 − ν²): the model
  overstates pure y-flexure by 3 % at L = 44 in (362.6 vs Euler 351.5 kips) and by ~10 % at L = 120 in (54.6 vs
  47.3 kips). Its global loads are therefore upper bounds by roughly that margin.

| L (in) | welds | analysis | mode | P_cr (kips) | character |
|---|---|---|---|---|---|
| 44 | 3 × 3 in at z = 4, 22, 40 | rigid-section | 1 | **298.9** | flexural-torsional: y-translation with twist (θ·r/v = 0.27) |
| 44 | | rigid-section | 2 | 686.5 | pure flexure in x |
| 44 | | rigid-section, pure torsion about the static shear center | | 778.4 | |
| 44 | | rigid-section, pure y-flexure | | 362.6 | Euler 351.5 |
| 44 | | unconstrained shell | 1 | 101.3 | local plate buckling of the 2.9 in flats |
| 120 | 7 × 3 in at z = 6 … 114 | unconstrained shell | 1 | **49.2** | global: weak-axis (y) flexure with slight twist (θ·r/v = 0.07), the coupled FT root |
| 120 | | rigid-section | 1 | 52.7 | same mode (E/(1 − ν²) artifact) |
| 120 | | unconstrained shell | 2 | 103.0 | local |
| 120 | | Euler, composite I_x | | 47.3 | |

Analytical comparison (Timoshenko, section symmetric about the horizontal axis so bending about x couples with
torsion: P_FT = [(P_ey + P_t) − √((P_ey + P_t)² − 4 β P_ey P_t)] / (2β), β = 1 − (x_o/r_o)²). The shear center
by statics on the rigid-section model (zero twist under a mid-length transverse force, welds included) is at
x = 0.39 in, x_o = −2.36 in from the composite centroid at x = 2.74 in (`section_centers.png`).

| L (in) | Timoshenko with shell P_ey and shell P_t | Timoshenko with Euler P_ey, P_t = G J_eff / r_o² (J_eff = 0.77 in⁴, C_w = 0) | shell |
|---|---|---|---|
| 44 | 280.3 kips (0.94 × shell) | 277.4 kips (P_t = 819 kips) | 298.9 (rigid-section) |
| 120 | 52.4 kips (0.99 × shell) | 45.8 kips (P_t = 798 kips) | 52.7 (rigid-section) / 49.2 (unconstrained) |

The shell's pure-torsional buckling load about the static shear center (778 kips at 44 in, 736 at 120 in) corresponds
to a torsional stiffness P_t r_o² of 8,280 and 8,030 kip-in², i.e. 92 to 95 % of G J_eff = 8,712 kip-in² with no
additional warping contribution: for buckling of the welded pair the static-twist J_eff is the right torsional
constant and the effective C_w is negligible. With a colleague's beam properties for this section (A = 1.219 in²,
r_x = 1.297 in, J = 0.319 in⁴ by the Tlumak equation, C_w = 4.358 in⁶, r_o = 3.272 in, x_o as above) the formula
gives P_cre = 202 kips at L = 44 in (P_t = 399 kips, stiffness 4,270 kip-in²), 68 % of the shell value; the
difference is the torsion constant, 0.319 vs 0.77 in⁴. At 120 in the buckling is flexure-dominated and every
approach lands within ±10 %.

Validation of the rigid-section buckling analysis on a single C (`single_c_buckling_check.jl`,
`single_c_Cw_check.jl`, `single_c_energy_check.jl`), where J, C_w = 2.385 in⁶ (cutwp, reproduced exactly by an
independent sectorial-coordinate calculation) and the shear center are known: pure torsion about the shear center
36.3 vs 34.0 kips theory at 44 in (1.07) and 6.17 vs 5.76 at 120 in (1.07); pure y-flexure 169.8 vs 175.7 (0.97);
pure x-flexure 145.3 vs 137.7 (1.06); mesh-independent. An energy decomposition of the torsional mode gives the
geometric term at 1.000 × P r_o² ∫θ′² and the elastic term at 1.07 (44 in) / 1.02 (400 in) × (G J ∫θ′² + E C_w ∫θ″²),
so the only bias is the E/(1 − ν²) stiffening of the rigid-section constraint. A static warping-fixed twist of the
shell gives C_w = 2.35 in⁶ with free sections (2.56 with rigid sections).

Pages: `mode_global_1_plotly.html` (L = 44 in, rigid-section mode) and `mode_all_1_L120_plotly.html` (L = 120 in,
unconstrained mode), generated by `make_mode_html.jl <mode csv> <results csv> <global|all>`; WGLMakie/CairoMakie
versions from `plot_mode_wglmakie.jl`. Results in `buckling_results.csv`, `buckling_results_L120.csv`,
`ft_analytical_check.csv`, `ft_analytical_check_L120.csv`.

## Files

| file | purpose |
|---|---|
| `calculate_J_built_up_member.jl` | model, welds as affine rigid ties, cases above; writes `J_built_up_results.csv`, `twist_profile_built_up.csv`, `run_log.txt` |
| `make_figure_built_up.jl` | `J_built_up_summary.png`: twist of each C along the member with welded zones shaded, and the J_eff comparison |
| `weld_spacing_study.jl` | weld spacing / length / member length sweep vs. the proposed equation; writes `weld_spacing_results.csv`, `weld_spacing_log.txt` |
| `make_figure_weld_spacing.jl` | `weld_spacing_summary.png` |
| `make_figure_spacing_comparison.jl` | `weld_spacing_comparison.png`: Ferrite J_eff vs. the Tlumak equation vs. weld spacing, 3 in welds |
| `buckling_built_up.jl` | eigenbuckling at L = 44 in (global rigid-section and unconstrained); writes `buckling_results.csv`, `mode_*.csv`, `buckling_log.txt` |
| `plot_mode_wglmakie.jl` | 3D mode shape: `mode_global_1_wglmakie.html` (WGLMakie) and `mode_global_1.png` |
| `ft_analytical_check.jl` | Timoshenko FT formula vs. the shell: constrained pure-torsion / pure-flexure eigen loads, static shear center; writes `ft_analytical_check.csv` (`_L120` for 120 in) |
| `perforated_weld_spacing_study.jl`, `make_figure_spacing_perforated.jl` | weld spacing sweep with and without perforations (Gmsh mesh); `weld_spacing_comparison_perforated.png` |
| `torsion_deformed_shape_perforated.jl`, `make_torsion_html_mesh.jl` | interactive page of the perforated pair twist from the Gmsh mesh (`torsion_perforated_plotly.html`) |
| `perforated_J_study.jl`, `perforated_J_gmsh.jl` | perforated single C and welded pair J (structured mesh with removed cells; Gmsh mixed quad/tri mesh with true hole outlines, TriShell + QuadShell assembly); `gmsh_strip_mesh.png` |
| `single_c_buckling_check.jl`, `single_c_Cw_check.jl`, `single_c_energy_check.jl` | validation of the rigid-section buckling machinery on a single C against classical P_t, Euler, C_w and an energy decomposition |
