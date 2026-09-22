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

## Global flexural-torsional buckling, L = 44 in (`buckling_built_up.jl`, `plot_mode_wglmakie.jl`)

Elastic eigenbuckling of the welded pair under uniform axial compression, both ends pinned and warping
free (all end nodes fixed in X and Y, rotations and u_z free), 3 in welds at 18 in spacing placed
symmetrically (centers z = 4, 22, 40 in). Reference load: self-equilibrated tributary end forces, 1 kip
total. K φ = P_cr (−K_g) φ solved on the Ferrite-condensed dofs (Cᵀ K C) by shift-invert Arnoldi
(ArnoldiMethod.jl, sparse LU). A rigid in-plane cross-section constraint at every station (shared u, v, θ
for all nodes of both C's, warping free, weld ties kept for u_z/θx/θy) isolates the global modes; the
unconstrained model gives the true lowest modes.

| analysis | mode | P_cr (kips) | σ_cr (ksi) | character |
|---|---|---|---|---|
| global (rigid section) | 1 | **160.6** | 109.4 | flexural-torsional: y-translation with twist, one half-wave |
| global | 2 | 413.3 | 281.5 | second FT branch (translation and twist out of phase) |
| global | 3 | 606.6 | 413.1 | FT, two half-waves |
| global | 4 | 637.6 | 434.2 | pure flexure in x (perpendicular to the webs) |
| unconstrained shell | 1 | 101.3 | 69.0 | local plate buckling of the 2.9 in flats (rigid-section participation 1 %) |
| Euler, fully composite section | P_ey / P_ex | 351.5 / 772.4 | | bending about the horizontal / vertical axis |

The built-up section is symmetric about its horizontal axis, so the shear center is offset horizontally
and, per Timoshenko, the vertical translation couples with twist while the horizontal translation buckles
as pure flexure — exactly the mode pattern found. The FT load (160.6 kips) is 46 % of the composite
Euler load in the same direction because the twist about the offset shear center is resisted only by the
intermittently welded pair (J_eff ≈ 0.77 in⁴) and its warping. Pure x-flexure at 637.6 kips is 17 % below
the composite Euler value, the loss of composite action between welds. At 44 in the member is governed
by local buckling (101 kips, 69 ksi), so the global FT load is well above the local one.

Mode shape: `mode_global_1_wglmakie.html` (interactive WGLMakie scene, standalone) and `mode_global_1.png`
(CairoMakie), both from `plot_mode_wglmakie.jl`; mode data in `mode_global_1.csv`, `mode_global_2.csv`,
`mode_all_1.csv`; loads in `buckling_results.csv`.

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
