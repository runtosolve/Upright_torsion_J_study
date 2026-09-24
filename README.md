# Upright torsion constant (J) study — Ferrite.jl

Saint-Venant torsion constant studies of a cold-formed steel upright made of two 3 × 3 × 0.074 in lipped C
sections welded back-to-lip with intermittent 3 in welds, using [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/)
shell models ([QuadShellFiniteElement.jl](https://github.com/runtosolve/QuadShellFiniteElement.jl)) and the static
twist method of Moen (2008), §4.2.7.3.2.3: J = T L / (G βo) with warping-free ends.

**Hosted page:** https://runtosolve.github.io/Upright_torsion_J_study/ — interactive 3D deformed shape of the
welded pair under the end twist (7 × 3 in welds at 18 in, L = 111 in), showing the welds, the fixed end
(z = 0: X and Y fixed, twist restrained, warping free) and the applied rigid twist at z = L (warping free).
Rotate, zoom and pan in the browser; the bottom-left corner shows the live camera vectors.

**Perforated upright:** https://runtosolve.github.io/Upright_torsion_J_study/perforated.html — the same twist test on the
welded pair with the teardrop web holes and square flange holes (Gmsh mixed quad/tri mesh), J_eff = 0.572 in⁴ vs 0.737 gross.

- `rmi_education_summit_2026/calculate_J_single_member/` — J of one C: static twist method vs. an exact 2D
  Saint-Venant solution; documents the drilling-dof finding that led to the Hughes–Brezzi option in
  QuadShellFiniteElement.jl / TriShellFiniteElement.jl.
- `rmi_education_summit_2026/calculate_J_built_up_member/` — effective J of the welded pair (J_eff = 0.768 in⁴
  for the r5 layout), weld spacing study vs. the Tlumak equation, eigenbuckling scripts, and the page generators
  (`make_torsion_html.jl`, `torsion_deformed_shape.jl`).
- `rmi_education_summit_2026/calculate_J_built_up_member/perforated_J_*.jl` — the same J study with the upright's
  teardrop/square perforation pattern: structured mesh with removed cells, and a Gmsh mixed quad/tri mesh with the true
  hole outlines assembled with QuadShellFiniteElement + TriShellFiniteElement (single C −9 %, welded pair −22 %).
- `docs/` — GitHub Pages site (`index.html` = torsion page) and summary figures.

The companion flexural-torsional buckling study is hosted at https://runtosolve.github.io/UprightLTB/
(repository [runtosolve/UprightLTB](https://github.com/runtosolve/UprightLTB)).

Run scripts from their folder with `julia --project=.` (Project/Manifest reference packages developed under `~/.julia/dev`).
