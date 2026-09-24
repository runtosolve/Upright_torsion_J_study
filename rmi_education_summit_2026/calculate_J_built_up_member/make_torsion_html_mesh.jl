# make_torsion_html_mesh.jl — browser-interactive 3D page (plotly.js) of the static twist of the PERFORATED welded pair
# from the Gmsh mixed mesh (torsion_perforated_{nodes,cells,summary}.csv). Same content and style as
# make_torsion_html.jl: deformed shell coloured by in-plane displacement, element edges, welds (red), fixed end
# outline (yellow), circular torque arrow at z = L (yellow), mid-length section inset.
# Run:  julia --project=. make_torsion_html_mesh.jl
using DelimitedFiles, Statistics, Printf
include(joinpath(@__DIR__, "..", "calculate_J_single_member", "calculate_J_single_member.jl"))

nodes = readdlm(joinpath(@__DIR__, "torsion_perforated_nodes.csv"), ','; skipstart = 1)
xyz = Float64.(nodes[:, 1:3]); u = Float64.(nodes[:, 4:6])
cells = [parse.(Int, split(strip(l), ',')) for l in readlines(joinpath(@__DIR__, "torsion_perforated_cells.csv"))]
sm = readdlm(joinpath(@__DIR__, "torsion_perforated_summary.csv"), ','; skipstart = 1)
L, T, J_eff, βo_, Xc, Yc, nq, nt = Float64.(sm[1, :])
weld_length = 3.0; weld_spacing = 18.0; B = 3.0; D = 3.0 - 0.074; R = 0.199
n_welds = 7; centers = range(weld_length / 2, L - weld_length / 2, n_welds)
sc = 0.30 / βo_; def = xyz .+ sc .* u
uin = hypot.(u[:, 1], u[:, 2]); cval = uin ./ maximum(uin)
# triangles (quads split), plot coords (Z, X, Y)
I = Int[]; J = Int[]; K = Int[]
for c in cells
    push!(I, c[1] - 1); push!(J, c[2] - 1); push!(K, c[3] - 1)
    length(c) == 4 && (push!(I, c[1] - 1); push!(J, c[3] - 1); push!(K, c[4] - 1))
end
# element edges; boundary edges (used once: hole outlines, ends, lip tips) are drawn, interior edges are not (file size)
edge_count = Dict{Tuple{Int,Int},Int}()
for c in cells, k in 1:length(c)
    a, b = c[k], c[mod1(k + 1, length(c))]; e = (min(a, b), max(a, b)); edge_count[e] = get(edge_count, e, 0) + 1
end
edges = keys(edge_count)
bedges = [e for (e, n) in edge_count if n == 1]
ex = Any[]; ey = Any[]; ez = Any[]
for (a, b) in bedges
    push!(ex, def[a, 3]); push!(ey, def[a, 1]); push!(ez, def[a, 2]); push!(ex, def[b, 3]); push!(ey, def[b, 1]); push!(ez, def[b, 2])
    push!(ex, nothing); push!(ey, nothing); push!(ez, nothing)
end
# welds: tie nodes (C1 corners, C2 web corners) within the weld windows
xp = xyz[:, 1] .- (xyz[:, 1] .> B + 0.5) .* B; s1 = xyz[:, 1] .<= B + 0.5
near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
corner = (s1 .& near.(xp, B, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R))) .| (.!s1 .& near.(xp, 0.0, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R)))
inweld = [any(near(z, c, weld_length / 2) for c in centers) for z in xyz[:, 3]]
w = findall(corner .& inweld .& (xyz[:, 3] .> 0) .& (xyz[:, 3] .< L))
# fixed-end outline: mesh edges with both nodes at z = 0
ox = Any[]; oy = Any[]; oz = Any[]
for (a, b) in edges
    (xyz[a, 3] < 1e-6 && xyz[b, 3] < 1e-6) || continue
    push!(ox, 0.0); push!(oy, xyz[a, 1]); push!(oz, xyz[a, 2]); push!(ox, 0.0); push!(oy, xyz[b, 1]); push!(oz, xyz[b, 2]); push!(ox, nothing); push!(oy, nothing); push!(oz, nothing)
end
# torque arrow (arc about the reference point at z = L + 1.5)
rarc = 3.6; zarc = L + 1.5; φs = range(-0.25π, 1.2π, 80)
ax_ = [zarc for φ in φs]; ay_ = [Xc + rarc * cos(φ) for φ in φs]; az_ = [Yc + rarc * sin(φ) for φ in φs]; φe = φs[end]
cx = [zarc]; cy = [Xc + rarc * cos(φe)]; cz = [Yc + rarc * sin(φe)]; cu = [-sin(φe)]; cv = [cos(φe)]
# mid-length section: nearest mesh node to each centerline point at z = L/2 (rigid section motion → good approximation)
X1, Y1 = centerline(shape); sec_traces = String[]
for s in 1:2
    xs0 = X1 .+ (s - 1) * B; ys0 = Y1
    rows = [argmin((xyz[:, 1] .- xs0[i]) .^ 2 .+ (xyz[:, 2] .- ys0[i]) .^ 2 .+ (xyz[:, 3] .- L / 2) .^ 2) for i in eachindex(xs0)]
    push!(sec_traces, "{type:'scatter',x:[$(join(round.(xs0; digits = 4), ','))],y:[$(join(round.(ys0; digits = 4), ','))],mode:'lines',line:{color:'#8a8a8a',width:4},showlegend:false,hoverinfo:'skip',xaxis:'x',yaxis:'y'}")
    push!(sec_traces, "{type:'scatter',x:[$(join(round.(xs0 .+ sc .* u[rows, 1]; digits = 4), ','))],y:[$(join(round.(ys0 .+ sc .* u[rows, 2]; digits = 4), ','))],mode:'lines',line:{color:'#2a78d6',width:4},showlegend:false,hoverinfo:'skip',xaxis:'x',yaxis:'y'}")
end
js(v) = join((x === nothing ? "null" : string(round(x, digits = 3)) for x in v), ',')
rx_ = maximum(def[:, 1]) - minimum(def[:, 1]); ry_ = maximum(def[:, 2]) - minimum(def[:, 2])
ar = round.([L, rx_, ry_] ./ max(rx_, ry_, L) .* 1.9; digits = 3)

html = """
<!doctype html>
<html><head><meta charset="utf-8"><title>Perforated built-up upright static twist (J_eff)</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>html,body{margin:0;height:100%;font-family:Helvetica,Arial,sans-serif;background:#fff}#plot{width:min(100vw,1600px);height:min(92vh,900px);margin:0 auto;overflow:hidden}</style></head>
<body><div id="plot"></div><div id="cam" style="position:fixed;left:8px;bottom:6px;font:11px Helvetica,Arial,sans-serif;color:#9a9a9a"></div>
<script>
const data = [
 {type:'mesh3d',x:[$(js(def[:, 3]))],y:[$(js(def[:, 1]))],z:[$(js(def[:, 2]))],i:[$(join(I, ','))],j:[$(join(J, ','))],k:[$(join(K, ','))],
  intensity:[$(js(cval))],colorscale:'Viridis',cmin:0,cmax:1,flatshading:true,lighting:{ambient:0.9,diffuse:0.2,specular:0.0},showscale:false,hoverinfo:'skip',showlegend:false,scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ex))],y:[$(js(ey))],z:[$(js(ez))],line:{color:'rgba(0,0,0,0.6)',width:2},hoverinfo:'skip',showlegend:false,scene:'scene'},
 {type:'scatter3d',mode:'markers',x:[$(js(def[w, 3]))],y:[$(js(def[w, 1]))],z:[$(js(def[w, 2]))],marker:{color:'#e34948',size:3.5},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ox))],y:[$(js(oy))],z:[$(js(oz))],line:{color:'#eda100',width:6},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ax_))],y:[$(js(ay_))],z:[$(js(az_))],line:{color:'#eda100',width:10},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'cone',x:[$(js(cx))],y:[$(js(cy))],z:[$(js(cz))],u:[$(js(zeros(1)))],v:[$(js(cu))],w:[$(js(cv))],sizemode:'absolute',sizeref:1.6,anchor:'tip',colorscale:[[0,'#eda100'],[1,'#eda100']],showscale:false,showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'markers',x:[$(L)],y:[$(Xc)],z:[$(Yc)],marker:{color:'#eda100',size:7,symbol:'diamond'},showlegend:false,hoverinfo:'skip',scene:'scene'},
 $(join(sec_traces, ",\n "))
];
const layout = {
 title:{text:'Static twist of the PERFORATED two-C welded upright, L = $(Int(L)) in, $(n_welds) welds × $(Int(weld_length)) in at $(Int(weld_spacing)) in;  J<sub>eff</sub> = T L / (G β<sub>o</sub>) = $(round(J_eff, digits = 3)) in⁴ (gross section 0.737)<br>teardrop web holes and square flange holes (Gmsh mesh, $(Int(nq)) quads + $(Int(nt)) triangles); z = 0 fixed in X, Y (yellow outline, warping free); rigid twist β<sub>o</sub> at z = L about (1.5, 1.5) (yellow arrow), warping free',x:0.02,y:0.97,xanchor:'left',yanchor:'top',font:{size:14}},
 scene:{domain:{x:[0,0.8],y:[0,1]},aspectmode:'manual',aspectratio:{x:$(ar[1]),y:$(ar[2]),z:$(ar[3])},xaxis:{title:{text:'Z (in)'},showbackground:false,showgrid:false,zeroline:false},yaxis:{visible:false},zaxis:{visible:false},
        camera:{projection:{type:'orthographic'},eye:{x:1.09,y:1.05,z:1.01},center:{x:0,y:0,z:0},up:{x:-0.399,y:-0.384,z:0.833}},dragmode:'orbit'},
 showlegend:false,
 xaxis:{domain:[0.83,0.99],title:{text:'X (in)'},scaleanchor:'y',scaleratio:1,zeroline:false},
 yaxis:{domain:[0.35,0.65],title:{text:'Y (in)'},zeroline:false},
 margin:{l:30,r:20,t:80,b:30},autosize:true,paper_bgcolor:'#fff'};
Plotly.newPlot('plot', data, layout, {responsive:true, displaylogo:false}).then(gd => {
  const show = c => { if(!c) return; const f = v => v.toFixed(2);
    document.getElementById('cam').textContent = 'camera eye (' + f(c.eye.x) + ', ' + f(c.eye.y) + ', ' + f(c.eye.z) + ')  up (' + f(c.up.x) + ', ' + f(c.up.y) + ', ' + f(c.up.z) + ')'; };
  show(gd.layout.scene.camera); gd.on('plotly_relayout', e => show(gd._fullLayout.scene.camera));
});
</script></body></html>
"""
tmp = joinpath(mktempdir(), "torsion_perforated_plotly.html"); write(tmp, html)
out = joinpath(@__DIR__, "torsion_perforated_plotly.html"); cp(tmp, out; force = true)
println("wrote ", out, "  (", round(filesize(out) / 1e6, digits = 1), " MB)")
