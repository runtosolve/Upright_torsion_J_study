# make_torsion_html.jl — standalone, browser-interactive 3D page of the static twist study of the two-C welded
# upright (plotly.js): deformed shell under the end twist, welds, the fixed end (z = 0: all nodes fixed in X and Y,
# twist and translation restrained, warping free) and the applied rigid twist βo at z = L about the reference
# point (1.5, 1.5) (Abaqus kinematic coupling), warping free. Inset: twist of each C along the member.
# Inputs: torsion_deformed_r5.csv, torsion_summary_r5.csv (torsion_deformed_shape.jl), twist_profile_built_up.csv.
# Run:  julia --project=. make_torsion_html.jl
using DelimitedFiles, Statistics, Printf

d = readdlm(joinpath(@__DIR__, "torsion_deformed_r5.csv"), ','; skipstart = 1)
sm = readdlm(joinpath(@__DIR__, "torsion_summary_r5.csv"), ','; skipstart = 1)
L, weld_length, n_welds, T, J_eff, βo, Xc, Yc = Float64.(sm[1, :])
weld_spacing = 18.0; B = 3.0; D = 3.0 - 0.074; R = 0.199
sid = Int.(d[:, 1]); iid = Int.(d[:, 2]); jid = Int.(d[:, 3])
xyz = Float64.(d[:, 4:6]); u = Float64.(d[:, 7:9])
nn = maximum(iid); nz = maximum(jid)
sc = 0.30 / βo                                     # amplify so the end section is twisted by ≈ 0.30 rad (17°)
def = xyz .+ sc .* u
uin = hypot.(u[:, 1], u[:, 2]); cval = uin ./ maximum(uin)
idx = Dict((sid[k], iid[k], jid[k]) => k for k in eachindex(sid))
I = Int[]; J = Int[]; K = Int[]
for jj in 1:nz-1, ss in 1:2, ii in 1:nn-1
    a = idx[(ss, ii, jj)]; b = idx[(ss, ii, jj + 1)]; c = idx[(ss, ii + 1, jj + 1)]; e = idx[(ss, ii + 1, jj)]
    push!(I, a - 1); push!(J, b - 1); push!(K, c - 1); push!(I, a - 1); push!(J, c - 1); push!(K, e - 1)
end
ex = Any[]; ey = Any[]; ez = Any[]
function seg!(a, b)
    push!(ex, def[a, 3]); push!(ey, def[a, 1]); push!(ez, def[a, 2]); push!(ex, def[b, 3]); push!(ey, def[b, 1]); push!(ez, def[b, 2])
    push!(ex, nothing); push!(ey, nothing); push!(ez, nothing)
end
for jj in 1:nz, ss in 1:2, ii in 1:nn-1; seg!(idx[(ss, ii, jj)], idx[(ss, ii + 1, jj)]); end
for jj in 1:nz-1, ss in 1:2, ii in 1:nn; seg!(idx[(ss, ii, jj)], idx[(ss, ii, jj + 1)]); end
# welds (r5 layout: centers 1.5, 19.5, …, flush with the ends)
xp = xyz[:, 1] .- (sid .- 1) .* B
centers = range(weld_length / 2, L - weld_length / 2, Int(n_welds))
near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
corner = ((sid .== 1) .& near.(xp, B, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R))) .|
         ((sid .== 2) .& near.(xp, 0.0, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R)))
inweld = [any(near(z, c, weld_length / 2) for c in centers) for z in xyz[:, 3]]
w = findall(corner .& inweld .& (xyz[:, 3] .> 0) .& (xyz[:, 3] .< L))
# fixed end z = 0 (all nodes fixed in X, Y) and twisted end z = L (rigid rotation about (Xc, Yc))
e0 = findall(xyz[:, 3] .≈ 0.0); eL = findall(xyz[:, 3] .≈ L)
ox = Any[]; oy = Any[]; oz = Any[]                 # section outline at the fixed end
for ss in 1:2
    rows = [idx[(ss, ii, 1)] for ii in 1:nn]
    append!(ox, xyz[rows, 3]); append!(oy, xyz[rows, 1]); append!(oz, xyz[rows, 2]); push!(ox, nothing); push!(oy, nothing); push!(oz, nothing)
end
rarc = 3.6; zarc = L + 1.5
φs = range(-0.25π, 1.2π, 80)
ax_ = [zarc for φ in φs]; ay_ = [Xc + rarc * cos(φ) for φ in φs]; az_ = [Yc + rarc * sin(φ) for φ in φs]
φe = φs[end]
cx = [zarc]; cy = [Xc + rarc * cos(φe)]; cz = [Yc + rarc * sin(φe)]; cu = [-sin(φe)]; cv = [cos(φe)]
jm = (nz + 1) ÷ 2
sec_traces = String[]
for ss in 1:2
    rows = [idx[(ss, ii, jm)] for ii in 1:nn]
    push!(sec_traces, "{type:'scatter',x:[$(join(round.(xyz[rows, 1]; digits = 4), ','))],y:[$(join(round.(xyz[rows, 2]; digits = 4), ','))],mode:'lines',line:{color:'#8a8a8a',width:4},showlegend:false,hoverinfo:'skip',xaxis:'x',yaxis:'y'}")
    push!(sec_traces, "{type:'scatter',x:[$(join(round.(def[rows, 1]; digits = 4), ','))],y:[$(join(round.(def[rows, 2]; digits = 4), ','))],mode:'lines',line:{color:'#2a78d6',width:4},showlegend:false,hoverinfo:'skip',xaxis:'x',yaxis:'y'}")
end
js(v) = join((x === nothing ? "null" : string(round(x, digits = 5)) for x in v), ',')
rx_ = maximum(def[:, 1]) - minimum(def[:, 1]); ry_ = maximum(def[:, 2]) - minimum(def[:, 2])
ar = round.([L, rx_, ry_] ./ max(rx_, ry_, L) .* 1.9; digits = 3)

html = """
<!doctype html>
<html><head><meta charset="utf-8"><title>Built-up upright static twist (J_eff)</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>html,body{margin:0;height:100%;font-family:Helvetica,Arial,sans-serif;background:#fff}#plot{width:min(100vw,1600px);height:min(92vh,900px);margin:0 auto;overflow:hidden}</style></head>
<body><div id="plot"></div><div id="cam" style="position:fixed;left:8px;bottom:6px;font:11px Helvetica,Arial,sans-serif;color:#9a9a9a"></div>
<script>
const data = [
 {type:'mesh3d',x:[$(js(def[:, 3]))],y:[$(js(def[:, 1]))],z:[$(js(def[:, 2]))],i:[$(join(I, ','))],j:[$(join(J, ','))],k:[$(join(K, ','))],
  intensity:[$(js(cval))],colorscale:'Viridis',cmin:0,cmax:1,flatshading:true,lighting:{ambient:0.9,diffuse:0.2,specular:0.0},
  showscale:false,hoverinfo:'skip',showlegend:false,scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ex))],y:[$(js(ey))],z:[$(js(ez))],line:{color:'rgba(0,0,0,0.35)',width:1},hoverinfo:'skip',showlegend:false,scene:'scene'},
 {type:'scatter3d',mode:'markers',x:[$(js(def[w, 3]))],y:[$(js(def[w, 1]))],z:[$(js(def[w, 2]))],marker:{color:'#e34948',size:3.5},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ox))],y:[$(js(oy))],z:[$(js(oz))],line:{color:'#eda100',width:6},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ax_))],y:[$(js(ay_))],z:[$(js(az_))],line:{color:'#eda100',width:10},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'cone',x:[$(js(cx))],y:[$(js(cy))],z:[$(js(cz))],u:[$(js(zeros(length(cx))))],v:[$(js(cu))],w:[$(js(cv))],
  sizemode:'absolute',sizeref:1.6,anchor:'tip',colorscale:[[0,'#eda100'],[1,'#eda100']],showscale:false,showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'markers',x:[$(L)],y:[$(Xc)],z:[$(Yc)],marker:{color:'#eda100',size:7,symbol:'diamond'},showlegend:false,hoverinfo:'skip',scene:'scene'},
 $(join(sec_traces, ",\n "))
];
const layout = {
 title:{text:'Static twist of the two-C welded upright, L = $(Int(L)) in, $(Int(n_welds)) welds × $(Int(weld_length)) in at $(Int(weld_spacing)) in;  J<sub>eff</sub> = T L / (G β<sub>o</sub>) = $(round(J_eff, digits = 3)) in⁴<br>z = 0 fixed in X, Y (yellow outline: twist and translation restrained, warping free); rigid twist β<sub>o</sub> applied at z = L about (1.5, 1.5) (yellow arrow), warping free',x:0.02,y:0.97,xanchor:'left',yanchor:'top',font:{size:14}},
 scene:{domain:{x:[0,0.8],y:[0,1]},aspectmode:'manual',aspectratio:{x:$(ar[1]),y:$(ar[2]),z:$(ar[3])},xaxis:{title:{text:'Z (in)'},showbackground:false,showgrid:false,zeroline:false},yaxis:{visible:false},zaxis:{visible:false},
        camera:{projection:{type:'orthographic'},eye:{x:1.09,y:1.05,z:1.01},center:{x:0,y:0,z:0},up:{x:-0.399,y:-0.384,z:0.833}},dragmode:'orbit'},
 showlegend:false,
 xaxis:{domain:[0.83,0.99],title:{text:'Y (in), cross-aisle'},scaleanchor:'y',scaleratio:1,zeroline:false},
 yaxis:{domain:[0.35,0.65],title:{text:'X (in), downaisle'},zeroline:false},
 margin:{l:30,r:20,t:80,b:30},autosize:true,paper_bgcolor:'#fff'};
Plotly.newPlot('plot', data, layout, {responsive:true, displaylogo:false}).then(gd => {
  const show = c => { if(!c) return; const f = v => v.toFixed(2);
    document.getElementById('cam').textContent = 'camera eye (' + f(c.eye.x) + ', ' + f(c.eye.y) + ', ' + f(c.eye.z) + ')  up (' + f(c.up.x) + ', ' + f(c.up.y) + ', ' + f(c.up.z) + ')'; };
  show(gd.layout.scene.camera);
  gd.on('plotly_relayout', e => show(gd._fullLayout.scene.camera));
});
</script></body></html>
"""
tmp = joinpath(mktempdir(), "torsion_r5_plotly.html"); write(tmp, html)
out = joinpath(@__DIR__, "torsion_r5_plotly.html"); cp(tmp, out; force = true)
println("wrote ", out, "  (", round(filesize(out) / 1e6, digits = 1), " MB)")
