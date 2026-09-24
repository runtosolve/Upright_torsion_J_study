# make_mode_html.jl — standalone, browser-interactive 3D buckling mode shape (plotly.js, no server needed) from
# mode_*.csv written by buckling_built_up.jl. WGLMakie's camera runs in Julia, so its static export is not
# interactive; this page is for hosting (GitHub Pages).
# Run:  julia --project=. make_mode_html.jl [mode_global_1.csv]
using DelimitedFiles, Statistics, Printf

path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "mode_global_1.csv")
respath = length(ARGS) >= 2 ? ARGS[2] : joinpath(@__DIR__, "buckling_results.csv")
res = readdlm(respath, ','; skipstart = 1)
analysis = length(ARGS) >= 3 ? ARGS[3] : "global"          # which row of buckling_results*.csv gives P_cre: "global" or "all"
row = findfirst(i -> strip(String(res[i, 1])) == analysis && Int(res[i, 2]) == 1, 1:size(res, 1))
P = Float64(res[row, 4])
modelnote = analysis == "all" ? "unconstrained shell" : "rigid-section (global) shell model"
weld_length = 3.0; weld_spacing = 18.0; B = 3.0; D = 3.0 - 0.074; R = 0.199; scale_frac = 0.035

d = readdlm(path, ','; skipstart = 1)
sid = Int.(d[:, 1]); iid = Int.(d[:, 2]); jid = Int.(d[:, 3])
xyz = Float64.(d[:, 4:6]); u = Float64.(d[:, 7:9])
nn = maximum(iid); nz = maximum(jid); L = maximum(xyz[:, 3])
umag = vec(sqrt.(sum(u .^ 2; dims = 2))); sc = scale_frac * L / maximum(umag)
def = xyz .+ sc .* u
idx = Dict((sid[k], iid[k], jid[k]) => k for k in eachindex(sid))
I = Int[]; J = Int[]; K = Int[]
for jj in 1:nz-1, ss in 1:2, ii in 1:nn-1
    a = idx[(ss, ii, jj)]; b = idx[(ss, ii, jj + 1)]; c = idx[(ss, ii + 1, jj + 1)]; e = idx[(ss, ii + 1, jj)]
    push!(I, a - 1); push!(J, b - 1); push!(K, c - 1)
    push!(I, a - 1); push!(J, c - 1); push!(K, e - 1)
end
# element edges as one polyline with null breaks
ex = Any[]; ey = Any[]; ez = Any[]
function seg!(a, b)
    push!(ex, def[a, 1]); push!(ey, def[a, 2]); push!(ez, def[a, 3])
    push!(ex, def[b, 1]); push!(ey, def[b, 2]); push!(ez, def[b, 3])
    push!(ex, nothing); push!(ey, nothing); push!(ez, nothing)
end
for jj in 1:nz, ss in 1:2, ii in 1:nn-1
    seg!(idx[(ss, ii, jj)], idx[(ss, ii + 1, jj)])
end
for jj in 1:nz-1, ss in 1:2, ii in 1:nn
    seg!(idx[(ss, ii, jj)], idx[(ss, ii, jj + 1)])
end
# welds
xp = xyz[:, 1] .- (sid .- 1) .* B
nw = floor(Int, (L - weld_length) / weld_spacing) + 1
centers = range(L / 2 - (nw - 1) * weld_spacing / 2, L / 2 + (nw - 1) * weld_spacing / 2, nw)
near(x, x0, tol) = abs(x - x0) <= tol + 1e-9
corner = ((sid .== 1) .& near.(xp, B, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R))) .|
         ((sid .== 2) .& near.(xp, 0.0, R) .& (near.(xyz[:, 2], D, R) .| near.(xyz[:, 2], 0.0, R)))
inweld = [any(near(z, c, weld_length / 2) for c in centers) for z in xyz[:, 3]]
w = findall(corner .& inweld .& (xyz[:, 3] .> 0) .& (xyz[:, 3] .< L))
# loads: cones at every 4th end node, pointing into the member
ends = findall((xyz[:, 3] .≈ 0.0) .| (xyz[:, 3] .≈ L))[1:4:end]
cx = [xyz[k, 1] for k in ends]; cy = [xyz[k, 2] for k in ends]
cz = [xyz[k, 3] for k in ends]                                   # arrow tips on the end nodes
cw = [xyz[k, 3] > L / 2 ? -1.0 : 1.0 for k in ends]
shaft = 2.0
ax_ = Any[]; ay_ = Any[]; az_ = Any[]
for (k, n) in enumerate(ends)
    push!(ax_, cx[k]); push!(ay_, cy[k]); push!(az_, cz[k]); push!(ax_, cx[k]); push!(ay_, cy[k]); push!(az_, cz[k] - cw[k] * shaft)
    push!(ax_, nothing); push!(ay_, nothing); push!(az_, nothing)
end
# mid-length section
jm = (nz + 1) ÷ 2
sec_traces = String[]
for ss in 1:2
    rows = [idx[(ss, ii, jm)] for ii in 1:nn]
    push!(sec_traces, "{type:'scatter',x:[$(join(xyz[rows, 1], ','))],y:[$(join(xyz[rows, 2], ','))],mode:'lines',line:{color:'#8a8a8a',width:4},showlegend:false,hoverinfo:'skip',xaxis:'x',yaxis:'y'}")
    push!(sec_traces, "{type:'scatter',x:[$(join(def[rows, 1], ','))],y:[$(join(def[rows, 2], ','))],mode:'lines',line:{color:'#2a78d6',width:4},showlegend:false,hoverinfo:'skip',xaxis:'x',yaxis:'y'}")
end
rx = maximum(def[:, 1]) - minimum(def[:, 1]); ry = maximum(def[:, 2]) - minimum(def[:, 2]); rz = L + 2 * shaft
ar = round.([rx, ry, rz] ./ max(rx, ry, rz) .* 1.9; digits = 3)
js(v) = join((x === nothing ? "null" : string(round(x, digits = 5)) for x in v), ',')

html = """
<!doctype html>
<html><head><meta charset="utf-8"><title>Built-up upright global FT buckling mode</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>html,body{margin:0;height:100%;font-family:Helvetica,Arial,sans-serif;background:#fff}#plot{width:min(100vw,1400px);height:min(90vh,760px);margin:0 auto;overflow:hidden}</style></head>
<body><div id="plot"></div>
<script>
const data = [
 {type:'mesh3d',x:[$(js(def[:, 1]))],y:[$(js(def[:, 2]))],z:[$(js(def[:, 3]))],i:[$(join(I, ','))],j:[$(join(J, ','))],k:[$(join(K, ','))],
  intensity:[$(js(umag ./ maximum(umag)))],colorscale:'Viridis',cmin:0,cmax:1,flatshading:true,
  lighting:{ambient:0.9,diffuse:0.2,specular:0.0},colorbar:{title:{text:'normalized<br>displacement'},len:0.5,x:0.5},
  hoverinfo:'skip',name:'deformed shell',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ex))],y:[$(js(ey))],z:[$(js(ez))],line:{color:'rgba(0,0,0,0.35)',width:1},hoverinfo:'skip',showlegend:false,scene:'scene'},
 {type:'scatter3d',mode:'markers',x:[$(js(def[w, 1]))],y:[$(js(def[w, 2]))],z:[$(js(def[w, 3]))],marker:{color:'#e34948',size:3.5},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'scatter3d',mode:'lines',x:[$(js(ax_))],y:[$(js(ay_))],z:[$(js(az_))],line:{color:'#eb6834',width:4},showlegend:false,hoverinfo:'skip',scene:'scene'},
 {type:'cone',x:[$(js(cx))],y:[$(js(cy))],z:[$(js(cz))],u:[$(js(zeros(length(cx))))],v:[$(js(zeros(length(cx))))],w:[$(js(cw))],
  sizemode:'absolute',sizeref:0.55,anchor:'tip',colorscale:[[0,'#eb6834'],[1,'#eb6834']],showscale:false,showlegend:false,hoverinfo:'skip',scene:'scene'},
 $(join(sec_traces, ",\n "))
];
const layout = {
 title:{text:'Global flexural-torsional buckling mode, two-C welded upright, L = $(Int(L)) in, pinned warping-free, 3 in welds at 18 in:  P<sub>cre</sub> = $(round(P, digits = 1)) kips ($(modelnote))',x:0.02,xanchor:'left',font:{size:15}},
 scene:{domain:{x:[0,0.58],y:[0,1]},aspectmode:'manual',aspectratio:{x:$(ar[1]),y:$(ar[2]),z:$(ar[3])},xaxis:{visible:false},yaxis:{visible:false},zaxis:{visible:false},
        camera:{projection:{type:'orthographic'},eye:{x:-2.0,y:-2.4,z:0.75},center:{x:0,y:0,z:0},up:{x:0,y:0,z:1}},dragmode:'orbit'},
 xaxis:{domain:[0.62,0.98],title:{text:'X (in)'},scaleanchor:'y',scaleratio:1,zeroline:false},
 yaxis:{domain:[0.2,0.8],title:{text:'Y (in)'},zeroline:false},
 showlegend:false,margin:{l:30,r:20,t:60,b:30},autosize:true,paper_bgcolor:'#fff'};
Plotly.newPlot('plot', data, layout, {responsive:true, displaylogo:false});
</script></body></html>
"""
base = splitext(basename(path))[1]
tmp = joinpath(mktempdir(), base * "_plotly.html")
write(tmp, html)
out = joinpath(@__DIR__, base * "_plotly.html")
cp(tmp, out; force = true)
println("wrote ", out, "  (", round(filesize(out) / 1e6, digits = 1), " MB)")
