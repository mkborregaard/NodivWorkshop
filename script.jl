# Node-based analysis of bird diversity, in environmental (birds_e) and geographic
# (birds_g) space. Run preprocess.jl first to build the cleaned inputs in
# data/clean/; this script loads them, builds the assemblages, computes and caches
# the node analysis, and explores the results.
#
# Plotting is Makie: NodivMakie for trees, maps and node panels, GLMakie for interactive
# windows. Run it in the REPL or VS Code; each figure is kept in a variable, so evaluate
# the variable to show it again.

using CSV, DataFrames, JLD2, Random
using Statistics, LogExpFunctions, GLM, MultivariateStats, StatsFuns
using Clustering, Graphs
using SpatialEcology, Phylo, Nodiv
import CairoMakie            # only for saving vector (PDF) files; GLMakie saves raster formats
using GLMakie, NodivMakie
GLMakie.activate!()          # loading a backend activates it, so make sure GLMakie is the one

include("functions.jl")

set_theme!(colormap = Reverse(:Spectral))

### Load the cleaned inputs (from preprocess.jl) and build the assemblages -----

tree = sort!(parsenewick(read("data/clean/tree.nwk", String)))
phylocom_e  = strsite!(CSV.read("data/clean/phylocom_e.csv", DataFrame))
coords_e    = strsite!(CSV.read("data/clean/coords_e.csv", DataFrame))
sitestats_e = CSV.read("data/clean/sitestats_e.csv", DataFrame)
phylocom_g  = strsite!(CSV.read("data/clean/phylocom_g.csv", DataFrame))
coords_g    = strsite!(CSV.read("data/clean/coords_g.csv", DataFrame))
sitestats_g = CSV.read("data/clean/sitestats_g.csv", DataFrame)
sitestats_g.ID_geo = string.(sitestats_g.ID_geo)
traits      = CSV.read("data/clean/traits.csv", DataFrame)

# coordinates were pre-aligned to each phylocom's site order in preprocessing, so
# they slot straight into the Assemblage (SpatialEcology aligns coords by row order).
birds_e = Assemblage(phylocom_e, coords_e)
addsitestats!(birds_e, sitestats_e, :ID_env)   # PC bins, area, occupancy, ...
addtraits!(birds_e, traits, :species)
richness_e = mapfigure(birds_e; title = "Environmental: species richness", label = "species")
# save("figures/Env species richness.png", richness_e)

birds_g = Assemblage(phylocom_g, coords_g)
addsitestats!(birds_g, sitestats_g, :ID_geo)   # CHELSA bioclim, PC1-3, area, ...
addtraits!(birds_g, traits, :species)
richness_g = mapfigure(birds_g; title = "Geographic: species richness", label = "species",
                       figure = (; size = (1000, 500)))
# save("figures/Geo species richness.png", richness_g)

### Heavy step: divergence metrics + SOS for every node, both spaces, cached to disk ----------
# `node_metrics` computes the divergence metrics (GND, RMS-SOS, ...) together with the
# per-cell SOS they are built on. This is the slow part - randomisations over the whole
# tree, and ~18k cells for the geographic scan - so cache it: re-running the script just
# reloads the results and jumps straight to the plotting below. The seeded `rng` makes a
# regenerated cache identical to this one, whatever the number of threads (`julia -t auto`).
cachefile = "data/node_analysis.jld2"
if !isfile(cachefile)
    res_e = node_metrics(birds_e, tree; nsims = 200, rng = Xoshiro(1))
    res_g = node_metrics(birds_g, tree; nsims = 200, rng = Xoshiro(2))
    jldsave(cachefile; res_e, res_g)
end
# each a NodeMetrics: nodes (in tree order), gnd, rms, sd, ses, pval, varying and sos
res_e, res_g = load(cachefile, "res_e", "res_g")

metric = :rms # alternatives are :pval and :gnd
threshold = 2
metric_e = getfield(res_e, metric)
metric_g = getfield(res_g, metric)

# Minimum number of shared occupied cells for the correlation behind `sos_distances`
# (SOS-pattern similarity, used by the explorers' ordination and the grouping section).
# Below it `sos_distances` pins the pair at distance 1 rather than trusting a correlation
# fit on a handful of cells - which is also what keeps disjoint pairs at the maximum. 3 is
# the smallest overlap where |r| is not trivially 1.
MINOVERLAP = 3
SIMCUT = 0.7

### ---- Exploratory plotting (from the cached NodeMetrics; `_e` vs `_g`) ----- ###

# strongly divergent nodes in each space (`metric` above `threshold`)
divergent_e = divergent_nodes(res_e; by = metric, threshold)
divergent_g = divergent_nodes(res_g; by = metric, threshold)
divergent = divergent_e ∩ divergent_g

# The metric of just the divergent nodes mapped onto the tree (markers only at `nodes`:
# a node missing from the Dict would get a transparent fill but still its outline).
# GND is a proportion, so it gets a fixed 0-1 colour range, as plot_gnd used.
function metric_tree(tree, values, nodes, title)
    colorrange = metric === :gnd ? (0, 1) : Makie.automatic
    fig, ax, p = treeplot(tree; treetype = :fan, showtips = false,
                          nodecolor = Dict(n => values[n] for n in nodes), shownodes = nodes,
                          markersize = 8, strokewidth = 0.5, colormap = :YlOrRd, colorrange, axis = (; title),
                          figure = (; size = (800, 700)))
    Colorbar(fig[1, 2], p; label = string(metric))
    fig
end
metric_tree_e = metric_tree(tree, metric_e, divergent_e, "Environmental: divergent nodes, $metric")
# save("figures/Env divergent nodes treeplot.png", metric_tree_e)
metric_tree_g = metric_tree(tree, metric_g, divergent_g, "Geographic: divergent nodes, $metric")
# save("figures/Geo divergent nodes treeplot.png", metric_tree_g)

# SOS of the most divergent node mapped onto each space (cached SOS, no recompute)
focal_e = argmax(n -> metric_e[n], divergent_e)
sosmap_e = mapfigure(res_e.sos[focal_e], birds_e; colormap = :RdYlBu, colorrange = (-8, 8),
                     title = "Environmental: SOS of $focal_e", label = "SOS")
# save("figures/Env SOS $focal_e.png", sosmap_e)
focal_g = argmax(n -> metric_g[n], divergent_g)
sosmap_g = mapfigure(res_g.sos[focal_g], birds_g; colormap = :RdYlBu, colorrange = (-8, 8),
                     title = "Geographic: SOS of $focal_g", label = "SOS", figure = (; size = (1000, 500)))
# save("figures/Geo SOS $focal_g.png", sosmap_g)

# The interactive entry point: the fan tree with the divergent nodes marked, next to the
# SOS map and the two child clades' maps, and an ordination of the divergent nodes by
# SOS-pattern similarity. Each opens on its space's most divergent node (focal_e, focal_g
# above). Click a node or a branch on the tree, or a point in the ordination, to show that
# node; hover for labels. The Birds of the World images in bow_images/ are private and not
# in the repo: they are used only if that folder is there.
imagedir = "bow_images/workshop_species"
explorer_options = (; metric, images = isdir(imagedir) ? imagedir : nothing,
                    imageoptions = (; whitebackground = true),
                    ordinationkw = (; minoverlap = MINOVERLAP))
explorer_fig_e, explorer_e = nodeexplorer(birds_e, tree, res_e; nodes = divergent_e, explorer_options...)
explorer_fig_g, explorer_g = nodeexplorer(birds_g, tree, res_g; nodes = divergent_g, explorer_options...)

# Link the two: a node picked in one space is shown in the other too, if it has an SOS there
for (from, to, res) in ((explorer_g, explorer_e, res_e), (explorer_e, explorer_g, res_g))
    on(from.panel.node) do n
        n != to.panel.node[] && hassos(tree, res.sos, n) && (to.panel.node[] = n)
    end
end

# each explorer in its own window (NODIVWORKSHOP_WINDOWS=false skips this, e.g. headless)
if isinteractive() && get(ENV, "NODIVWORKSHOP_WINDOWS", "true") != "false"
    display(GLMakie.Screen(), explorer_fig_e)
    display(GLMakie.Screen(), explorer_fig_g)
end

# ordinate the divergent nodes of both spaces by SOS-pattern similarity (cached SOS ->
# `sos_distances` from Nodiv -> classical MDS in NodivMakie's `sosordination`)
D_g = sos_distances(res_g, divergent; minoverlap = MINOVERLAP)
D_e = sos_distances(res_e, divergent; minoverlap = MINOVERLAP)
sos_mds_plot(D, nodes, title) =
    ordinationplot(sosordination(D, nodes); nodelabels = true, axis = (; title)).figure
mds_e = sos_mds_plot(D_e, divergent, "Environmental: SOS-pattern similarity")
# save("figures/Env SOS-pattern similarity.png", mds_e)
mds_g = sos_mds_plot(D_g, divergent, "Geographic: SOS-pattern similarity")
# save("figures/Geo SOS-pattern similarity.png", mds_g)

# parent/SOS/children panel for one node (4th arg = cached SOS, no recompute); also
# `explorer_e.panel.node[] = focal` shows it in the explorer
focal = "Node 17672"   # node names are numbered by data/clean/tree.nwk: re-running preprocess.jl renumbers them
panel_e, _ = nodepanel(birds_e, tree, focal, res_e)
# save("figures/Env node panel $focal.png", panel_e)
panel_g, _ = nodepanel(birds_g, tree, focal, res_g)
# save("figures/Geo node panel $focal.png", panel_g)

allnodes = collect(keys(metric_e))
dat = DataFrame(
    :log_g => [log(metric_g[n]) for n in allnodes], #NB logit or log, depends on metric
    :log_e => [log(metric_e[n]) for n in allnodes]
)
dat = filter(row -> all(x -> !ismissing(x) && isfinite(x), row), dat)

metric_scatter = scatter(dat.log_g, dat.log_e;
                         axis = (; xlabel = "log geo $metric", ylabel = "log env $metric"))
ablines!(metric_scatter.axis, 0, 1; color = :red)   # the 1:1 line
# save("figures/Env vs geo $metric.png", metric_scatter)

rms_fit = lm(@formula(log_e ~ log_g), dat)


occupied_e = Dict(node => noccupied(get_clade(birds_e, tree, node)) for node in allnodes)
occupied_hist = hist(collect(values(occupied_e)); axis = (; xlabel = "occupied env sites", ylabel = "nodes"))
# save("figures/Env occupied sites histogram.png", occupied_hist)

occupied_tree = let (fig, ax, p) = treeplot(tree; treetype = :fan, nodecolor = occupied_e, showtips = false,
                                            markersize = 5, figure = (; size = (800, 700)))
    Colorbar(fig[1, 2], p; label = "occupied env sites")
    fig
end
# save("figures/Env occupied sites treeplot.png", occupied_tree)

occupied_scatter = scatter([occupied_e[n] for n in allnodes], [metric_e[n] for n in allnodes];
                           axis = (; xlabel = "occupied env sites", ylabel = "env $metric"))
# save("figures/Env $metric vs occupied sites.png", occupied_scatter)


nspecies_g = Dict(node => nspecies(get_clade(birds_g, tree, node)) for node in allnodes)
nspecies_scatter = scatter([log(nspecies_g[n]) for n in allnodes], [metric_g[n] for n in allnodes];
                           axis = (; xlabel = "log number of species in clade", ylabel = "geo $metric"))
# save("figures/Geo $metric vs clade species.png", nspecies_scatter)


### ===========================================================================
### Grouping divergent nodes by SOS-pattern similarity
### The MDS scatter (`sos_mds_plot`) above is read together with its eigenvalue diagnostic
### (1); the primary read is the complete-linkage clustered heatmap (2), with a thresholded
### similarity graph (3) as the confirmatory secondary. Distances come from `sos_distances`
### in Nodiv: 1 - |r| over the shared occupied cells, with the minimum-overlap floor
### `MINOVERLAP`. Each space's distances are computed once and every view is derived from
### them. The two spaces are run separately and their magnitudes are NOT compared
### (environmental "occupancy" is over tens of PC bins, geographic over ~18k cells).
### ===========================================================================

# (2) PRIMARY VIEW. Complete-linkage hierarchical clustering of the SOS-pattern distances `D`
# of `nodes`. Genuine groups (if any) appear as dense blocks on the diagonal of the heatmap
# below; the orthogonal/disjoint background stays uniform. Complete linkage is the
# conservative default - it groups only all-pairs-similar nodes and will not chain marginal
# pairs. Cut the tree at |r| >= `simcut` (i.e. distance height 1 - simcut).
# Returns (hclust, groups::Dict node=>cluster).
function sos_clusters(D, nodes; simcut, linkage = :complete)
    hc  = hclust(D; linkage)
    grp = cutree(hc; h = 1 - simcut)
    hc, Dict(node => grp[i] for (i, node) in enumerate(nodes))
end

# The clustering from `sos_clusters` drawn as a dendrogram-ordered |r| heatmap.
# Layout is the standard clustermap: a horizontal dendrogram on the LEFT and the heatmap on
# the RIGHT, sharing the y-axis (linked), so the leaves stay aligned with the heatmap rows;
# the colourbar only takes width from the heatmap. The dendrogram is Makie's (experimental)
# recipe, turned to put the leaves (height 0) against the heatmap; its heights are drawn at
# negative x, so the ticks show them as positive. The heatmap's y labels print in the gap
# between the two panels, so every dendrogram tip can be read straight off as a node name. The
# x-axis carries the same names (rotated). The cut clusters with more than one node are
# outlined on the diagonal with their number and colour from `plot_cluster_tree`. Tune
# `labelsize`/`figsize` for iterative exploration. The leaf order, bottom-to-top, shown on the
# axes is `nodes[hc.order]`.
function sos_cluster_heatmap(D, nodes, hc, groups, title; labelsize = 9, figsize = (1300, 1150))
    ord  = hc.order
    labs = nodes[ord]                             # node names in dendrogram-leaf order
    S    = (1 .- D)[ord, ord]                     # similarity |r|, reordered to match
    n    = length(labs)
    grp  = [groups[node] for node in labs]        # cluster of each leaf, in leaf order

    fig  = Figure(; size = figsize)
    Label(fig[0, 1:3], title; fontsize = 18, font = :bold)
    dend = Axis(fig[1, 1]; xlabel = "1 - |r|", xtickformat = xs -> string.(round.(abs.(xs); digits = 2)),
                xgridvisible = false,
                ygridvisible = false, yticksvisible = false, yticklabelsvisible = false,
                leftspinevisible = false, topspinevisible = false, rightspinevisible = false)
    hm   = Axis(fig[1, 2]; xticks = (1:n, labs), yticks = (1:n, labs),
                xticklabelrotation = pi / 2, xticklabelsize = labelsize,
                yticklabelsize = labelsize)
    # leaf hc.order[k] at y = k, lined up with heatmap row k; each merge at its height
    dendrogram!(dend, Makie.hcl_nodes(hc; useheight = true); absolute = true, rotation = :right,
                color = :black)
    h = heatmap!(hm, 1:n, 1:n, S; colormap = :viridis, colorrange = (0, 1))
    Colorbar(fig[1, 3], h; label = "|r|")

    idmap = cluster_idmap(groups)
    colors = clustercolors(length(idmap))
    for (c, i) in idmap                           # cutree clusters are contiguous in `ord`
        rows = findall(==(c), grp)
        lo, hi = extrema(rows)
        box = Rect2d(lo - 0.5, lo - 0.5, hi - lo + 1, hi - lo + 1)
        poly!(hm, box; color = :transparent, strokecolor = :black, strokewidth = 4)
        poly!(hm, box; color = :transparent, strokecolor = colors[i], strokewidth = 2)
        textlabel!(hm, Point2d(hi + 0.5, hi + 0.5); text = string(i), fontsize = 11,
                   background_color = colors[i], strokewidth = 1, padding = 2,
                   cornerradius = 2,   # inside the corner at the top edge, not clipped
                   text_align = hi == n ? (:right, :top) : (:left, :bottom))
    end

    linkyaxes!(dend, hm)
    xlims!(dend, -1.02maximum(hc.heights), 0)
    limits!(hm, 0.5, n + 0.5, 0.5, n + 0.5)
    colsize!(fig.layout, 1, Relative(0.20))
    fig
end

# SECONDARY / CONFIRMATORY. Thresholded similarity graph + modularity communities, from the
# SOS-pattern distance matrix `D` of `nodes`. Build an edge only between pairs with
# |r| >= `simthresh`; the overlap floor is already enforced inside `sos_distances`
# (thin-overlap and disjoint pairs sit at distance 1, |r| = 0, so they never become edges).
# Co-patterned nodes then fall out as the communities that maximise the graph's modularity Q
# (greedy agglomeration, Clauset-Newman-Moore): groups more densely linked inside than a graph
# with the same degrees would be by chance - unlike connected components, which chain any
# path of links into one group. Returns the non-singleton communities as vectors of node
# names, and Q.
function sos_similarity_communities(D, nodes; simthresh)
    n = length(nodes)
    g = Graphs.SimpleGraph(n)
    for i in 1:n, j in i+1:n
        (1 - D[i, j]) >= simthresh && Graphs.add_edge!(g, i, j)
    end
    m = Graphs.ne(g)
    comm = collect(1:n)
    m == 0 && return (; communities = Vector{String}[], modularity = 0.0)
    deg = Graphs.degree(g)
    while true
        best, merge = 0.0, nothing
        for a in unique(comm), b in unique(comm)
            a < b || continue
            links = count(e -> minmax(comm[Graphs.src(e)], comm[Graphs.dst(e)]) == (a, b), Graphs.edges(g))
            links == 0 && continue
            dQ = links / m - sum(deg[comm .== a]) * sum(deg[comm .== b]) / (2m^2)
            dQ > best && ((best, merge) = (dQ, (a, b)))
        end
        merge === nothing && break
        comm[comm .== merge[2]] .= merge[1]
    end
    groups = filter(c -> length(c) > 1, [findall(==(c), comm) for c in unique(comm)])
    (; communities = [[nodes[i] for i in c] for c in groups], modularity = Graphs.modularity(g, comm))
end

function print_clusters(space, groups, communities)
    sizes = sos_cluster_sizes(groups)
    idmap = cluster_idmap(groups)
    println("$space: $(sizes.nclusters) clusters, $(sizes.nsingletons) singletons")
    for (c, i) in sort(collect(idmap); by = last)
        println("  cluster $i: ", join(sizes.members[c], ", "))
    end
    println("  modularity Q = ", round(communities.modularity; digits = 3))
    for (i, members) in enumerate(communities.communities)
        println("  community $i: ", join(members, ", "))
    end
end

# --- Run both spaces on the divergent set (no cross-space magnitude comparison) ----------
# (1) ONE-TIME DIAGNOSTIC, not the analysis. Fit MDS at a higher dimension and look at the
# eigenvalue spectrum: if axes 3+ carry weight comparable to axes 1-2, the 2-D scatter is a
# projection artefact and the "ring" is the honest report of near-equidistance.
mds_eig_g = eigenvalueplot(sosordination(D_g, divergent; maxoutdim = 10);
                           axis = (; title = "Geographic: MDS eigenvalues (n = $(length(divergent)))")).figure
# save("figures/Geo MDS eigenvalues.png", mds_eig_g)
mds_eig_e = eigenvalueplot(sosordination(D_e, divergent; maxoutdim = 10);
                           axis = (; title = "Environmental: MDS eigenvalues (n = $(length(divergent)))")).figure
# save("figures/Env MDS eigenvalues.png", mds_eig_e)

hc_g, groups_g = sos_clusters(D_g, divergent; simcut = SIMCUT)
hc_e, groups_e = sos_clusters(D_e, divergent; simcut = SIMCUT)

heat_g = sos_cluster_heatmap(D_g, divergent, hc_g, groups_g, "Geographic: SOS clusters")
# save("figures/Geo SOS clusters.png", heat_g)
heat_e = sos_cluster_heatmap(D_e, divergent, hc_e, groups_e, "Environmental: SOS clusters")
# save("figures/Env SOS clusters.png", heat_e)

communities_g = sos_similarity_communities(D_g, divergent; simthresh = SIMCUT)
communities_e = sos_similarity_communities(D_e, divergent; simthresh = SIMCUT)

# Clusters mapped back onto the phylogeny (numbers and colours match the heatmap outlines)
tree_clusters_g = plot_cluster_tree(tree, groups_g, "Geographic: SOS clusters on the phylogeny")
# save("figures/Geo SOS clusters on the phylogeny.png", tree_clusters_g)
tree_clusters_e = plot_cluster_tree(tree, groups_e, "Environmental: SOS clusters on the phylogeny")
# save("figures/Env SOS clusters on the phylogeny.png", tree_clusters_e)

# Report the grouping result directly (this IS the scientific output):
# mostly singletons with a few multi-node clusters = "largely idiosyncratic, with named
# co-patterned exceptions"; substantial blocks = real groups to map onto the phylogeny.
print_clusters("Geographic", groups_g, communities_g)
print_clusters("Environmental", groups_e, communities_e)


### --- Two node-level views of the divergent set ------------------------------------------

# Fan tree showing ONLY the divergent nodes, each labelled with its name in a small pale box
# so it stays readable over the branches; the rest of the tree is a plain grey skeleton. The
# boxes will overlap if packed too tightly, so widen the figure size (or drop
# `nodelabelsize`) until they clear; the redundant "Node " is dropped so the boxes stay small.
divergent_tree = treeplot(tree; treetype = :fan, showtips = false, branchcolor = :gray75,
                          nodelabels = Dict(n => replace(n, "Node " => "") for n in divergent),
                          nodelabelbackground = (:lightyellow, 0.85), nodelabelsize = 11,
                          nodelabelalign = (:center, :center), nodelabeloffset = (0, 0),
                          figure = (; size = (1600, 1600))).figure
# save("figures/Divergent nodes labelled.png", divergent_tree)
# plot_node_pdf(birds_g, tree, divergent, res_g, "figures/divergent_node_panels_geo.pdf")
# plot_node_pdf(birds_e, tree, divergent, res_e, "figures/divergent_node_panels_env.pdf")


### --- Traits: PCA of the AVONET morphometrics ---------------------------------------------

# Shapiro-Francia W': the squared correlation of the sorted values with the normal quantiles
function normality(x)
    n = length(x)
    q = norminvcdf.(((1:n) .- 0.375) ./ (n + 0.25))
    cor(sort(x), q)^2
end

# PCA of the columns of `df`: a column is logged where that raises its W' by more than
# `tol`, then all are z-transformed. Returns the PCA, the logged columns and the
# z-transformed matrix (species x traits) it was fit on.
function traitpca(df; tol = 0.01)
    X = Matrix{Float64}(df)
    logged = [all(>(0), x) && normality(log.(x)) - normality(x) > tol for x in eachcol(X)]
    X[:, logged] .= log.(X[:, logged])
    Z = (X .- mean(X; dims = 1)) ./ std(X; dims = 1)
    (; pca = fit(PCA, permutedims(Z); pratio = 1), logged = names(df)[logged], z = Z)
end

speciestraits = SpatialEcology.traits(birds_g)
trait_pca = traitpca(speciestraits[:, 11:21])
pca_explained = principalvars(trait_pca.pca) ./ var(trait_pca.pca)

pcs = DataFrame(permutedims(predict(trait_pca.pca, permutedims(trait_pca.z))[1:4, :]), ["pca$i" for i in 1:4])
pcs.species = speciestraits.name
addtraits!(birds_e, pcs, :species)
addtraits!(birds_g, pcs, :species)

cross2(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])

# Convex hull of 2-d points (Andrew's monotone chain), counter-clockwise, first point not repeated
function convexhull(pts)
    ps = sort(unique(pts); by = p -> (p[1], p[2]))
    length(ps) < 3 && return ps
    function half(ps)
        h = eltype(ps)[]
        for p in ps
            while length(h) >= 2 && cross2(h[end-1], h[end], p) <= 0
                pop!(h)
            end
            push!(h, p)
        end
        h
    end
    lower, upper = half(ps), half(reverse(ps))
    [lower[1:end-1]; upper[1:end-1]]
end

function polyarea(h)
    n = length(h)
    n < 3 && return 0.0
    abs(sum(h[i][1] * h[mod1(i + 1, n)][2] - h[mod1(i + 1, n)][1] * h[i][2] for i in 1:n)) / 2
end

# The intersection of two convex polygons (Sutherland-Hodgman: clip `a` by each edge of `b`)
function clippolygon(a, b)
    out = a
    for i in eachindex(b)
        isempty(out) && break
        p, q = b[i], b[mod1(i + 1, length(b))]
        inside(x) = cross2(p, q, x) >= 0
        cut(s, e) = s + (e - s) * (cross2(p, q, s) / (cross2(p, q, s) - cross2(p, q, e)))
        input, out = out, eltype(a)[]
        for j in eachindex(input)
            cur, prev = input[j], input[mod1(j - 1, length(input))]
            if inside(cur)
                inside(prev) || push!(out, cut(prev, cur))
                push!(out, cur)
            elseif inside(prev)
                push!(out, cut(prev, cur))
            end
        end
    end
    out
end

# The overlap of two convex hulls as a proportion of the smaller one; NaN if either has no area
function hulloverlap(h1, h2)
    a = min(polyarea(h1), polyarea(h2))
    a > 0 || return NaN
    polyarea(clippolygon(h1, h2)) / a
end

# Species => point in trait space, from the columns `x` and `y` of an assemblage's traits
function traitpoints(asm, x, y)
    t = SpatialEcology.traits(asm)
    Dict(zip(t.name, Point2d.(t[!, x], t[!, y])))
end

# The trait-space points of the species of each of `node`'s two child clades
childpoints(tree, node, pts) =
    [[pts[sp] for sp in nodespecies(tree, getnodename(tree, c)) if haskey(pts, sp)]
     for c in getchildren(tree, node)[1:2]]

closedhull(pts) = (h = convexhull(pts); length(h) < 3 ? Point2d[] : [h; h[1:1]])

# All species in trait space in grey, with the two child clades of `node` (an Observable)
# in the explorer's clade colours, the smaller clade on top, each outlined by its convex hull
function traitpanel!(gp, asm, tree, node, x, y; axis = (;))
    pts = traitpoints(asm, x, y)
    colors = cladecolors(:RdYlBu)
    ax = Axis(gp; xgridvisible = false, ygridvisible = false, axis...)
    scatter!(ax, collect(values(pts)); color = :gray80, markersize = 3, inspectable = false)
    clades = lift(n -> childpoints(tree, n, pts), node)
    for (k, color) in enumerate(colors)
        cladepts = lift(c -> c[k], clades)
        sc = scatter!(ax, cladepts; color, markersize = 5, inspectable = false)
        on(c -> translate!(sc, 0, 0, length(c[k]) <= length(c[3 - k])), clades; update = true)
        hull = lines!(ax, lift(closedhull, cladepts); color, linewidth = 2, inspectable = false)
        translate!(hull, 0, 0, 2)
    end
    ax
end

# A node explorer of the two spaces and trait space: the tree, the SOS of the node shown in
# geographic and environmental space, and its two child clades on PCA axes 1-2 and 3-4
function traitexplorer(tree, marked, values, birds_g, res_g, birds_e, res_e, explained;
                       images = nothing, imageoptions = (;))
    fig = Figure(; size = (1600, 850))
    node = Observable(argmax(n -> marked[n], keys(marked)))
    tr = explorertree!(fig[1, 1], tree, node, marked; values, label = "geo $metric",
                       selectable = n -> hassos(tree, res_g.sos, n) && hassos(tree, res_e.sos, n),
                       unselectable = "no SOS in both spaces", images, imageoptions,
                       rangesize = birds_g)
    panels = fig[1, 2] = GridLayout()
    sosmap!(panels[1, 1], birds_g, node, res_g; title = "Geographic SOS")
    sosmap!(panels[1, 2], birds_e, node, res_e; title = "Environmental SOS")
    pclabel(i) = "pca$i ($(round(100explained[i]; digits = 1))%)"
    for (col, (i, j)) in enumerate(((1, 2), (3, 4)))
        traitpanel!(panels[2, col], birds_g, tree, node, Symbol("pca$i"), Symbol("pca$j");
                    axis = (; xlabel = pclabel(i), ylabel = pclabel(j)))
    end
    colsize!(fig.layout, 1, Relative(0.45))
    DataInspector(fig)
    fig, tr
end

trait_marked = Dict(n => metric_g[n] for n in divergent_e ∪ divergent_g)
trait_fig, trait_explorer = traitexplorer(tree, trait_marked, metric_g, birds_g, res_g, birds_e, res_e,
                                          pca_explained; explorer_options.images, explorer_options.imageoptions)
if isinteractive() && get(ENV, "NODIVWORKSHOP_WINDOWS", "true") != "false"
    display(GLMakie.Screen(), trait_fig)
end

pts12 = traitpoints(birds_g, :pca1, :pca2)
trait_overlap = Dict(n => hulloverlap(convexhull.(childpoints(tree, n, pts12))...) for n in allnodes)

overlap_dat = DataFrame(
    :overlap => [trait_overlap[n] for n in allnodes],
    :log_g => [log(metric_g[n]) for n in allnodes],
    :log_e => [log(metric_e[n]) for n in allnodes]
)
overlap_dat = filter(row -> all(isfinite, row), overlap_dat)
overlap_fit_g = lm(@formula(log_g ~ overlap), overlap_dat)
overlap_fit_e = lm(@formula(log_e ~ overlap), overlap_dat)

overlap_scatter = let fig = Figure(; size = (1100, 500))
    for (col, (y, lmfit, title)) in enumerate(((:log_g, overlap_fit_g, "Geographic: trait overlap and $metric"),
                                             (:log_e, overlap_fit_e, "Environmental: trait overlap and $metric")))
        ax = Axis(fig[1, col]; title, xlabel = "trait overlap (pca1-2)", ylabel = "log $metric")
        scatter!(ax, overlap_dat.overlap, overlap_dat[!, y]; markersize = 5)
        ablines!(ax, coef(lmfit)...; color = :red)
    end
    fig
end
# save("figures/Trait overlap vs $metric.png", overlap_scatter)
