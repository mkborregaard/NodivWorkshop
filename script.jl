# Node-based analysis of bird diversity, in environmental (birds_e) and geographic
# (birds_g) space. Run preprocess.jl first to build the cleaned inputs in
# data/clean/; this script loads them, builds the assemblages, computes and caches
# the node analysis, and explores the results.
#
# Plotting is Makie: NodivMakie for trees, maps and node panels, GLMakie for interactive
# windows. Run it in the REPL or VS Code; each figure is kept in a variable, so evaluate
# the variable to show it again.

using CSV, DataFrames, JLD2
using MultivariateStats, Statistics, LogExpFunctions, GLM
using SpatialEcology, Phylo, Nodiv
import CairoMakie            # only for saving vector (PDF) files; GLMakie saves raster formats
using GLMakie, NodivMakie
GLMakie.activate!()          # loading a backend activates it, so make sure GLMakie is the one

set_theme!(colormap = Reverse(:Spectral))

# A map of one value per site, with a colour bar beside it
function mapfigure(args...; title = "", label = "", kw...)
    fig, ax, p = sitemap(args...; axis = (; title), kw...)
    Colorbar(fig[1, 2], p; label)
    fig
end

### Load the cleaned inputs (from preprocess.jl) and build the assemblages -----

tree = parsenewick(read("data/clean/tree.nwk", String))
strsite!(df) = (df.site = string.(df.site); df)   # site ids stay strings after CSV
phylocom_e  = strsite!(CSV.read("data/clean/phylocom_e.csv", DataFrame))
coords_e    = strsite!(CSV.read("data/clean/coords_e.csv", DataFrame))
sitestats_e = CSV.read("data/clean/sitestats_e.csv", DataFrame)
phylocom_g  = strsite!(CSV.read("data/clean/phylocom_g.csv", DataFrame))
coords_g    = strsite!(CSV.read("data/clean/coords_g.csv", DataFrame))
sitestats_g = CSV.read("data/clean/sitestats_g.csv", DataFrame)
sitestats_g.ID_geo = string.(sitestats_g.ID_geo)

# coordinates were pre-aligned to each phylocom's site order in preprocessing, so
# they slot straight into the Assemblage (SpatialEcology aligns coords by row order).
birds_e = Assemblage(phylocom_e, coords_e)
addsitestats!(birds_e, sitestats_e, :ID_env)   # PC bins, area, occupancy, ...
richness_e = mapfigure(birds_e; title = "Species richness (environmental)", label = "species")

birds_g = Assemblage(phylocom_g, coords_g)
addsitestats!(birds_g, sitestats_g, :ID_geo)   # CHELSA bioclim, PC1-3, area, ...
richness_g = mapfigure(birds_g; title = "Species richness (geographic)", label = "species",
                       figure = (; size = (1000, 500)))

### Heavy step: GND + SOS for every node, both spaces, cached to disk ----------
# `node_analysis` computes GND and the per-cell SOS together (SOS is needed for
# GND anyway). This is the slow part - randomisations over the whole tree, and
# ~18k cells for the geographic scan - so cache it: re-running the script just
# reloads the results and jumps straight to the plotting below.
cachefile = "data/node_analysis.jld2"
if !isfile(cachefile)
    res_e = node_metrics(birds_e, tree; nsims = 200)
    res_g = node_metrics(birds_g, tree; nsims = 200)
    jldsave(cachefile; res_e, res_g)
end
res_e, res_g = load(cachefile, "res_e", "res_g")   # each a NodeMetrics (gnd/rms/spatial/ses/pval + sos)

metric = :rms # alternatives are :pval, :gnd and :ses
e_metric = getfield(res_e, metric)
g_metric = getfield(res_g, metric)

### ---- Exploratory plotting (from the cached NodeMetrics; `_e` vs `_g`) ----- ###

# strongly divergent nodes in each space. by = :gnd keeps the original GND > 0.8
# selection; switch to the default (RMS-SOS > 1.5) or by = :pval to use the
# effect-size / null-calibrated scores instead.
divergent_e = divergent_nodes(res_e; by = metric, threshold = 2)
divergent_g = divergent_nodes(res_g; by = metric, threshold = 2)
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
metric_tree_e = metric_tree(tree, e_metric, divergent_e, "Divergent nodes, $metric (environmental)")
metric_tree_g = metric_tree(tree, g_metric, divergent_g, "Divergent nodes, $metric (geographic)")

# SOS of the most divergent node mapped onto each space (cached SOS, no recompute)
focal_e = argmax(n -> e_metric[n], divergent_e)
sosmap_e = mapfigure(res_e.sos[focal_e], birds_e; colormap = :RdYlBu, colorrange = (-8, 8),
                     title = "env SOS - $focal_e", label = "SOS")
focal_g = argmax(n -> g_metric[n], divergent_g)
sosmap_g = mapfigure(res_g.sos[focal_g], birds_g; colormap = :RdYlBu, colorrange = (-8, 8),
                     title = "geo SOS - $focal_g", label = "SOS", figure = (; size = (1000, 500)))

# The interactive entry point: the fan tree with the divergent nodes marked, next to the
# parent/SOS/children panel. Each opens on its space's most divergent node (focal_e,
# focal_g above). Click a node or a branch to show that node; hover for labels. The
# Birds of the World images in bow_images/ are private and not in the repo: they are
# used only if that folder is there.
imagedir = "bow_images/workshop_species"
explorer_options = (; metric, images = isdir(imagedir) ? imagedir : nothing,
                    imageoptions = (; whitebackground = true))
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

# ordinate the divergent nodes by SOS-pattern similarity (cached SOS -> distances
# from Nodiv -> MDS; presentation stays here)
function sos_mds_plot(res, nodes, title)
    coords = predict(fit(MDS, sos_distances(res, nodes); distances = true, maxoutdim = 2))
    fig, ax, _ = scatter(coords[1, :], coords[2, :];
                         axis = (; xlabel = "MDS axis 1", ylabel = "MDS axis 2", title))
    text!(ax, coords[1, :], coords[2, :]; text = nodes, fontsize = 8,
          align = (:center, :bottom), offset = (0, 4))
    fig
end
mds_e = sos_mds_plot(res_e, divergent, "SOS-pattern similarity (environmental)")
mds_g = sos_mds_plot(res_g, divergent, "SOS-pattern similarity (geographic)")

# parent/SOS/children panel for one node (4th arg = cached SOS, no recompute); also
# `explorer_e.panel.node[] = focal` shows it in the explorer
focal = "Node 17672"
panel_e, _ = nodepanel(birds_e, tree, focal, res_e)
panel_g, _ = nodepanel(birds_g, tree, focal, res_g)

same = divergent_e ∩ divergent_g

nodes = collect(keys(e_metric))
dat = DataFrame(
    :logit_g => [log(g_metric[n]) for n in nodes], #NB logit or log, depends on metric
    :logit_e => [log(e_metric[n]) for n in nodes]
)
dat = filter(row -> all(x -> !ismissing(x) && isfinite(x), row), dat)

metric_scatter = scatter(dat.logit_g, dat.logit_e;
                         axis = (; xlabel = "geo $metric", ylabel = "env $metric"))
ablines!(metric_scatter.axis, 0, 1; color = :red)   # the 1:1 line

mod = lm(@formula(logit_e ~ logit_g), dat)


sizes = Dict(node => noccupied(get_clade(birds_e, tree, node)) for node in nodes)
sizes_hist = hist(collect(values(sizes)); axis = (; xlabel = "occupied env sites", ylabel = "nodes"))

sizes_tree = let (fig, ax, p) = treeplot(tree; treetype = :fan, nodecolor = sizes, showtips = false,
                                         markersize = 5, figure = (; size = (800, 700)))
    Colorbar(fig[1, 2], p; label = "occupied env sites")
    fig
end

sizes_scatter = scatter([sizes[n] for n in nodes], [e_metric[n] for n in nodes];
                        axis = (; xlabel = "occupied env sites", ylabel = "env $metric"))


nodesizes = Dict(node => nspecies(get_clade(birds_g, tree, node)) for node in nodes)
nodesizes_scatter = scatter([log(nodesizes[n]) for n in nodes], [g_metric[n] for n in nodes];
                            axis = (; xlabel = "number of species in clade", ylabel = "geo $metric"))


### ===========================================================================
### Grouping divergent nodes by SOS-pattern similarity
### (docs/sos_pattern_grouping_design.md). Appended below the existing exploration;
### nothing above is changed. The MDS scatter (`sos_mds_plot`) is kept ONLY as the
### eigenvalue diagnostic in (1); the primary read is the complete-linkage clustered
### heatmap (2), with a thresholded similarity graph (3) as the confirmatory secondary.
### Distances come from the corrected `sos_distances` in Nodiv: 1 - |r| over the shared
### occupied cells, with a per-space minimum-overlap floor. The two spaces are run
### separately and their magnitudes are NOT compared (environmental "occupancy" is over
### tens of PC bins, geographic over ~18k cells).
### ===========================================================================

using Clustering, Graphs, LinearAlgebra

# Per-space minimum-overlap floors for the correlation's own sample size. The geographic
# scan has ~18k cells, so a floor of ~8 shared occupied cells is defensible; the
# environmental scan has only tens of PC bins, so its floor must be much smaller. Below the
# floor `sos_distances` pins the pair at distance 1 rather than trusting a correlation fit
# on a handful of cells - which is also what keeps disjoint pairs at the maximum.
const MINOVERLAP_G = 8
const MINOVERLAP_E = 3

# (1) ONE-TIME DIAGNOSTIC, not the analysis. Fit MDS at a higher dimension and look at the
# eigenvalue spectrum: if axes 3+ carry weight comparable to axes 1-2, the 2-D scatter above
# is a projection artefact and the "ring" is the honest report of near-equidistance.
function sos_mds_eigenvalues(res, nodes, title; minoverlap, method = :pearson)
    D = sos_distances(res, nodes; minoverlap, method)
    M = fit(MDS, D; distances = true, maxoutdim = min(10, length(nodes) - 1))
    λ = eigvals(M)
    fig, ax, _ = barplot(1:length(λ), λ;
                         axis = (; xlabel = "MDS axis", ylabel = "eigenvalue",
                                 title = "$title  (n = $(length(nodes)))"))
    fig
end

# The clusters worth showing: the cut clusters with more than one node, relabelled 1..m in
# the order of their `cutree` ids. The heatmap and the tree both number and colour the
# clusters by this, so the two views cross-reference directly.
function cluster_idmap(groups)
    ids = sort(collect(keys(sos_cluster_sizes(groups).members)))
    Dict(c => i for (i, c) in enumerate(ids))
end
clustercolors(m) = (c = Makie.to_colormap(:tab20); [c[mod1(i, length(c))] for i in 1:m])

# The branches of a `Hclust` dendrogram as line segments, with the merge height on x and the
# leaves on y: leaf `hc.order[k]` sits at y = k, and each merge midway between its two
# branches. So leaf k lines up with row k of a matrix reordered by `hc.order`.
function dendrogram_segments(hc)
    y = zeros(size(hc.merges, 1))
    pos = invperm(hc.order)
    segs = Point2d[]
    for i in axes(hc.merges, 1)
        h = hc.heights[i]
        ends = map(hc.merges[i, :]) do c
            c < 0 ? (0.0, Float64(pos[-c])) : (hc.heights[c], y[c])
        end
        for (hchild, yc) in ends
            push!(segs, Point2d(hchild, yc), Point2d(h, yc))            # branch to the merge
        end
        push!(segs, Point2d(h, ends[1][2]), Point2d(h, ends[2][2]))     # the merge itself
        y[i] = (ends[1][2] + ends[2][2]) / 2
    end
    segs
end

# (2) PRIMARY VIEW. Complete-linkage hierarchical clustering of the SOS-pattern distances,
# drawn as a dendrogram-ordered |r| heatmap. Genuine groups (if any) appear as dense blocks
# on the diagonal; the orthogonal/disjoint background stays uniform. Complete linkage is the
# conservative default - it groups only all-pairs-similar nodes and will not chain marginal
# pairs. Cut the tree at |r| >= `simcut` (default 0.7, i.e. distance height 1 - simcut).
# Layout is the standard clustermap: a horizontal dendrogram on the LEFT and the heatmap on
# the RIGHT, sharing the y-axis (linked), so the leaves stay aligned with the heatmap rows;
# the colourbar only takes width from the heatmap. The dendrogram's x-axis is reversed to put
# the leaves (height 0) against the heatmap, and the heatmap's y labels print in the gap
# between the two panels, so every dendrogram tip can be read straight off as a node name. The
# x-axis carries the same names (rotated). The cut clusters with more than one node are
# outlined on the diagonal with their number and colour from `plot_cluster_tree`. Tune
# `labelsize`/`figsize` for iterative exploration.
# Returns (figure, hclust, groups::Dict node=>cluster, ordered_nodes) - `ordered_nodes` is the
# leaf order, bottom-to-top, shown on the axes.
function sos_cluster_heatmap(res, nodes, title; minoverlap, method = :pearson,
                             simcut = 0.7, linkage = :complete,
                             labelsize = 9, figsize = (1300, 1150))
    D    = sos_distances(res, nodes; minoverlap, method)
    hc   = hclust(D; linkage)
    ord  = hc.order
    grp  = cutree(hc; h = 1 - simcut)
    labs = nodes[ord]                             # node names in dendrogram-leaf order
    S    = (1 .- D)[ord, ord]                     # similarity |r|, reordered to match
    n    = length(labs)
    groups = Dict(node => grp[i] for (i, node) in enumerate(nodes))

    fig  = Figure(; size = figsize)
    Label(fig[0, 1:3], title; fontsize = 18, font = :bold)
    dend = Axis(fig[1, 1]; xlabel = "1 - |r|", xreversed = true, xgridvisible = false,
                ygridvisible = false, yticksvisible = false, yticklabelsvisible = false,
                leftspinevisible = false, topspinevisible = false, rightspinevisible = false)
    hm   = Axis(fig[1, 2]; xticks = (1:n, labs), yticks = (1:n, labs),
                xticklabelrotation = pi / 2, xticklabelsize = labelsize,
                yticklabelsize = labelsize)
    linesegments!(dend, dendrogram_segments(hc); color = :black)
    h = heatmap!(hm, 1:n, 1:n, S; colormap = :viridis, colorrange = (0, 1))
    Colorbar(fig[1, 3], h; label = "|r|")

    idmap = cluster_idmap(groups)
    colors = clustercolors(length(idmap))
    for (c, i) in idmap                           # cutree clusters are contiguous in `ord`
        rows = findall(==(c), grp[ord])
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
    xlims!(dend, 0, 1.02maximum(hc.heights))
    limits!(hm, 0.5, n + 0.5, 0.5, n + 0.5)
    colsize!(fig.layout, 1, Relative(0.20))
    fig, hc, groups, labs
end

# Cluster-size table for a `groups` Dict: how many nodes fall in each cut cluster, and which
# clusters are non-trivial (size > 1). A near-flat table of singletons is the spec's "largely
# idiosyncratic" null result; a few multi-node clusters are the co-patterned exceptions.
function sos_cluster_sizes(groups)
    counts = Dict{Int,Int}()
    for c in values(groups); counts[c] = get(counts, c, 0) + 1; end
    nontrivial = [c for (c, k) in counts if k > 1]
    members = Dict(c => sort([n for (n, g) in groups if g == c]) for c in nontrivial)
    (; nclusters = length(counts), nsingletons = count(==(1), values(counts)), members)
end

# Map the heatmap clusters onto the phylogeny: a marker at each node of a non-trivial
# cluster (size > 1), one colour per cluster, and nothing elsewhere, so the co-patterned
# groups stand out against the tree; idiosyncratic singletons are left unmarked. Pass the
# `groups` Dict that `sos_cluster_heatmap` returned. The cluster numbers and colours are
# those outlined on the heatmap (`cluster_idmap`).
function plot_cluster_tree(tree, groups, title; markersize = 12, kw...)
    idmap = cluster_idmap(groups)
    shown = Dict(n => idmap[c] for (n, c) in groups if haskey(idmap, c))
    fig, ax, p = treeplot(tree; treetype = :fan, showtips = false, nodegroup = shown,
                          groupcolors = clustercolors(length(idmap)), markersize,
                          strokewidth = 0.5, strokecolor = :gray30, axis = (; title),
                          figure = (; size = (900, 800)), kw...)
    Legend(fig[1, 2], ax, "cluster")
    fig
end

# (3) SECONDARY / CONFIRMATORY. Thresholded similarity graph + connected components. Build an
# edge only between pairs with |r| >= `simthresh`; the overlap floor is already enforced
# inside `sos_distances` (thin-overlap and disjoint pairs sit at distance 1, |r| = 0, so they
# never become edges). Co-patterned nodes then fall out as connected components. Run this AFTER
# the heatmap, once it shows there is block structure worth resolving. Returns the non-singleton
# components as vectors of node names.
function sos_similarity_communities(res, nodes; minoverlap, method = :pearson, simthresh = 0.7)
    D = sos_distances(res, nodes; minoverlap, method)
    n = length(nodes)
    g = Graphs.SimpleGraph(n)
    for i in 1:n, j in i+1:n
        (1 - D[i, j]) >= simthresh && Graphs.add_edge!(g, i, j)
    end
    comps = filter(c -> length(c) > 1, Graphs.connected_components(g))
    [[nodes[i] for i in c] for c in comps]
end

# --- Run both spaces on the divergent set (no cross-space magnitude comparison) ----------
mds_eig_g = sos_mds_eigenvalues(res_g, divergent, "MDS eigenvalues (geographic)"; minoverlap = MINOVERLAP_G)
mds_eig_e = sos_mds_eigenvalues(res_e, divergent, "MDS eigenvalues (environmental)"; minoverlap = MINOVERLAP_E)

heat_g, hc_g, groups_g, order_g = sos_cluster_heatmap(res_g, divergent, "SOS clusters (geographic)"; minoverlap = MINOVERLAP_G)
heat_e, hc_e, groups_e, order_e = sos_cluster_heatmap(res_e, divergent, "SOS clusters (environmental)"; minoverlap = MINOVERLAP_E)

communities_g = sos_similarity_communities(res_g, divergent; minoverlap = MINOVERLAP_G)
communities_e = sos_similarity_communities(res_e, divergent; minoverlap = MINOVERLAP_E)

# Clusters mapped back onto the phylogeny (numbers and colours match the heatmap outlines)
tree_clusters_g = plot_cluster_tree(tree, groups_g, "Geographic SOS clusters on the phylogeny")
tree_clusters_e = plot_cluster_tree(tree, groups_e, "Environmental SOS clusters on the phylogeny")

# Report the grouping result directly (this IS the scientific output - see spec "Reporting"):
# mostly singletons with a few multi-node clusters = "largely idiosyncratic, with named
# co-patterned exceptions"; substantial blocks = real groups to map onto the phylogeny.
sizes_g = sos_cluster_sizes(groups_g)
sizes_e = sos_cluster_sizes(groups_e)
@info "geographic SOS clusters"    sizes_g.nclusters sizes_g.nsingletons sizes_g.members communities_g
@info "environmental SOS clusters" sizes_e.nclusters sizes_e.nsingletons sizes_e.members communities_e


### --- Two node-level views of the divergent set ------------------------------------------

# (Task 1) Fan tree showing ONLY the divergent `nodes`, each labelled with its name in a small
# pale box so it stays readable over the branches; the rest of the tree is a plain grey
# skeleton. The boxes are centred on the nodes (`node_points` of the tree plot) and sized to
# their labels. They will overlap if packed too tightly, so widen `figsize` (or drop
# `fontsize`) until they clear; `stripprefix` drops the redundant "Node " so the boxes stay small.
function plot_divergent_tree(tree, nodes; fontsize = 11, figsize = (1600, 1600), stripprefix = true)
    labels = stripprefix ? replace.(string.(nodes), "Node " => "") : string.(nodes)
    fig, ax, p = treeplot(tree; treetype = :fan, showtips = false, branchcolor = :gray75,
                          figure = (; size = figsize))
    idx = [p.tree_layout[].index[n] for n in nodes]
    textlabel!(ax, p.node_points[][idx]; text = labels, fontsize,
               background_color = (:lightyellow, 0.85), strokecolor = :gray40,
               strokewidth = 0.5, padding = 2)
    fig
end

# (Task 2) One 4-panel node panel (parent / SOS / child 1 / child 2) per node, written as a
# single multi-page PDF, one node per page. SOS panel comes from the cached `res` - no
# recompute. The panel is built once and switched from node to node; each page is saved with
# CairoMakie (GLMakie cannot write PDF) and the pages merged with `pdfunite` (poppler).
function plot_node_pdf(assemblage, tree, nodes, res, outfile)
    Sys.which("pdfunite") === nothing &&
        error("plot_node_pdf needs `pdfunite` (poppler) on PATH - install it (e.g. `brew install poppler`)")
    tmp = mktempdir()
    pages = String[]
    fig, panel = nodepanel(assemblage, tree, first(nodes), res)
    for (i, node) in enumerate(nodes)
        panel.node[] = node
        page = joinpath(tmp, string(lpad(i, 3, '0'), ".pdf"))
        save(page, fig; backend = CairoMakie)
        push!(pages, page)
    end
    run(`pdfunite $pages $outfile`)
    rm(tmp; recursive = true)
    @info "wrote node-panel PDF" outfile npages = length(pages)
    outfile
end

divergent_tree = plot_divergent_tree(tree, divergent)
plot_node_pdf(birds_g, tree, divergent, res_g, "figures/divergent_node_panels_geo.pdf")
plot_node_pdf(birds_e, tree, divergent, res_e, "figures/divergent_node_panels_env.pdf")
