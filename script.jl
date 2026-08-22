# Node-based analysis of bird diversity, in environmental (birds_e) and geographic

# (birds_g) space. Run preprocess.jl first to build the cleaned inputs in
# data/clean/; this script loads them, builds the assemblages, computes and caches
# the node analysis, and explores the results.

using CSV, DataFrames, SpatialEcology, Phylo, Plots
using MultivariateStats, Statistics, JLD2, LogExpFunctions, GLM
using Nodiv

default(color = cgrad(:Spectral, rev = true))

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
plot(birds_e)

birds_g = Assemblage(phylocom_g, coords_g)
addsitestats!(birds_g, sitestats_g, :ID_geo)   # CHELSA bioclim, PC1-3, area, ...
plot(birds_g)

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

# GND of just the divergent nodes mapped onto the tree (plot_gnd marks every node
# in the Dict it is given, so pass the divergent subset rather than the full result)
plot_gnd(tree, Dict(n => e_metric[n] for n in divergent_e))
plot_gnd(tree, Dict(n => g_metric[n] for n in divergent_g))

# SOS of the most divergent node mapped onto each space (cached SOS, no recompute)
focal_e = argmax(n -> e_metric[n], divergent_e)
plot(res_e.sos[focal_e], birds_e, fillcolor = :RdYlBu, clim = (-8, 8), title = "env SOS - $focal_e")
focal_g = argmax(n -> g_metric[n], divergent_g)
plot(res_g.sos[focal_g], birds_g, fillcolor = :RdYlBu, clim = (-8, 8), title = "geo SOS - $focal_g")

# parent/SOS/children panel for that node (4th arg = cached SOS, no recompute)
plot_node(birds_e, tree, focal_e, res_e)
plot_node(birds_g, tree, focal_g, res_g)

# ordinate the divergent nodes by SOS-pattern similarity (cached SOS -> distances
# from Nodiv -> MDS; presentation stays here)
function sos_mds_plot(res, nodes, title)
    coords = predict(fit(MDS, sos_distances(res, nodes); distances = true, maxoutdim = 2))
    scatter(coords[1, :], coords[2, :], label = "",
            series_annotations = text.(nodes, 6, :bottom),
            xlabel = "MDS axis 1", ylabel = "MDS axis 2", title = title)
end
sos_mds_plot(res_e, divergent, "SOS-pattern similarity (environmental)")
sos_mds_plot(res_g, divergent, "SOS-pattern similarity (geographic)")

focal = "Node 17672"
plot_node(birds_e, tree, focal, res_e)
plot_node(birds_g, tree, focal, res_g)

same = divergent_e ∩ divergent_g

nodes = collect(keys(e_metric))
dat = DataFrame(
    :logit_g => [log(g_metric[n]) for n in nodes], #NB logit or log, depends on metric
    :logit_e => [log(e_metric[n]) for n in nodes]
)
dat = filter(row -> all(x -> !ismissing(x) && isfinite(x), row), dat)

scatter(dat.logit_g, dat.logit_e,
        xlabel = "geo $metric", ylabel = "env $metric", label = "")

mod = lm(@formula(logit_e ~ logit_g), dat)


sizes = Dict(node => noccupied(get_clade(birds_e, tree, node)) for node in nodes)
histogram(collect(values(sizes)))

plot(tree, treetype = :fan, marker_z = sizes, showtips = false, msw = 0)

scatter([sizes[n] for n in nodes], [e_metric[n] for n in nodes],
        xlabel = "occupied env sites", ylabel = "env $metric", label = "")


nodesizes = Dict(node => nspecies(get_clade(birds_g, tree, node)) for node in nodes)
scatter([log(nodesizes[n]) for n in nodes], [g_metric[n] for n in nodes],
        xlabel = "number of species in clade", ylabel = "geo $metric", label = "")



plot!([-2, 4], [-2, 4], c = :red, label = "")


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

using Clustering, StatsPlots, Graphs, LinearAlgebra
using Plots.PlotMeasures: mm     # margin units (mm) for the labelled heatmap

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
    bar(1:length(λ), λ, label = "", xlabel = "MDS axis", ylabel = "eigenvalue",
        title = "$title  (n = $(length(nodes)))")
end

# (2) PRIMARY VIEW. Complete-linkage hierarchical clustering of the SOS-pattern distances,
# drawn as a dendrogram-ordered |r| heatmap. Genuine groups (if any) appear as dense blocks
# on the diagonal; the orthogonal/disjoint background stays uniform. Complete linkage is the
# conservative default - it groups only all-pairs-similar nodes and will not chain marginal
# pairs. Cut the tree at |r| >= `simcut` (default 0.7, i.e. distance height 1 - simcut).
# Layout is the standard clustermap: a horizontal dendrogram on the LEFT and the heatmap on
# the RIGHT, sharing the y-axis. This keeps the leaves aligned with the heatmap rows whatever
# the colourbar does - the colourbar only steals width from the heatmap and so can never shift
# the row correspondence (the failure mode of stacking the dendrogram on top). `xflip` puts the
# leaves (height 0) hard against the heatmap, and the heatmap's y labels print in the gap
# between the two panels, so every dendrogram tip can be read straight off as a node name. The
# x-axis carries the same names (rotated). Tune `labelsize`/`figsize` for iterative exploration.
# Returns (plot, hclust, groups::Dict node=>cluster, ordered_nodes) - `ordered_nodes` is the
# leaf order, top-to-bottom, shown on the axes.
function sos_cluster_heatmap(res, nodes, title; minoverlap, method = :pearson,
                             simcut = 0.7, linkage = :complete,
                             labelsize = 5, figsize = (1300, 1150))
    D    = sos_distances(res, nodes; minoverlap, method)
    hc   = hclust(D; linkage)
    ord  = hc.order
    grp  = cutree(hc; h = 1 - simcut)
    labs = nodes[ord]                             # node names in dendrogram-leaf order
    S    = (1 .- D)[ord, ord]                     # similarity |r|, reordered to match
    n    = length(labs)
    # `orientation = :horizontal` drives StatsPlots' dendrogram recipe (it emits a harmless
    # "orientation is deprecated" notice from Plots - the recipe still consumes it). That recipe
    # also defaults the height axis to the SUM of merge heights, squashing the tree into a sliver,
    # so set `xlims` to the actual root height instead. `xflip` then puts the leaves (height 0)
    # against the heatmap.
    dend = plot(hc; orientation = :horizontal, xflip = true, yticks = false,
                xlims = (0, 1.02maximum(hc.heights)),
                xlabel = "1 - |r|", linecolor = :black, legend = false,
                left_margin = 4mm, bottom_margin = (6 + 1.6labelsize)mm)
    hm   = heatmap(S; c = :viridis, clims = (0, 1), colorbar_title = "|r|",
                   xticks = (1:n, labs), yticks = (1:n, labs), xrotation = 90,
                   tickfontsize = labelsize, left_margin = (4 + 1.4labelsize)mm,
                   bottom_margin = (6 + 1.6labelsize)mm)
    p = plot(dend, hm; layout = Plots.grid(1, 2, widths = [0.20, 0.80]), link = :y,
             size = figsize, plot_title = title)
    groups = Dict(node => grp[i] for (i, node) in enumerate(nodes))
    p, hc, groups, labs
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

# Map the heatmap clusters onto the phylogeny. Reuses Nodiv's `plot_gnd`, which draws a marker
# at every node in the Dict it is given and nothing elsewhere. Only the non-trivial clusters
# (size > 1) are drawn, each a distinct colour; idiosyncratic singletons are left unmarked so
# the co-patterned groups stand out against the tree. Pass the `groups` Dict that
# `sos_cluster_heatmap` returned. Cluster ids are relabelled 1..m for a compact colour scale -
# the same id labels the heatmap cut and this tree, so the two views cross-reference directly.
function plot_cluster_tree(tree, groups, title; markersize = 9, kw...)
    members = sos_cluster_sizes(groups).members          # cluster_id => member nodes (size > 1)
    ids     = sort(collect(keys(members)))
    idmap   = Dict(c => i for (i, c) in enumerate(ids))   # compact 1..m for the categorical scale
    shown   = Dict(n => float(idmap[c]) for c in ids for n in members[c])
    m       = length(ids)
    plot_gnd(tree, shown; color = cgrad(:tab20, max(m, 2); categorical = true),
             clim = (0.5, m + 0.5), colorbar_title = "cluster", title = title,
             markersize, kw...)                            # markersize now overridable (recipe fix)
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

# Clusters mapped back onto the phylogeny (colours match the heatmap cut ids above)
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

# (Task 1) Fan tree showing ONLY the divergent `nodes`, each labelled with its name on a small
# pale background so it stays readable over the branches; the rest of the tree is a plain
# skeleton. Labels are squares centred on each node (positions computed in the same fan
# coordinates the Phylo recipe uses: radius = node height, angle from node depth). The labels
# will overlap if packed too tightly, so widen `figsize` (or drop `fontsize`) until they clear;
# `stripprefix` drops the redundant "Node " so the boxes stay small.
function plot_divergent_tree(tree, nodes; fontsize = 6, boxsize = 18,
                             figsize = (1600, 1600), stripprefix = true)
    height, depth, _ = Phylo._findxy(tree)               # radius = height, angle from depth
    nleaves = length(getleafnames(tree))
    ang(node) = 2pi * depth[node] / (nleaves + 1)
    xs = [height[node] * cos(ang(node)) for node in nodes]
    ys = [height[node] * sin(ang(node)) for node in nodes]
    labels = stripprefix ? replace.(string.(nodes), "Node " => "") : string.(nodes)
    plt = plot(tree; treetype = :fan, showtips = false, linecolor = :gray75,
               legend = false, colorbar = false, size = figsize)
    scatter!(plt, xs, ys; markershape = :rect, markersize = boxsize,
             markercolor = :lightyellow, markeralpha = 0.85, markerstrokecolor = :gray40,
             markerstrokewidth = 0.3, label = "",
             series_annotations = text.(labels, fontsize, :center, :black))
    plt
end

# (Task 2) One 4-panel `plot_node` (parent / SOS / child 1 / child 2) per node, written as a
# single multi-page PDF, one node per page. SOS panel comes from the cached `res` - no
# recompute. Pages are rendered individually and merged with `pdfunite` (poppler).
function plot_node_pdf(assemblage, tree, nodes, res, outfile)
    Sys.which("pdfunite") === nothing &&
        error("plot_node_pdf needs `pdfunite` (poppler) on PATH - install it (e.g. `brew install poppler`)")
    tmp = mktempdir()
    pages = String[]
    for (i, node) in enumerate(nodes)
        page = joinpath(tmp, string(lpad(i, 3, '0'), ".pdf"))
        savefig(plot_node(assemblage, tree, node, res), page)
        push!(pages, page)
    end
    run(`pdfunite $pages $outfile`)
    rm(tmp; recursive = true)
    @info "wrote node-panel PDF" outfile npages = length(pages)
    outfile
end

divergent_tree = plot_divergent_tree(tree, divergent)
plot_node_pdf(birds_g, tree, divergent, res_g, "divergent_node_panels_geo.pdf")
plot_node_pdf(birds_e, tree, divergent, res_e, "divergent_node_panels_env.pdf")