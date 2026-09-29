# Node-based analysis of bird diversity, in environmental (birds_e) and geographic
# (birds_g) space. Run preprocess.jl and then create_objects.jl first; this script loads
# the objects, computes and caches the node analysis, and explores the results.
#
# Plotting is Makie: NodivMakie for trees, maps and node panels, GLMakie for interactive
# windows. Run it in the REPL or VS Code; each figure is kept in a variable, so evaluate
# the variable to show it again.

# CairoMakie only for saving vector (PDF) files; GLMakie saves raster formats
using CairoMakie: CairoMakie
using DataFrames
using GLM
using GLMakie
using JLD2
using LogExpFunctions
using Nodiv
using NodivMakie
using Phylo
using SpatialEcology
using Statistics
# Loading a backend activates it, so make sure GLMakie is the one
GLMakie.activate!()

include("functions.jl")

set_theme!(; colormap=Reverse(:Spectral))

### ---- Load the objects (from create_objects.jl) ---- ###

const OBJECTFILE = "data/objects.jld2"
tree, birds_e, birds_g = load(OBJECTFILE, "tree", "birds_e", "birds_g")

### ---- Heavy step: divergence metrics + SOS for every node, cached to disk ---- ###
# `node_metrics` computes the divergence metrics (GND, RMS-SOS, ...) together with the
# per-cell SOS they are built on. This is the slow part - randomisations over the whole
# tree, and ~18k cells for the geographic scan - so cache it: re-running the script just
# reloads the results and jumps straight to the plotting below.
const CACHEFILE = "data/node_analysis.jld2"
if !isfile(CACHEFILE)
    res_e = node_metrics(birds_e, tree; nsims=200)
    res_g = node_metrics(birds_g, tree; nsims=200)
    jldsave(CACHEFILE; res_e, res_g)
end
# Each a NodeMetrics: nodes (in tree order), gnd, rms, sd, ses, pval, varying and sos
res_e, res_g = load(CACHEFILE, "res_e", "res_g")

const METRIC = :rms  # Alternatives are :pval and :gnd
const THRESHOLD = 2
metric_e = getfield(res_e, METRIC)
metric_g = getfield(res_g, METRIC)

# Minimum number of shared occupied cells for the correlation behind `sos_distances`
# (SOS-pattern similarity, used by the explorers' ordination and the grouping section).
# Below it `sos_distances` pins the pair at distance 1 rather than trusting a correlation
# fit on a handful of cells - which is also what keeps disjoint pairs at the maximum. 3 is
# the smallest overlap where |r| is not trivially 1.
const MINOVERLAP = 3
const SIMCUT = 0.6

### ---- Exploratory plotting (from the cached NodeMetrics; `_e` vs `_g`) ---- ###

richness_e = map_figure(birds_e; title="Environmental: species richness", label="species")
richness_e
# save("figures/Env species richness.png", richness_e)

richness_g = map_figure(
    birds_g;
    title="Geographic: species richness",
    label="species",
    figure=(; size=(1000, 500)),
)
richness_g
# save("figures/Geo species richness.png", richness_g)

# Strongly divergent nodes in each space (`METRIC` above `THRESHOLD`)
divergent_e = divergent_nodes(res_e; by=METRIC, threshold=THRESHOLD)
divergent_g = divergent_nodes(res_g; by=METRIC, threshold=THRESHOLD)

# The metric of just the divergent nodes mapped onto the tree (markers only at `nodes`:
# a node missing from the Dict would get a transparent fill but still its outline).
# GND is a proportion, so it gets a fixed 0-1 colour range, as plot_gnd used.
metric_tree_e = metric_tree(
    tree,
    res_e;
    metric=METRIC,
    nodes=divergent_e,
    title="Environmental: divergent nodes, $METRIC",
)
metric_tree_e
# save("figures/Env divergent nodes treeplot.png", metric_tree_e)
metric_tree_g = metric_tree(
    tree,
    res_g;
    metric=METRIC,
    nodes=divergent_g,
    title="Geographic: divergent nodes, $METRIC",
)
metric_tree_g
# save("figures/Geo divergent nodes treeplot.png", metric_tree_g)

# SOS of the most divergent node mapped onto each space (cached SOS, no recompute)
focal_e = argmax(n -> metric_e[n], divergent_e)
sosmap_e = map_figure(
    res_e.sos[focal_e],
    birds_e;
    colormap=:RdYlBu,
    colorrange=(-8, 8),
    title="Environmental: SOS of $focal_e",
    label="SOS",
)
sosmap_e
# save("figures/Env SOS $focal_e.png", sosmap_e)
focal_g = argmax(n -> metric_g[n], divergent_g)
sosmap_g = map_figure(
    res_g.sos[focal_g],
    birds_g;
    colormap=:RdYlBu,
    colorrange=(-8, 8),
    title="Geographic: SOS of $focal_g",
    label="SOS",
    figure=(; size=(1000, 500)),
)
sosmap_g
# save("figures/Geo SOS $focal_g.png", sosmap_g)

# Parent/SOS/children panel for one node (4th arg = cached SOS, no recompute); also
# `explorer_e.panel.node[] = focal` shows it in the explorer
# Node names are numbered by data/clean/tree.nwk: re-running preprocess.jl renumbers them
focal = "Node 17672"
panel_e, _ = node_panel(birds_e, tree, focal, res_e)
panel_e
# save("figures/Env node panel $focal.png", panel_e)
panel_g, _ = node_panel(birds_g, tree, focal, res_g)
panel_g
# save("figures/Geo node panel $focal.png", panel_g)

allnodes = collect(keys(metric_e))
dat = DataFrame(;
    log_g=[log(metric_g[n]) for n in allnodes],  # NB logit or log, depends on metric
    log_e=[log(metric_e[n]) for n in allnodes],
)
dat = filter(row -> all(isfinite, row), dat)

metric_scatter = scatter(
    dat.log_g, dat.log_e; axis=(; xlabel="log geo $METRIC", ylabel="log env $METRIC")
)
ablines!(metric_scatter.axis, 0, 1; color=:red)  # The 1:1 line
metric_scatter
# save("figures/Env vs geo $METRIC.png", metric_scatter)

rms_fit = lm(@formula(log_e ~ log_g), dat)

occupied_of = clade_richness(birds_e, tree)
occupied_e = Dict(node => count(>(0), occupied_of(node)) for node in allnodes)
occupied_hist = hist(
    collect(values(occupied_e)); axis=(; xlabel="occupied env sites", ylabel="nodes")
)
occupied_hist
# save("figures/Env occupied sites histogram.png", occupied_hist)

occupied_tree = let
    fig, ax, p = treeplot(
        tree;
        treetype=:fan,
        nodecolor=occupied_e,
        showtips=false,
        markersize=5,
        figure=(; size=(800, 700)),
    )
    Colorbar(fig[1, 2], p; label="occupied env sites")
    fig
end
occupied_tree
# save("figures/Env occupied sites treeplot.png", occupied_tree)

occupied_scatter = scatter(
    [occupied_e[n] for n in allnodes],
    [metric_e[n] for n in allnodes];
    axis=(; xlabel="occupied env sites", ylabel="env $METRIC"),
)
occupied_scatter
# save("figures/Env $METRIC vs occupied sites.png", occupied_scatter)

nspecies_g = Dict(node => length(nodespecies(tree, node)) for node in allnodes)
nspecies_scatter = scatter(
    [log(nspecies_g[n]) for n in allnodes],
    [metric_g[n] for n in allnodes];
    axis=(; xlabel="log number of species in clade", ylabel="geo $METRIC"),
)
nspecies_scatter
# save("figures/Geo $METRIC vs clade species.png", nspecies_scatter)

### ---- Grouping divergent nodes by SOS-pattern similarity ---- ###
# Complete-linkage clustering of each space's divergent nodes, cut at |r| >= SIMCUT: every
# pair of nodes in a cluster has SOS maps correlated at |r| >= SIMCUT, every node is in
# exactly one cluster, and a node with no such partner is on its own. Distances come from
# `sos_distances` in Nodiv: 1 - |r| over the cells where both SOS maps are defined, with
# the minimum-overlap floor `MINOVERLAP`. The two spaces' |r| are NOT compared
# (environmental space has ~500 PC bins, geographic ~18k cells). See
# docs/sos_pattern_grouping_design.md for why the number of clusters is not chosen from
# the data.
D_g = sos_distances(res_g, divergent_g; minoverlap=MINOVERLAP)
D_e = sos_distances(res_e, divergent_e; minoverlap=MINOVERLAP)
clusters_g = sos_clusters(D_g, divergent_g; simcut=SIMCUT)
clusters_e = sos_clusters(D_e, divergent_e; simcut=SIMCUT)
display(clusters_g)
display(clusters_e)

# The divergent nodes ordinated by SOS-pattern similarity (classical MDS of the distances
# above, `sos_ordination` in Nodiv), numbered and coloured by their cluster (grey: on its
# own). A 2-D projection, so its distances are approximate; the clusters come from the full
# distances.
function sos_mds_scatter(D, nodes, clusters, title)
    return ordinationplot(
        sos_ordination(D, nodes);
        nodelabels=true,
        nodecolor=cluster_nodecolors(clusters),
        axis=(; title),
    ).figure
end
mds_g = sos_mds_scatter(
    D_g, divergent_g, clusters_g, "Geographic: SOS-pattern similarity"
)
mds_g
# save("figures/Geo SOS-pattern similarity.png", mds_g)
mds_e = sos_mds_scatter(
    D_e, divergent_e, clusters_e, "Environmental: SOS-pattern similarity"
)
mds_e
# save("figures/Env SOS-pattern similarity.png", mds_e)

# Clusters mapped back onto the phylogeny
tree_clusters_g = cluster_tree(
    tree, clusters_g; title="Geographic: SOS clusters on the phylogeny"
)
tree_clusters_g
# save("figures/Geo SOS clusters on the phylogeny.png", tree_clusters_g)
tree_clusters_e = cluster_tree(
    tree, clusters_e; title="Environmental: SOS clusters on the phylogeny"
)
tree_clusters_e
# save("figures/Env SOS clusters on the phylogeny.png", tree_clusters_e)

# The |r| of every pair in dendrogram order, clusters outlined with the numbers and colours
# of the tree above
heat_g = sos_cluster_heatmap(clusters_g; title="Geographic: SOS clusters")
heat_g
# save("figures/Geo SOS clusters.png", heat_g)
heat_e = sos_cluster_heatmap(clusters_e; title="Environmental: SOS clusters")
heat_e
# save("figures/Env SOS clusters.png", heat_e)

### ---- The node explorers ---- ###

# The interactive entry point: the fan tree with the divergent nodes marked, next to the
# SOS map and the two child clades' maps, and an ordination of the divergent nodes by
# SOS-pattern similarity, coloured by their SOS cluster above (grey: on its own). The
# ordination is a 2-D projection for browsing, so its distances are approximate; the
# clusters come from the full distances. Each opens on its space's most divergent node
# (focal_e, focal_g above). Click a node or a branch on the tree, or a point in the
# ordination, to show that node; hover for labels. The Birds of the World images in
# bow_images/ are private and not in the repo: they are used only if that folder is there.
const IMAGEDIR = "bow_images/workshop_species"
explorer_options = (;
    metric=METRIC,
    images=isdir(IMAGEDIR) ? IMAGEDIR : nothing,
    imageoptions=(; whitebackground=true),
    ordinationkw=(; minoverlap=MINOVERLAP),
)
explorer_fig_e, explorer_e = node_explorer(
    birds_e, tree, res_e; nodes=divergent_e, explorer_options...
)
explorer_fig_g, explorer_g = node_explorer(
    birds_g, tree, res_g; nodes=divergent_g, explorer_options...
)

# Link the two: a node picked in one space is shown in the other too, if it has an SOS there
link_explorers!(tree, explorer_e, explorer_g)

for (explorer, birds, clusters) in
    ((explorer_e, birds_e, clusters_e), (explorer_g, birds_g, clusters_g))
    taxa_hover!(explorer, birds, tree)
    taxa_image_hover!(explorer, birds)
    color_by_clusters!(explorer, clusters)
end

# Each explorer in its own window (NODIVWORKSHOP_WINDOWS=false skips this, e.g. headless)
const SHOW_WINDOWS = isinteractive() && get(ENV, "NODIVWORKSHOP_WINDOWS", "true") != "false"
if SHOW_WINDOWS
    display(GLMakie.Screen(), explorer_fig_e)
    display(GLMakie.Screen(), explorer_fig_g)
end

### ---- Two node-level views of each space's divergent nodes ---- ###

# Fan tree showing ONLY the divergent nodes, each labelled with its name in a small pale box
# so it stays readable over the branches; the rest of the tree is a plain grey skeleton. The
# boxes will overlap if packed too tightly, so widen the figure size (or drop
# `nodelabelsize`) until they clear; the redundant "Node " is dropped so the boxes stay
# small.
function labelled_tree(nodes, title)
    return treeplot(
        tree;
        treetype=:fan,
        showtips=false,
        branchcolor=:gray75,
        nodelabels=Dict(n => replace(n, "Node " => "") for n in nodes),
        nodelabelbackground=(:lightyellow, 0.85),
        nodelabelsize=11,
        nodelabelalign=(:center, :center),
        nodelabeloffset=(0, 0),
        axis=(; title),
        figure=(; size=(1600, 1600)),
    ).figure
end
divergent_tree_g = labelled_tree(divergent_g, "Geographic: divergent nodes")
# save("figures/Geo divergent nodes labelled.png", divergent_tree_g)
divergent_tree_e = labelled_tree(divergent_e, "Environmental: divergent nodes")
# save("figures/Env divergent nodes labelled.png", divergent_tree_e)
# node_panel_pdf(
#     birds_g, tree, divergent_g, res_g, "figures/divergent_node_panels_geo.pdf";
#     backend=CairoMakie,
# )
# node_panel_pdf(
#     birds_e, tree, divergent_e, res_e, "figures/divergent_node_panels_env.pdf";
#     backend=CairoMakie,
# )
