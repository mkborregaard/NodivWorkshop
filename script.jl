# Node-based analysis of bird diversity, in environmental (birds_e) and geographic
# (birds_g) space. Run preprocess.jl first to build the cleaned inputs in
# data/clean/; this script loads them, builds the assemblages, computes and caches
# the node analysis, and explores the results.
#
# Plotting is Makie: NodivMakie for trees, maps and node panels, GLMakie for interactive
# windows. Run it in the REPL or VS Code; each figure is kept in a variable, so evaluate
# the variable to show it again.

using CSV
# CairoMakie only for saving vector (PDF) files; GLMakie saves raster formats
using CairoMakie: CairoMakie
using DataFrames
using GLM
using GLMakie
using GeoInterface: GeoInterface
using GeometryOps: GeometryOps
using JLD2
using LogExpFunctions
using MultivariateStats
using Nodiv
using NodivMakie
using Phylo
using SpatialEcology
using Statistics
using StatsFuns
# Loading a backend activates it, so make sure GLMakie is the one
GLMakie.activate!()

include("functions.jl")

set_theme!(; colormap=Reverse(:Spectral))

### ---- Load the cleaned inputs (from preprocess.jl) and build the assemblages ---- ###

tree = sort!(parsenewick(read("data/clean/tree.nwk", String)))
phylocom_e = string_sites!(CSV.read("data/clean/phylocom_e.csv", DataFrame))
coords_e = string_sites!(CSV.read("data/clean/coords_e.csv", DataFrame))
sitestats_e = CSV.read("data/clean/sitestats_e.csv", DataFrame)
phylocom_g = string_sites!(CSV.read("data/clean/phylocom_g.csv", DataFrame))
coords_g = string_sites!(CSV.read("data/clean/coords_g.csv", DataFrame))
sitestats_g = CSV.read("data/clean/sitestats_g.csv", DataFrame)
sitestats_g.ID_geo = string.(sitestats_g.ID_geo)
avonet = CSV.read("data/clean/traits.csv", DataFrame)

# coordinates were pre-aligned to each phylocom's site order in preprocessing, so
# they slot straight into the Assemblage (SpatialEcology aligns coords by row order).
birds_e = Assemblage(phylocom_e, coords_e)
addsitestats!(birds_e, sitestats_e, :ID_env)   # PC bins, area, occupancy, ...
addtraits!(birds_e, avonet, :species)
richness_e = map_figure(birds_e; title="Environmental: species richness", label="species")
# save("figures/Env species richness.png", richness_e)

birds_g = Assemblage(phylocom_g, coords_g)
addsitestats!(birds_g, sitestats_g, :ID_geo)   # CHELSA bioclim, PC1-3, area, ...
addtraits!(birds_g, avonet, :species)
richness_g = map_figure(
    birds_g;
    title="Geographic: species richness",
    label="species",
    figure=(; size=(1000, 500)),
)
# save("figures/Geo species richness.png", richness_g)

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
const SIMCUT = 0.7

### ---- Exploratory plotting (from the cached NodeMetrics; `_e` vs `_g`) ---- ###

# Strongly divergent nodes in each space (`METRIC` above `THRESHOLD`)
divergent_e = divergent_nodes(res_e; by=METRIC, threshold=THRESHOLD)
divergent_g = divergent_nodes(res_g; by=METRIC, threshold=THRESHOLD)
divergent = divergent_e ∩ divergent_g

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
# save("figures/Env divergent nodes treeplot.png", metric_tree_e)
metric_tree_g = metric_tree(
    tree,
    res_g;
    metric=METRIC,
    nodes=divergent_g,
    title="Geographic: divergent nodes, $METRIC",
)
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
# save("figures/Geo SOS $focal_g.png", sosmap_g)

# The interactive entry point: the fan tree with the divergent nodes marked, next to the
# SOS map and the two child clades' maps, and an ordination of the divergent nodes by
# SOS-pattern similarity. Each opens on its space's most divergent node (focal_e, focal_g
# above). Click a node or a branch on the tree, or a point in the ordination, to show that
# node; hover for labels. The Birds of the World images in bow_images/ are private and not
# in the repo: they are used only if that folder is there.
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

# Each explorer in its own window (NODIVWORKSHOP_WINDOWS=false skips this, e.g. headless)
const SHOW_WINDOWS = isinteractive() && get(ENV, "NODIVWORKSHOP_WINDOWS", "true") != "false"
if SHOW_WINDOWS
    display(GLMakie.Screen(), explorer_fig_e)
    display(GLMakie.Screen(), explorer_fig_g)
end

# Ordinate the divergent nodes of both spaces by SOS-pattern similarity (cached SOS ->
# `sos_distances` -> classical MDS in `sos_ordination`, both from Nodiv)
D_g = sos_distances(res_g, divergent; minoverlap=MINOVERLAP)
D_e = sos_distances(res_e, divergent; minoverlap=MINOVERLAP)
function sos_mds_scatter(D, nodes, title)
    return ordinationplot(sos_ordination(D, nodes); nodelabels=true, axis=(; title)).figure
end
mds_e = sos_mds_scatter(D_e, divergent, "Environmental: SOS-pattern similarity")
# save("figures/Env SOS-pattern similarity.png", mds_e)
mds_g = sos_mds_scatter(D_g, divergent, "Geographic: SOS-pattern similarity")
# save("figures/Geo SOS-pattern similarity.png", mds_g)

# Parent/SOS/children panel for one node (4th arg = cached SOS, no recompute); also
# `explorer_e.panel.node[] = focal` shows it in the explorer
# Node names are numbered by data/clean/tree.nwk: re-running preprocess.jl renumbers them
focal = "Node 17672"
panel_e, _ = node_panel(birds_e, tree, focal, res_e)
# save("figures/Env node panel $focal.png", panel_e)
panel_g, _ = node_panel(birds_g, tree, focal, res_g)
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
# save("figures/Env vs geo $METRIC.png", metric_scatter)

rms_fit = lm(@formula(log_e ~ log_g), dat)

occupied_of = clade_richness(birds_e, tree)
occupied_e = Dict(node => count(>(0), occupied_of(node)) for node in allnodes)
occupied_hist = hist(
    collect(values(occupied_e)); axis=(; xlabel="occupied env sites", ylabel="nodes")
)
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
# save("figures/Env occupied sites treeplot.png", occupied_tree)

occupied_scatter = scatter(
    [occupied_e[n] for n in allnodes],
    [metric_e[n] for n in allnodes];
    axis=(; xlabel="occupied env sites", ylabel="env $METRIC"),
)
# save("figures/Env $METRIC vs occupied sites.png", occupied_scatter)

nspecies_g = Dict(node => length(nodespecies(tree, node)) for node in allnodes)
nspecies_scatter = scatter(
    [log(nspecies_g[n]) for n in allnodes],
    [metric_g[n] for n in allnodes];
    axis=(; xlabel="log number of species in clade", ylabel="geo $METRIC"),
)
# save("figures/Geo $METRIC vs clade species.png", nspecies_scatter)

### ---- Grouping divergent nodes by SOS-pattern similarity ---- ###
# The MDS scatter (`sos_mds_scatter`) above is read together with its eigenvalue diagnostic
# (1); the primary read is the complete-linkage clustered heatmap (2), with a thresholded
# similarity graph (3) as the confirmatory secondary. Distances come from `sos_distances`
# in Nodiv: 1 - |r| over the shared occupied cells, with the minimum-overlap floor
# `MINOVERLAP`. Each space's distances are computed once and every view is derived from
# them. The two spaces are run separately and their magnitudes are NOT compared
# (environmental "occupancy" is over tens of PC bins, geographic over ~18k cells).

### ---- Run both spaces on the divergent set (no cross-space magnitude comparison) ---- ###
# (1) ONE-TIME DIAGNOSTIC, not the analysis. Fit MDS at a higher dimension and look at the
# eigenvalue spectrum: if axes 3+ carry weight comparable to axes 1-2, the 2-D scatter is a
# projection artefact and the "ring" is the honest report of near-equidistance.
mds_eig_g = eigenvalueplot(
    sos_ordination(D_g, divergent; maxoutdim=10);
    axis=(; title="Geographic: MDS eigenvalues (n = $(length(divergent)))"),
).figure
# save("figures/Geo MDS eigenvalues.png", mds_eig_g)
mds_eig_e = eigenvalueplot(
    sos_ordination(D_e, divergent; maxoutdim=10);
    axis=(; title="Environmental: MDS eigenvalues (n = $(length(divergent)))"),
).figure
# save("figures/Env MDS eigenvalues.png", mds_eig_e)

# (2) PRIMARY VIEW. Complete-linkage hierarchical clustering, cut at |r| >= SIMCUT, drawn
# as a dendrogram-ordered |r| heatmap
clusters_g = sos_clusters(D_g, divergent; simcut=SIMCUT)
clusters_e = sos_clusters(D_e, divergent; simcut=SIMCUT)

heat_g = sos_cluster_heatmap(clusters_g; title="Geographic: SOS clusters")
# save("figures/Geo SOS clusters.png", heat_g)
heat_e = sos_cluster_heatmap(clusters_e; title="Environmental: SOS clusters")
# save("figures/Env SOS clusters.png", heat_e)

# (3) SECONDARY / CONFIRMATORY. Modularity communities of the similarity graph with edges at
# |r| >= SIMCUT
communities_g = sos_similarity_communities(D_g, divergent; simthresh=SIMCUT)
communities_e = sos_similarity_communities(D_e, divergent; simthresh=SIMCUT)

# Clusters mapped back onto the phylogeny (numbers and colours match the heatmap outlines)
tree_clusters_g = cluster_tree(
    tree, clusters_g; title="Geographic: SOS clusters on the phylogeny"
)
# save("figures/Geo SOS clusters on the phylogeny.png", tree_clusters_g)
tree_clusters_e = cluster_tree(
    tree, clusters_e; title="Environmental: SOS clusters on the phylogeny"
)
# save("figures/Env SOS clusters on the phylogeny.png", tree_clusters_e)

# Report the grouping result directly (this IS the scientific output):
# mostly singletons with a few multi-node clusters = "largely idiosyncratic, with named
# co-patterned exceptions"; substantial blocks = real groups to map onto the phylogeny.
display(clusters_g)
display(communities_g)
display(clusters_e)
display(communities_e)

### ---- Two node-level views of the divergent set ---- ###

# Fan tree showing ONLY the divergent nodes, each labelled with its name in a small pale box
# so it stays readable over the branches; the rest of the tree is a plain grey skeleton. The
# boxes will overlap if packed too tightly, so widen the figure size (or drop
# `nodelabelsize`) until they clear; the redundant "Node " is dropped so the boxes stay
# small.
divergent_tree = treeplot(
    tree;
    treetype=:fan,
    showtips=false,
    branchcolor=:gray75,
    nodelabels=Dict(n => replace(n, "Node " => "") for n in divergent),
    nodelabelbackground=(:lightyellow, 0.85),
    nodelabelsize=11,
    nodelabelalign=(:center, :center),
    nodelabeloffset=(0, 0),
    figure=(; size=(1600, 1600)),
).figure
# save("figures/Divergent nodes labelled.png", divergent_tree)
# node_panel_pdf(
#     birds_g, tree, divergent, res_g, "figures/divergent_node_panels_geo.pdf";
#     backend=CairoMakie,
# )
# node_panel_pdf(
#     birds_e, tree, divergent, res_e, "figures/divergent_node_panels_env.pdf";
#     backend=CairoMakie,
# )

### ---- Traits: PCA of the AVONET morphometrics ---- ###

# Shapiro-Francia W': the squared correlation of the sorted values with the normal quantiles
function normality(x)
    n = length(x)
    q = norminvcdf.(((1:n) .- 0.375) ./ (n + 0.25))
    return cor(sort(x), q)^2
end

# PCA of the columns of `df`: a column is logged where that raises its W' by more than
# `tol`, then all are z-transformed. Returns the PCA, the logged columns and the
# z-transformed matrix (species x traits) it was fit on.
function trait_pca(df; tol=0.01)
    X = Matrix{Float64}(df)
    logged = [all(>(0), x) && normality(log.(x)) - normality(x) > tol for x in eachcol(X)]
    X[:, logged] .= log.(X[:, logged])
    Z = (X .- mean(X; dims=1)) ./ std(X; dims=1)
    return (; pca=fit(PCA, permutedims(Z); pratio=1), logged=names(df)[logged], z=Z)
end

speciestraits = traits(birds_g)
trait_pca_fit = trait_pca(speciestraits[:, Between(:Beak_Length_Culmen, :Mass)])
pca_explained = principalvars(trait_pca_fit.pca) ./ var(trait_pca_fit.pca)

pcs = DataFrame(
    permutedims(predict(trait_pca_fit.pca, permutedims(trait_pca_fit.z))[1:4, :]),
    ["pca$i" for i in 1:4],
)
pcs.species = speciestraits.name
addtraits!(birds_e, pcs, :species)
addtraits!(birds_g, pcs, :species)

# Convex hull of 2-d points (a GeoInterface polygon), or `nothing` for fewer than three
# distinct points
function convex_hull(pts)
    return length(unique(pts)) < 3 ? nothing : GeometryOps.convex_hull(pts)
end

# The overlap of two convex hulls as a proportion of the smaller one; NaN if either has no
# area
function hull_overlap(h1, h2)
    (h1 === nothing || h2 === nothing) && return NaN
    a = min(GeometryOps.area(h1), GeometryOps.area(h2))
    a > 0 || return NaN
    alg = GeometryOps.ConvexConvexSutherlandHodgman()
    return GeometryOps.intersection_area(alg, h1, h2) / a
end

# Species => point in trait space, from the columns `x` and `y` of an assemblage's traits
function trait_points(asm, x, y)
    t = traits(asm)
    return Dict(zip(t.name, Point2d.(t[!, x], t[!, y])))
end

# The trait-space points of the species of each of `node`'s two child clades
function child_points(tree, node, pts)
    return [
        [pts[sp] for sp in nodespecies(tree, getnodename(tree, c)) if haskey(pts, sp)] for
        c in getchildren(tree, node)[1:2]
    ]
end

function closed_hull(pts)
    h = convex_hull(pts)
    (h === nothing || GeometryOps.area(h) == 0) && return Point2d[]
    return Point2d.(GeoInterface.getpoint(GeoInterface.getexterior(h)))
end

# All species in trait space in grey, with the two child clades of `node` (an Observable) in
# the explorer's clade colours, the smaller clade on top, each outlined by its convex hull
function trait_panel!(gp, asm, tree, node, x, y; axis=(;))
    pts = trait_points(asm, x, y)
    colors = clade_colors(:RdYlBu)
    ax = Axis(gp; xgridvisible=false, ygridvisible=false, axis...)
    scatter!(ax, collect(values(pts)); color=:gray80, markersize=3, inspectable=false)
    clades = lift(n -> child_points(tree, n, pts), node)
    for (k, color) in enumerate(colors)
        cladepts = lift(c -> c[k], clades)
        sc = scatter!(ax, cladepts; color, markersize=5, inspectable=false)
        on(c -> translate!(sc, 0, 0, length(c[k]) <= length(c[3 - k])), clades; update=true)
        hull = lines!(
            ax, lift(closed_hull, cladepts); color, linewidth=2, inspectable=false
        )
        translate!(hull, 0, 0, 2)
    end
    return ax
end

# A node explorer of the two spaces and trait space: the tree, the SOS of the node shown in
# geographic and environmental space, and its two child clades on PCA axes 1-2 and 3-4.
# Each space is passed as an (assemblage, NodeMetrics) pair.
function trait_explorer(
    tree,
    marked,
    nodevalues,
    (birds_g, res_g),
    (birds_e, res_e),
    explained;
    metric,
    images=nothing,
    imageoptions=(;),
)
    fig = Figure(; size=(1600, 850))
    node = Observable(argmax(n -> marked[n], keys(marked)))
    tr = explorer_tree!(
        fig[1, 1],
        tree,
        node,
        marked;
        values=nodevalues,
        label="geo $metric",
        selectable=n -> has_sos(tree, res_g.sos, n) && has_sos(tree, res_e.sos, n),
        unselectable="no SOS in both spaces",
        images,
        imageoptions,
        rangesize=birds_g,
    )
    panels = fig[1, 2] = GridLayout()
    sos_map!(panels[1, 1], birds_g, node, res_g; title="Geographic SOS")
    sos_map!(panels[1, 2], birds_e, node, res_e; title="Environmental SOS")
    pc_label(i) = "pca$i ($(round(100explained[i]; digits = 1))%)"
    for (col, (i, j)) in enumerate(((1, 2), (3, 4)))
        trait_panel!(
            panels[2, col],
            birds_g,
            tree,
            node,
            Symbol("pca$i"),
            Symbol("pca$j");
            axis=(; xlabel=pc_label(i), ylabel=pc_label(j)),
        )
    end
    colsize!(fig.layout, 1, Relative(0.45))
    DataInspector(fig)
    return fig, tr
end

trait_marked = Dict(n => metric_g[n] for n in divergent_e ∪ divergent_g)
trait_fig, trait_explorer_tree = trait_explorer(
    tree,
    trait_marked,
    metric_g,
    (birds_g, res_g),
    (birds_e, res_e),
    pca_explained;
    metric=METRIC,
    explorer_options.images,
    explorer_options.imageoptions,
)
if SHOW_WINDOWS
    display(GLMakie.Screen(), trait_fig)
end

pts12 = trait_points(birds_g, :pca1, :pca2)
trait_overlap = Dict(
    n => hull_overlap(convex_hull.(child_points(tree, n, pts12))...) for n in allnodes
)

overlap_dat = DataFrame(;
    overlap=[trait_overlap[n] for n in allnodes],
    log_g=[log(metric_g[n]) for n in allnodes],
    log_e=[log(metric_e[n]) for n in allnodes],
)
overlap_dat = filter(row -> all(isfinite, row), overlap_dat)
overlap_fit_g = lm(@formula(log_g ~ overlap), overlap_dat)
overlap_fit_e = lm(@formula(log_e ~ overlap), overlap_dat)

overlap_scatter = let fig = Figure(; size=(1100, 500))
    panels = (
        (:log_g, overlap_fit_g, "Geographic: trait overlap and $METRIC"),
        (:log_e, overlap_fit_e, "Environmental: trait overlap and $METRIC"),
    )
    for (col, (y, lmfit, title)) in enumerate(panels)
        ax = Axis(fig[1, col]; title, xlabel="trait overlap (pca1-2)", ylabel="log $METRIC")
        scatter!(ax, overlap_dat.overlap, overlap_dat[!, y]; markersize=5)
        ablines!(ax, coef(lmfit)...; color=:red)
    end
    fig
end
# save("figures/Trait overlap vs $METRIC.png", overlap_scatter)
