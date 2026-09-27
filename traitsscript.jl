# The traits of the clades in the node analysis: a PCA of the AVONET morphometrics, a node
# explorer that shows each node's two child clades in trait space, and the overlap of the
# child clades in trait space against the divergence metrics. Loads the objects from
# create_objects.jl and the cached node analysis from script.jl; like script.jl, it is
# worked through interactively.

using DataFrames
using GLM
using GLMakie
using GeoInterface: GeoInterface
using GeometryOps: GeometryOps
using JLD2
using MultivariateStats
using Nodiv
using NodivMakie
using Phylo
using SpatialEcology
using Statistics
using StatsFuns

include("functions.jl")

set_theme!(; colormap=Reverse(:Spectral))

### ---- Load the objects and the node analysis ---- ###

# data/objects.jld2 is written by create_objects.jl, data/node_analysis.jld2 by script.jl
const OBJECTFILE = "data/objects.jld2"
const CACHEFILE = "data/node_analysis.jld2"
tree, birds_e, birds_g = load(OBJECTFILE, "tree", "birds_e", "birds_g")
res_e, res_g = load(CACHEFILE, "res_e", "res_g")

# The same metric and divergent nodes as in script.jl
const METRIC = :rms
const THRESHOLD = 2
metric_e = getfield(res_e, METRIC)
metric_g = getfield(res_g, METRIC)
divergent_e = divergent_nodes(res_e; by=METRIC, threshold=THRESHOLD)
divergent_g = divergent_nodes(res_g; by=METRIC, threshold=THRESHOLD)
allnodes = collect(keys(metric_e))

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

# The species images, as in script.jl's explorers, only if that folder is there
const IMAGEDIR = "bow_images/workshop_species"
trait_marked = Dict(n => metric_g[n] for n in divergent_e ∪ divergent_g)
trait_fig, trait_explorer_tree = trait_explorer(
    tree,
    trait_marked,
    metric_g,
    (birds_g, res_g),
    (birds_e, res_e),
    pca_explained;
    metric=METRIC,
    images=isdir(IMAGEDIR) ? IMAGEDIR : nothing,
    imageoptions=(; whitebackground=true),
)
taxa_hover!(trait_explorer_tree, birds_g, tree)
taxa_image_hover!(trait_explorer_tree, birds_g)

# The explorer in its own window (NODIVWORKSHOP_WINDOWS=false skips this, e.g. headless)
const SHOW_WINDOWS = isinteractive() && get(ENV, "NODIVWORKSHOP_WINDOWS", "true") != "false"
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
