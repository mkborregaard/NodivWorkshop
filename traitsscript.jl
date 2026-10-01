# The traits of the clades in the node analysis: a PCA of the AVONET morphometrics, a node
# explorer that shows each node's two child clades in trait space, and the overlap of the
# child clades in trait space against the divergence metrics. Loads the objects from
# create_objects.jl and the cached node analysis from script.jl; like script.jl, it is
# worked through interactively.

using DataFrames
using GLM
using GLMakie
using JLD2
using MultivariateStats
using Nodiv
using NodivMakie
using Phylo
using SpatialEcology
using Statistics
using StatsFuns

include("functions.jl")
include("traitfunctions.jl")

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

pcs, trait_pca_fit, pca_explained = trait_pcs(birds_g)
addtraits!(birds_e, pcs, :species)
addtraits!(birds_g, pcs, :species)

### ---- Trait overlap (TPD) and the trait explorer ---- ###

# The two PCA axes of trait space
const TRAIT_PCS = (2, 3)
trait_pts = trait_points(birds_g, (Symbol("pca$i") for i in TRAIT_PCS)...)
trait_overlap = trait_overlaps(tree, trait_pts, allnodes)

# The species images, as in script.jl's explorers, only if that folder is there
const IMAGEDIR = "bow_images/workshop_species"
trait_divergent(t) = divergent_nodes(res_e; by=METRIC, threshold=t) ∪
    divergent_nodes(res_g; by=METRIC, threshold=t)
trait_marked(t) = Dict(n => metric_g[n] for n in trait_divergent(t))
trait_fig, trait_explorer_tree = trait_explorer(
    tree,
    trait_marked,
    metric_g,
    (birds_g, res_g),
    (birds_e, res_e),
    pca_explained,
    trait_overlap;
    metric=METRIC,
    pcs=TRAIT_PCS,
    threshold=THRESHOLD,
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
        xlabel = "TPD overlap (pca$(TRAIT_PCS[1])-$(TRAIT_PCS[2]))"
        ax = Axis(fig[1, col]; title, xlabel, ylabel="log $METRIC")
        scatter!(ax, overlap_dat.overlap, overlap_dat[!, y]; markersize=5)
        ablines!(ax, coef(lmfit)...; color=:red)
    end
    fig
end
# save("figures/Trait overlap vs $METRIC.png", overlap_scatter)
