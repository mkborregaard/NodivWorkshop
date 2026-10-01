# Open the trait explorer from the cached node analysis, in its own window: the tree, the
# SOS of the node shown in geographic and environmental space, its two child clades in trait
# space (PCA axes 2-3 of the AVONET morphometrics) and in environmental space, the latter
# with the overlap of their kernel densities. The nodes of the tree are coloured by which of
# geography and environment (rms > 1.5) and traits (TPD overlap < 0.2) their child clades
# diverge in; nodes that diverge in none of them are not marked. The script ends when the
# window is closed. Needs the same setup and caches as explorer.jl.
#
#     julia traitsexplorer.jl

using Pkg: Pkg
Pkg.activate(@__DIR__; io=devnull)
cd(@__DIR__)

using DataFrames
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

tree, birds_e, birds_g = load("data/objects.jld2", "tree", "birds_e", "birds_g")
res_e, res_g = load("data/node_analysis.jld2", "res_e", "res_g")

# The same settings as the trait explorer in traitsscript.jl
pcs, _, pca_explained = trait_pcs(birds_g)
addtraits!(birds_g, pcs, :species)
overlap = trait_overlaps(tree, trait_points(birds_g, :pca2, :pca3), collect(keys(res_e.rms)))
env_overlap = site_overlaps(tree, birds_e, collect(keys(res_g.rms)))
divergence = divergence_classes(
    res_g.rms, res_e.rms, overlap; threshold=1.5, overlap_threshold=0.2
)
const IMAGEDIR = "bow_images/workshop_species"
fig, explorer_tree = trait_explorer(
    tree,
    filter(p -> last(p) > 1, divergence),
    res_g.rms,
    (birds_g, res_g),
    (birds_e, res_e),
    pca_explained,
    overlap;
    metric=:rms,
    pcs=(2, 3),
    env_overlap,
    classes=DIVERGENCE_CLASSES,
    classlabel="divergent in (rms > 1.5, trait overlap < 0.2)",
    images=isdir(IMAGEDIR) ? IMAGEDIR : nothing,
    imageoptions=(; whitebackground=true),
)
taxa_hover!(explorer_tree, birds_g, tree)
taxa_image_hover!(explorer_tree, birds_g)

wait(display(GLMakie.Screen(; title="Trait space"), fig))
