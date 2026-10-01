# Open the trait explorer from the cached node analysis, in its own window: the tree, the
# SOS of the node shown in geographic and environmental space, its two child clades in trait
# space (PCA axes 1-2 of the AVONET morphometrics), and every node's geographic rms against
# the trait overlap (TPD) of its two child clades, where clicking a node shows it. The
# script ends when the window is closed. Needs the same setup and caches as explorer.jl.
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
overlap = trait_overlaps(tree, trait_points(birds_g, :pca1, :pca2), collect(keys(res_e.rms)))
divergent = divergent_nodes(res_e; by=:rms, threshold=2) ∪
    divergent_nodes(res_g; by=:rms, threshold=2)
const IMAGEDIR = "bow_images/workshop_species"
fig, explorer_tree = trait_explorer(
    tree,
    Dict(n => res_g.rms[n] for n in divergent),
    res_g.rms,
    (birds_g, res_g),
    (birds_e, res_e),
    pca_explained,
    overlap;
    metric=:rms,
    images=isdir(IMAGEDIR) ? IMAGEDIR : nothing,
    imageoptions=(; whitebackground=true),
)
taxa_hover!(explorer_tree, birds_g, tree)
taxa_image_hover!(explorer_tree, birds_g)

wait(display(GLMakie.Screen(; title="Trait space"), fig))
