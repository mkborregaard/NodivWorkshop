# Open the interactive node explorers for the geographic and the environmental space from
# the cached node analysis, each in its own window and linked: a node picked in one is shown
# in the other too, if it has an SOS there. The script ends when both windows are closed.
# Needs the setup in the README and the caches data/objects.jld2, written by
# create_objects.jl, and data/node_analysis.jld2, written by script.jl.
#
#     julia explorer.jl

using Pkg: Pkg
Pkg.activate(@__DIR__; io=devnull)
cd(@__DIR__)

using DataFrames
using GLMakie
using JLD2
using Nodiv
using NodivMakie
using Phylo
using SpatialEcology

include("functions.jl")

set_theme!(; colormap=Reverse(:Spectral))

tree, birds_e, birds_g = load("data/objects.jld2", "tree", "birds_e", "birds_g")
res_e, res_g = load("data/node_analysis.jld2", "res_e", "res_g")

# The same settings as the explorers in script.jl
const IMAGEDIR = "bow_images/workshop_species"
function explorer(birds, res)
    return node_explorer(
        birds,
        tree,
        res;
        metric=:rms,
        nodes=divergent_nodes(res; by=:rms, threshold=2),
        images=isdir(IMAGEDIR) ? IMAGEDIR : nothing,
        imageoptions=(; whitebackground=true),
        ordinationkw=(; minoverlap=3),
    )
end
explorer_fig_e, explorer_e = explorer(birds_e, res_e)
explorer_fig_g, explorer_g = explorer(birds_g, res_g)
link_explorers!(tree, explorer_e, explorer_g)
for (ex, birds) in ((explorer_e, birds_e), (explorer_g, birds_g))
    taxa_hover!(ex, birds, tree)
    taxa_image_hover!(ex, birds)
end

screens = [
    display(GLMakie.Screen(; title="Environmental space"), explorer_fig_e),
    display(GLMakie.Screen(; title="Geographic space"), explorer_fig_g),
]
foreach(wait, screens)
