# Build the analysis objects from the cleaned inputs in data/clean/ (written by
# preprocess.jl) and cache them in data/objects.jld2, which script.jl, traitsscript.jl
# and explorer.jl load:
#   - tree               the phylogeny, ladderized
#   - birds_e, birds_g   the assemblages of each space (e = env, g = geo), with their site
#                        covariates, the AVONET traits and the morphometric PCA axes
#                        (pca1 ... pca9)
#   - trait_pca          the variance on each of those axes
# Run this once, and again whenever preprocess.jl has been re-run.
#
#     julia --project=. create_objects.jl

using CSV
using DataFrames
using JLD2
using Phylo
using SpatialEcology

const INDIR = "data/clean"
const OBJECTFILE = "data/objects.jld2"

# Read a cleaned input with plain `String` columns rather than CSV's own string type, so
# that loading the cache does not need CSV
read_clean(file) = CSV.read(joinpath(INDIR, file), DataFrame; stringtype=String)

"""
    string_sites!(df)

Convert the `site` column of `df` to strings, so site ids read by CSV stay strings.
Returns `df`.
"""
function string_sites!(df)
    df.site = string.(df.site)
    return df
end

# The assemblage of one space ("e" or "g"), with the site covariates keyed by `siteid`
# and the traits keyed by species
function assemblage(suffix, siteid, avonet)
    phylocom = string_sites!(read_clean("phylocom_$suffix.csv"))
    coords = string_sites!(read_clean("coords_$suffix.csv"))
    sitestats = read_clean("sitestats_$suffix.csv")
    sitestats[!, siteid] = string.(sitestats[!, siteid])
    # coordinates were pre-aligned to each phylocom's site order in preprocessing, so
    # they slot straight into the Assemblage (SpatialEcology aligns coords by row order).
    birds = Assemblage(phylocom, coords)
    addsitestats!(birds, sitestats, siteid)
    addtraits!(birds, avonet, :species)
    return birds
end

tree = sort!(parsenewick(read(joinpath(INDIR, "tree.nwk"), String)))

# The raw AVONET traits with the supplied morphometric PCA joined on, so the assemblages
# carry both. The trait scripts read the axes from here instead of fitting a PCA of their
# own; `trait_pca` carries the variance on each axis.
avonet = read_clean("traits.csv")
trait_pca = read_clean("trait_pca_variance.csv")
let scores = read_clean("trait_pca_scores.csv")
    idx = indexin(avonet.species, scores.species)
    any(isnothing, idx) &&
        error("a species in traits.csv is missing from trait_pca_scores.csv")
    rows = Int.(idx)
    for ax in trait_pca.axis
        avonet[!, ax] = scores[rows, ax]
    end
end

birds_e = assemblage("e", :ID_env, avonet)  # PC bins, area, occupancy, ...
birds_g = assemblage("g", :ID_geo, avonet)  # CHELSA bioclim, PC1-3, area, ...

jldsave(OBJECTFILE; tree, birds_e, birds_g, trait_pca)
