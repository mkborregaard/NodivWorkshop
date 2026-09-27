# Build the analysis objects from the cleaned inputs in data/clean/ (written by
# preprocess.jl) and cache them in data/objects.jld2, which script.jl, traitsscript.jl
# and explorer.jl load:
#   - tree               the phylogeny, ladderized
#   - birds_e, birds_g   the assemblages of each space (e = env, g = geo), with their site
#                        covariates and the AVONET traits
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
avonet = read_clean("traits.csv")
birds_e = assemblage("e", :ID_env, avonet)  # PC bins, area, occupancy, ...
birds_g = assemblage("g", :ID_geo, avonet)  # CHELSA bioclim, PC1-3, area, ...

jldsave(OBJECTFILE; tree, birds_e, birds_g)
