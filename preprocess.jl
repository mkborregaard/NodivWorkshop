# Preprocessing for the Nodiv bird analysis. Reads the matched raw data from
# data/data_birds_matched_simplified.rds (an R list: phylogeny, e_space, g_space,
# traits) and writes cleaned inputs to data/clean/ that `script.jl` then loads:
#   - tree.nwk                                the pruned phylogeny, taxon-named nodes
#   - phylocom_e/g.csv, coords_e/g.csv,       per-space occurrences, coordinates,
#     sitestats_e/g.csv                       and site covariates (e = env, g = geo)
#   - traits.csv                              AVONET traits, one row per tree tip
# Run this once (or whenever the raw data changes); it is the slow I/O step.
# Needs R with the sf and ape packages.

using CSV
using DataFrames
using Phylo
using RCall

const RDSFILE = "data/data_birds_matched_simplified.rds"
const OUTDIR = "data/clean"

# The species names are already matched across tree, presences and traits in the
# RDS; only swap spaces for underscores, as Newick tip labels need.
underscore(s) = replace(string(s), " " => "_")

# Long-format presence/absence table [site, abundance, species]
function make_phylocom(sitevals, species)
    return DataFrame(; site=string.(sitevals), abundance=1, species=underscore.(species))
end

# reorder a per-site (site, x, y) lookup to the assemblage's site order (unique
# appearance in the phylocom), since SpatialEcology aligns coords by row order.
function align_coords(phylo, lookup)
    sites = unique(phylo.site)
    idx = indexin(sites, string.(lookup.site))
    return DataFrame(; site=sites, x=lookup.x[idx], y=lookup.y[idx])
end

### ---- Read the RDS in R ---- ###
# The geographic grid is a Behrmann equal-area grid stored as (clipped, simplified)
# lon/lat polygons. Projected back to Behrmann (ESRI:54017) it is exactly regular:
# square cells one degree of longitude wide, with the lattice anchored at 0 (the
# equator is a cell edge). So each cell's row and column follow from where its
# centroid falls - a clipped coastal cell's centroid still lies inside the cell.
R"""
suppressMessages({library(sf); library(ape)})
x <- readRDS($RDSFILE)
g <- x$g_space$grid_sf
gb <- st_transform(st_geometry(g), "ESRI:54017")
cellsize <- diff(sf_project("EPSG:4326", "ESRI:54017", rbind(c(0, 0), c(1, 0)))[, 1])
cb <- st_coordinates(st_centroid(gb))
ll <- st_coordinates(st_transform(st_centroid(gb), "EPSG:4326"))
gdf <- st_drop_geometry(g)
gdf$col <- floor(cb[, 1] / cellsize)
gdf$row <- floor(cb[, 2] / cellsize)
gdf$lon <- ll[, 1]
gdf$lat <- ll[, 2]
edf <- st_drop_geometry(x$e_space$grid_sf)
tr <- x$phylogeny
# Support values; Phylo rejects them as duplicate node names
tr$node.label <- NULL
tr$tip.label <- gsub(" ", "_", tr$tip.label)
"""
cellsize = rcopy(R"cellsize")

### ---- Environmental space ---- ###
env = rcopy(DataFrame, R"edf")
pres_e = rcopy(DataFrame, R"x$e_space$presence")
phylocom_e = make_phylocom(pres_e.ID_env, pres_e.Species)
sitestats_e = rename(env, :mn_g_d_ => :mean_geo_dist_km, :occupid => :occupied)
sitestats_e.occupied = sitestats_e.occupied .== 1

### ---- Geographic space ---- ###
geo = rcopy(DataFrame, R"gdf")
pres_g = rcopy(DataFrame, R"x$g_space$presence")
phylocom_g = make_phylocom(pres_g.ID_geo, pres_g.Species)
allunique(zip(geo.col, geo.row)) || error("two grid cells fall in the same Behrmann cell")
sitestats_g = select(geo, Not([:col, :row]))
sitestats_g.area_m = parse.(Float64, sitestats_g.area_m)
sitestats_g.ID_env = [s == "NA" ? missing : s for s in sitestats_g.ID_env]

### ---- Phylogeny ---- ###
# Keep only the taxa present in both spaces, which also drops the few species found in
# geographic space alone.
tree = rcopy(RootedTree, R"tr")
shared = intersect(
    getleafnames(tree), unique(phylocom_e.species), unique(phylocom_g.species)
)
keeptips!(tree, shared)
sort!(tree)
sharedset = Set(shared)
filter!(r -> r.species in sharedset, phylocom_e)
filter!(r -> r.species in sharedset, phylocom_g)

# Coordinates: PC bin midpoints for env; Behrmann cell centres in km for geo
coords_e = align_coords(
    phylocom_e, DataFrame(; site=env.ID_env, x=env.pc1_mid, y=env.pc2_mid)
)
coords_g = align_coords(
    phylocom_g,
    DataFrame(;
        site=geo.ID_geo,
        x=(geo.col .+ 0.5) .* cellsize ./ 1000,
        y=(geo.row .+ 0.5) .* cellsize ./ 1000,
    ),
)

### ---- Traits ---- ###
# One row per tree tip, keyed by `species` for addtraits!
avonet = rcopy(DataFrame, R"x$traits")
avonet = select(avonet, :Species1 => ByRow(underscore) => :species, Not(:Species1))
filter!(r -> r.species in sharedset, avonet)

### ---- Name the internal nodes that are exactly a genus, family or order ---- ###
# Round-trip the tree through ape's Newick first: that is the file script.jl used to
# read, so parsing it back gives the same auto-generated "Node N" names that the
# cached node analysis (data/node_analysis.jld2) is keyed by.
tree = parsenewick(rcopy(String, R"write.tree($tree)"))

# For each genus (from the species name), family and order with more than one
# species, rename its MRCA to the taxon if the taxon is monophyletic - its species
# are exactly the tips below that node. Non-monophyletic taxa stay unnamed. Where one
# clade is several taxa at once (e.g. a family of a single genus) the highest rank wins.
function taxonnodes(tree, avonet)
    genus = String.(first.(split.(avonet.species, "_")))
    taxonnames = Dict{String,String}()
    # Low to high rank: higher overwrites
    for taxa in (genus, avonet.Family1, avonet.Order1)
        for taxon in unique(taxa)
            sp = avonet.species[taxa .== taxon]
            length(sp) > 1 || continue
            node = getnodename(tree, mrca(tree, sp))
            ntips = count(n -> isleaf(tree, n), getdescendants(tree, node))
            if ntips == length(sp)
                taxonnames[node] = taxon
            end
        end
    end
    return taxonnames
end
for (node, taxon) in taxonnodes(tree, avonet)
    renamenode!(tree, node, taxon)
end
# Ladderize (order each node's clades by size) for plotting. parsenewick does not keep
# the file's child order, so script.jl ladderizes again after reading the tree.
sort!(tree)

### ---- Write the cleaned inputs ---- ###
mkpath(OUTDIR)
CSV.write(joinpath(OUTDIR, "phylocom_e.csv"), phylocom_e)
CSV.write(joinpath(OUTDIR, "coords_e.csv"), coords_e)
CSV.write(joinpath(OUTDIR, "sitestats_e.csv"), sitestats_e)
CSV.write(joinpath(OUTDIR, "phylocom_g.csv"), phylocom_g)
CSV.write(joinpath(OUTDIR, "coords_g.csv"), coords_g)
CSV.write(joinpath(OUTDIR, "sitestats_g.csv"), sitestats_g)
CSV.write(joinpath(OUTDIR, "traits.csv"), avonet)
# Phylo's own Newick writer keeps every internal node name, taxon or "Node N", so
# re-reading the file gives back exactly these names.
write(joinpath(OUTDIR, "tree.nwk"), tree)
