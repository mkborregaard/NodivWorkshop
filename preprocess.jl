# Preprocessing for the Nodiv bird analysis. Reads the matched raw data from
# data/data_birds_matched_simplified.rds (an R list: phylogeny, e_space, g_space,
# traits) and writes cleaned inputs to data/clean/ that `script.jl` then loads:
#   - tree.nwk                                the pruned phylogeny
#   - phylocom_e/g.csv, coords_e/g.csv,       per-space occurrences, coordinates,
#     sitestats_e/g.csv                       and site covariates (e = env, g = geo)
#   - traits.csv                              AVONET traits, one row per tree tip
# Run this once (or whenever the raw data changes); it is the slow I/O step.
# Needs R with the sf and ape packages.

using CSV, DataFrames, Phylo, RCall

rdsfile = "data/data_birds_matched_simplified.rds"
outdir = "data/clean"

# The species names are already matched across tree, presences and traits in the
# RDS; only swap spaces for underscores, as Newick tip labels need.
underscore(s) = replace(string(s), " " => "_")

# long-format presence/absence table [site, abundance, species]
make_phylocom(sitevals, species) =
    DataFrame(site = string.(sitevals), abundance = 1, species = underscore.(species))

# reorder a per-site (site, x, y) lookup to the assemblage's site order (unique
# appearance in the phylocom), since SpatialEcology aligns coords by row order.
function align_coords(phylo, lookup)
    sites = unique(phylo.site)
    idx = indexin(sites, string.(lookup.site))
    DataFrame(site = sites, x = lookup.x[idx], y = lookup.y[idx])
end

### Read the RDS in R ----------------------------------------------------------
# The geographic grid is a Behrmann equal-area grid stored as (clipped, simplified)
# lon/lat polygons. Projected back to Behrmann (ESRI:54017) it is exactly regular:
# square cells one degree of longitude wide, with the lattice anchored at 0 (the
# equator is a cell edge). So each cell's row and column follow from where its
# centroid falls - a clipped coastal cell's centroid still lies inside the cell.
R"""
suppressMessages({library(sf); library(ape)})
x <- readRDS($rdsfile)
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
tr$node.label <- NULL                    # support values; Phylo rejects them as duplicate node names
tr$tip.label <- gsub(" ", "_", tr$tip.label)
"""
cellsize = rcopy(R"cellsize")

### Environmental space --------------------------------------------------------
env = rcopy(DataFrame, R"edf")
pres_e = rcopy(DataFrame, R"x$e_space$presence")
phylocom_e = make_phylocom(pres_e.ID_env, pres_e.Species)
sitestats_e = rename(env, :mn_g_d_ => :mean_geo_dist_km, :occupid => :occupied)
sitestats_e.occupied = sitestats_e.occupied .== 1

### Geographic space -----------------------------------------------------------
geo = rcopy(DataFrame, R"gdf")
pres_g = rcopy(DataFrame, R"x$g_space$presence")
phylocom_g = make_phylocom(pres_g.ID_geo, pres_g.Species)
allunique(zip(geo.col, geo.row)) || error("two grid cells fall in the same Behrmann cell")
sitestats_g = select(geo, Not([:col, :row]))
sitestats_g.area_m = parse.(Float64, sitestats_g.area_m)
sitestats_g.ID_env = [s == "NA" ? missing : s for s in sitestats_g.ID_env]

### Phylogeny: keep only the taxa present in both spaces, which also drops the few
### species found in geographic space alone.
tree = rcopy(RootedTree, R"tr")
shared = intersect(getleafnames(tree), unique(phylocom_e.species), unique(phylocom_g.species))
keeptips!(tree, shared)
sort!(tree)
filter!(r -> r.species in shared, phylocom_e)
filter!(r -> r.species in shared, phylocom_g)

# coordinates: PC bin midpoints for env; Behrmann cell centres in km for geo
coords_e = align_coords(phylocom_e,
    DataFrame(site = env.ID_env, x = env.pc1_mid, y = env.pc2_mid))
coords_g = align_coords(phylocom_g,
    DataFrame(site = geo.ID_geo, x = (geo.col .+ 0.5) .* cellsize ./ 1000,
                                 y = (geo.row .+ 0.5) .* cellsize ./ 1000))

### Traits: one row per tree tip, keyed by `species` for addtraits!
traits = rcopy(DataFrame, R"x$traits")
traits = select(traits, :Species1 => ByRow(underscore) => :species, Not(:Species1))
filter!(r -> r.species in shared, traits)

### Write the cleaned inputs ---------------------------------------------------
mkpath(outdir)
CSV.write(joinpath(outdir, "phylocom_e.csv"), phylocom_e)
CSV.write(joinpath(outdir, "coords_e.csv"), coords_e)
CSV.write(joinpath(outdir, "sitestats_e.csv"), sitestats_e)
CSV.write(joinpath(outdir, "phylocom_g.csv"), phylocom_g)
CSV.write(joinpath(outdir, "coords_g.csv"), coords_g)
CSV.write(joinpath(outdir, "sitestats_g.csv"), sitestats_g)
CSV.write(joinpath(outdir, "traits.csv"), traits)
# write the tree as Newick via R's ape::write.tree (Phylo has no Newick writer).
# ape drops internal node labels, so re-reading renumbers internal nodes - fine
# here, the labels are just auto-generated "Node N" placeholders anyway.
treefile = joinpath(outdir, "tree.nwk")
R"write.tree($tree, file = $treefile)"
