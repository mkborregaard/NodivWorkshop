# Preprocessing for the Nodiv bird analysis. Reads the matched raw data from
# data/Bird_data_PCA30.rds (an R list: g_space, e_space, phylogeny, traits,
# env_PCA_loadings) and writes cleaned inputs to data/clean/ that `create_objects.jl`
# then loads:
#   - tree.nwk                                the pruned phylogeny, taxon-named nodes
#   - phylocom_e/g.csv, coords_e/g.csv,       per-space occurrences, coordinates,
#     sitestats_e/g.csv                       and site covariates (e = env, g = geo)
#   - traits.csv                              AVONET traits, one row per tree tip
#   - trait_pca_scores.csv,                   the supplied morphometric PCA: scores per
#     trait_pca_loadings.csv,                 species, loadings, and variance per axis
#     trait_pca_variance.csv
#   - env_pca_loadings.csv                    loadings of the environmental PCA
#   - invented_nodes.csv                      the internal nodes that exist only because a
#                                             polytomy was resolved, so a divergent node
#                                             can be checked against them
# Run this once (or whenever the raw data changes); it is the slow I/O step.
# Needs R with the sf and ape packages.
#
# See docs/data_migration.md for how this file differs from the earlier
# data_birds_matched_simplified.rds and why each choice below was made.

using CSV
using DataFrames
using Phylo
using RCall
using Statistics

const RDSFILE = "data/Bird_data_PCA30.rds"
# The raw CHELSA / MODIS / GLWD / elevation layers per geographic cell exist only in the
# earlier file; the new one keeps just their PCA. The cells are identical (same ids, same
# order, centroids agreeing to 5e-13 degrees), so they are carried across by id below.
# Optional: without it, sitestats_g.csv holds only the new file's columns.
const RAWLAYERFILE = "data/data_birds_matched_simplified.rds"
const OUTDIR = "data/clean"

# Which of the three supplied phylogenies to use, and the seed for resolving its
# polytomies. McTavish_2025 has 49 polytomies, so resolving them stochastically invents
# 63 of 10,353 nodes (0.6%) and leaves ~3 of ~1,914 analysable nodes arbitrary; the
# alternatives cost far more (Claramund_2025: 1,684 invented nodes, 9% arbitrary) or are
# worse resolved (Jetz_2012, the tree used before this migration). The resolved tree is
# written to tree.nwk, so results are reproducible from that file whatever the seed.
const TREE = "McTavish_2025"
const TREESEED = 20261124

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
g <- x$g_space$g_space_sf
gb <- st_transform(st_geometry(g), "ESRI:54017")
cellsize <- diff(sf_project("EPSG:4326", "ESRI:54017", rbind(c(0, 0), c(1, 0)))[, 1])
cb <- st_coordinates(st_centroid(gb))
ll <- st_coordinates(st_transform(st_centroid(gb), "EPSG:4326"))
gdf <- st_drop_geometry(g)
gdf$col <- floor(cb[, 1] / cellsize)
gdf$row <- floor(cb[, 2] / cellsize)
gdf$lon <- ll[, 1]
gdf$lat <- ll[, 2]
edf <- st_drop_geometry(x$e_space$e_space_sf)
pres_e <- x$e_space$e_space_PAM
pres_g <- x$g_space$g_space_PAM
trdf <- x$traits$avonet_raw_traits
scores <- x$traits$PCA_scores
score_species <- rownames(x$traits$PCA_scores)
tld <- x$traits$PCA_loadings
eld <- x$env_PCA_loadings
"""
cellsize = rcopy(R"cellsize")

### ---- Phylogeny: select, resolve, make ultrametric ---- ###
# The supplied trees carry no node labels and already use underscores, but strip and
# swap defensively. `multi2di` resolves each polytomy into a random ladder of
# zero-length branches: Nodiv only analyses a node with exactly two children
# (`Nodiv._analysis_prepass`), so an unresolved node is dropped from the analysis
# entirely. The tree is then made ultrametric by extending terminal branches, which on
# this tree closes a tip-depth spread of ~1e-5 on a root depth of ~117 (8.6e-8
# relative) - a rounding artifact, not a real deviation. Both steps are asserted below.
R"""
tr <- x$phylogeny[[$TREE]]
tr$node.label <- NULL
tr$tip.label <- gsub(" ", "_", tr$tip.label)
polytomies_before <- sum(table(tr$edge[, 1]) > 2)
nodes_before <- Nnode(tr)
tr_unresolved <- tr          # kept to tell the invented nodes apart further down
set.seed($TREESEED)
tr <- multi2di(tr, random = TRUE)
d <- node.depth.edgelength(tr)[seq_len(Ntip(tr))]
spread_before <- max(d) - min(d)
term <- match(seq_len(Ntip(tr)), tr$edge[, 2])
tr$edge.length[term] <- tr$edge.length[term] + (max(d) - d)
d2 <- node.depth.edgelength(tr)[seq_len(Ntip(tr))]
spread_after <- max(d2) - min(d2)
"""
@info(
    "phylogeny $TREE",
    tips = rcopy(Int, R"Ntip(tr)"),
    polytomies_resolved = rcopy(Int, R"polytomies_before"),
    nodes_invented = rcopy(Int, R"Nnode(tr) - nodes_before"),
    root_depth = rcopy(R"max(d2)"),
    tip_depth_spread_before = rcopy(R"spread_before"),
    tip_depth_spread_after = rcopy(R"spread_after"),
)
rcopy(Bool, R"is.ultrametric(tr)") || error("the tree is not ultrametric after extension")
rcopy(Bool, R"is.binary(tr)") || error("the tree still has polytomies after multi2di")

### ---- Environmental space ---- ###
# Already named as the cleaned inputs want (`mean_geo_dist_km`, `occupied::Bool`), and
# `n_cells` and a fully populated `area_km` are new here. The bin edges xmin/xmax/
# ymin/ymax of the earlier file are gone; nothing downstream read them.
env = rcopy(DataFrame, R"edf")
pres_e = rcopy(DataFrame, R"pres_e")
phylocom_e = make_phylocom(pres_e.ID_env, pres_e.Species)
sitestats_e = env

### ---- Geographic space ---- ###
geo = rcopy(DataFrame, R"gdf")
pres_g = rcopy(DataFrame, R"pres_g")
# 218 of the PAM's site ids are absent from the grid, so they have no coordinates. They
# carry ~0.15% of the records and dropping them loses no species outright.
gridsites = Set(string.(geo.ID_geo))
unmapped = Set(setdiff(unique(string.(pres_g.ID_geo)), gridsites))
dropped = count(s -> string(s) in unmapped, pres_g.ID_geo)
filter!(r -> !(string(r.ID_geo) in unmapped), pres_g)
@info(
    "geographic presences",
    unmapped_sites = length(unmapped),
    dropped_records = dropped,
    kept = nrow(pres_g),
)
phylocom_g = make_phylocom(pres_g.ID_geo, pres_g.Species)
allunique(zip(geo.col, geo.row)) || error("two grid cells fall in the same Behrmann cell")
sitestats_g = select(geo, Not([:col, :row]))
sitestats_g.ID_env = [
    ismissing(s) || s == "NA" ? missing : String(s) for s in sitestats_g.ID_env
]

# Carry the raw environmental layers over from the earlier file, matched by cell id.
if isfile(RAWLAYERFILE)
    R"""
    rawg <- st_drop_geometry(readRDS($RAWLAYERFILE)$g_space$grid_sf)
    """
    rawg = rcopy(DataFrame, R"rawg")
    rawg.area_m = parse.(Float64, rawg.area_m)
    extracols = setdiff(names(rawg), names(sitestats_g))
    idx = indexin(string.(sitestats_g.ID_geo), string.(rawg.ID_geo))
    any(isnothing, idx) && error("a grid cell is missing from $RAWLAYERFILE")
    for n in extracols
        sitestats_g[!, n] = rawg[idx, n]
    end
    @info "raw environmental layers carried over" count = length(extracols)
else
    @warn "no $RAWLAYERFILE: sitestats_g.csv will hold only the new file's columns"
end

### ---- Phylogeny: prune to the taxa the data share ---- ###
tree = rcopy(RootedTree, R"tr")
shared = intersect(
    getleafnames(tree), unique(phylocom_e.species), unique(phylocom_g.species)
)
keeptips!(tree, shared)
sort!(tree)
sharedset = Set(shared)
filter!(r -> r.species in sharedset, phylocom_e)
filter!(r -> r.species in sharedset, phylocom_g)
@info(
    "species",
    tree_tips = length(shared),
    env_records = nrow(phylocom_e),
    geo_records = nrow(phylocom_g),
)

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
avonet = rcopy(DataFrame, R"trdf")
avonet = select(avonet, :Species1 => ByRow(underscore) => :species, Not(:Species1))
filter!(r -> r.species in sharedset, avonet)

# The supplied morphometric PCA, used in place of one recomputed in the trait scripts.
# It is fit on nine non-redundant traits (the earlier in-workshop PCA spanned eleven
# columns, including Kipps_Distance and the hand-wing index, which are derived from the
# others). All nine axes are kept, so the proportions below are exact.
trait_scores = DataFrame(rcopy(Matrix{Float64}, R"scores"), ["pca$i" for i in 1:9])
trait_scores.species = underscore.(rcopy(Vector{String}, R"score_species"))
filter!(r -> r.species in sharedset, trait_scores)
variances = [var(trait_scores[!, "pca$i"]) for i in 1:9]
trait_variance = DataFrame(;
    axis=["pca$i" for i in 1:9],
    variance=variances,
    proportion=variances ./ sum(variances),
    cumulative=cumsum(variances) ./ sum(variances),
)
trait_loadings = rename(rcopy(DataFrame, R"tld"), :X => :trait)
# Spell the traits as traits.csv's columns do, so the two can be joined
trait_loadings.trait = replace.(trait_loadings.trait, "." => "_")
env_loadings = rename(rcopy(DataFrame, R"eld"), :X => :variable)

### ---- Name the internal nodes that are exactly a genus, family or order ---- ###
# Round-trip the tree through ape's Newick first, so the auto-generated "Node N" names
# are the ones a reader of the file sees.
tree = parsenewick(rcopy(String, R"write.tree($tree)"))

# For each genus (from the species name), family and order with more than one
# species, rename its MRCA to the taxon if the taxon is monophyletic - its species
# are exactly the tips below that node. Non-monophyletic taxa stay unnamed. Where one
# clade is several taxa at once (e.g. a family of a single genus) the highest rank wins.
function taxon_nodes(tree, avonet)
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
named = taxon_nodes(tree, avonet)
for (node, taxon) in named
    renamenode!(tree, node, taxon)
end
@info "named internal nodes" count = length(named)
# Ladderize (order each node's clades by size) for plotting. parsenewick does not keep
# the file's child order, so script.jl ladderizes again after reading the tree.
sort!(tree)

### ---- Record the nodes that resolution invented ---- ###
# `multi2di` split each polytomy into a random ladder, so some of the tree's clades are an
# artifact of that draw rather than anything the source tree asserts. A node is invented
# if its clade is not a clade of the unresolved tree. Written out so that a divergent node
# can be checked against the list: a high rms on an invented node is not evidence about
# the real phylogeny. On the McTavish_2025 draw used here only one of the 63 is analysable
# at all, and it is well below the divergence threshold in both spaces.
R"""
tr0 <- keep.tip(tr_unresolved, $shared)
p0 <- prop.part(tr0)
lab0 <- attr(p0, "labels")
source_clades <- vapply(p0, function(i) paste(sort(lab0[i]), collapse = "\u0001"), character(1))
"""
source_clades = Set(rcopy(Vector{String}, R"source_clades"))
# getdescendants of a leaf is empty, so a leaf counts as the one species it is
cladesize(node) =
    isleaf(tree, node) ? 1 : count(n -> isleaf(tree, n), getdescendants(tree, node))
clade_signature(node) =
    join(sort!([getnodename(tree, d) for d in getdescendants(tree, node) if isleaf(tree, d)]), '\x01')

invented = DataFrame(; node=String[], species=Int[], child1=Int[], child2=Int[])
for node in getnodenames(tree)
    isleaf(tree, node) && continue
    clade_signature(node) in source_clades && continue
    cs = [cladesize(getnodename(tree, c)) for c in getchildren(tree, node)]
    push!(invented, (node, cladesize(node), get(cs, 1, 0), get(cs, 2, 0)))
end
sort!(invented, :node)
@info "nodes invented by resolving polytomies" count = nrow(invented)
nrow(invented) == rcopy(Int, R"Nnode(tr) - nodes_before") ||
    error("the invented-node count disagrees with what multi2di added")

### ---- Write the cleaned inputs ---- ###
mkpath(OUTDIR)
CSV.write(joinpath(OUTDIR, "phylocom_e.csv"), phylocom_e)
CSV.write(joinpath(OUTDIR, "coords_e.csv"), coords_e)
CSV.write(joinpath(OUTDIR, "sitestats_e.csv"), sitestats_e)
CSV.write(joinpath(OUTDIR, "phylocom_g.csv"), phylocom_g)
CSV.write(joinpath(OUTDIR, "coords_g.csv"), coords_g)
CSV.write(joinpath(OUTDIR, "sitestats_g.csv"), sitestats_g)
CSV.write(joinpath(OUTDIR, "traits.csv"), avonet)
CSV.write(joinpath(OUTDIR, "trait_pca_scores.csv"), trait_scores)
CSV.write(joinpath(OUTDIR, "trait_pca_variance.csv"), trait_variance)
CSV.write(joinpath(OUTDIR, "trait_pca_loadings.csv"), trait_loadings)
CSV.write(joinpath(OUTDIR, "env_pca_loadings.csv"), env_loadings)
CSV.write(joinpath(OUTDIR, "invented_nodes.csv"), invented)
# Phylo's own Newick writer keeps every internal node name, taxon or "Node N", so
# re-reading the file gives back exactly these names.
write(joinpath(OUTDIR, "tree.nwk"), tree)
