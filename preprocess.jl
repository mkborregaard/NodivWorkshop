# Preprocessing for the Nodiv bird analysis. Reads the raw data, harmonises the
# taxonomy, and writes cleaned inputs to data/clean/ that `script.jl` then loads:
#   - tree.jld2                               the pruned phylogeny
#   - phylocom_e/g.csv, coords_e/g.csv,       per-space occurrences, coordinates,
#     sitestats_e/g.csv                       and site covariates (e = env, g = geo)
# Run this once (or whenever the raw data changes); it is the slow I/O step.

using CSV, DataFrames, Shapefile, Phylo, JLD2

# Cross.csv maps the PAM/BirdLife taxonomy (Species1) to the tree/BirdTree
# taxonomy (Species3); relabel PAM species to the tree names so the datasets
# share as many taxa as possible.
cross = CSV.read("data/Cross.csv", DataFrame)
namemap = Dict(string(r.Species1) => replace(string(r.Species3), " " => "_")
               for r in eachrow(cross) if !ismissing(r.Species1) && !ismissing(r.Species3))

# long-format presence/absence table [site, abundance, species], species relabelled
make_phylocom(sitevals, species) =
    DataFrame(site = string.(sitevals), abundance = 1,
              species = [get(namemap, s, replace(s, " " => "_")) for s in species])

# reorder a per-site (site, x, y) lookup to the assemblage's site order (unique
# appearance in the phylocom), since SpatialEcology aligns coords by row order.
function align_coords(phylo, lookup)
    sites = unique(phylo.site)
    idx = indexin(sites, string.(lookup.site))
    DataFrame(site = sites, x = lookup.x[idx], y = lookup.y[idx])
end

### Environmental space --------------------------------------------------------
env = CSV.read("data/Env.csv", DataFrame)
pam_e = CSV.read("data/PAM_E.csv", DataFrame)
phylocom_e = make_phylocom(pam_e.ID_env, pam_e.Species)
coords_e = align_coords(phylocom_e,                       # PC bin midpoints
    DataFrame(site = string.(env.ID_env), x = (env.xmin .+ env.xmax) ./ 2, y = (env.ymin .+ env.ymax) ./ 2))
sitestats_e = env

### Geographic space -----------------------------------------------------------
pam_g = CSV.read("data/PAM_G.csv", DataFrame)
phylocom_g = make_phylocom(pam_g.ID_geo, pam_g.Species)

shp = Shapefile.Table(joinpath("data", "g_space", "BehrmannMeterGrid_WGS84_land_PCA_30.shp"))
centroid(g) = (ex = extrema(p.x for p in g.points); ey = extrema(p.y for p in g.points);
               ((ex[1] + ex[2]) / 2, (ey[1] + ey[2]) / 2))
cents = centroid.(Shapefile.shapes(shp))
# Behrmann equal-area cells -> regular grid: rank longitude into columns and
# sin(latitude) into rows, scaling the sin axis by 1/(dlon*cos^2(30deg)) ~ 76.4 so
# cells come out square; the -0.5 centres bins on the bands (avoids an equator gap).
gridindex(v) = (u = sort(unique(v)); pos = Dict(u .=> eachindex(u)); Float64[pos[x] for x in v])
behrmann = 1 / (deg2rad(1) * cosd(30)^2)
coords_g = align_coords(phylocom_g,
    DataFrame(site = string.(shp.ID_geo),
              x = gridindex(round.(Int, first.(cents))),
              y = gridindex(round.(Int, sin.(deg2rad.(last.(cents))) .* behrmann .- 0.5))))
sitestats_g = select(DataFrame(shp), Not(:geometry))      # CHELSA bioclim, PC1-3, area, ...
sitestats_g.ID_geo = string.(sitestats_g.ID_geo)

### Phylogeny: strip numeric internal (support) labels Phylo rejects as duplicate
### node names, then keep only the taxa shared by both spaces.
tree = parsenewick(replace(read("data/birds.nwk", String), r"\)[0-9.]+" => ")"))
keeptips!(tree, intersect(getleafnames(tree), unique(phylocom_e.species), unique(phylocom_g.species)))
sort!(tree)

### Write the cleaned inputs ---------------------------------------------------
mkpath("data/clean")
CSV.write("data/clean/phylocom_e.csv", phylocom_e)
CSV.write("data/clean/coords_e.csv", coords_e)
CSV.write("data/clean/sitestats_e.csv", sitestats_e)
CSV.write("data/clean/phylocom_g.csv", phylocom_g)
CSV.write("data/clean/coords_g.csv", coords_g)
CSV.write("data/clean/sitestats_g.csv", sitestats_g)
jldsave("data/clean/tree.jld2"; tree)
