# Node-based analysis of bird diversity, in environmental (birds_e) and geographic

# (birds_g) space. Run preprocess.jl first to build the cleaned inputs in
# data/clean/; this script loads them, builds the assemblages, computes and caches
# the node analysis, and explores the results.

using CSV, DataFrames, SpatialEcology, Phylo, Plots
using MultivariateStats, Statistics, JLD2, LogExpFunctions, GLM
using Nodiv

default(color = cgrad(:Spectral, rev = true))

### Load the cleaned inputs (from preprocess.jl) and build the assemblages -----

tree = parsenewick(read("data/clean/tree.nwk", String))
strsite!(df) = (df.site = string.(df.site); df)   # site ids stay strings after CSV
phylocom_e  = strsite!(CSV.read("data/clean/phylocom_e.csv", DataFrame))
coords_e    = strsite!(CSV.read("data/clean/coords_e.csv", DataFrame))
sitestats_e = CSV.read("data/clean/sitestats_e.csv", DataFrame)
phylocom_g  = strsite!(CSV.read("data/clean/phylocom_g.csv", DataFrame))
coords_g    = strsite!(CSV.read("data/clean/coords_g.csv", DataFrame))
sitestats_g = CSV.read("data/clean/sitestats_g.csv", DataFrame)
sitestats_g.ID_geo = string.(sitestats_g.ID_geo)

# coordinates were pre-aligned to each phylocom's site order in preprocessing, so
# they slot straight into the Assemblage (SpatialEcology aligns coords by row order).
birds_e = Assemblage(phylocom_e, coords_e)
addsitestats!(birds_e, sitestats_e, :ID_env)   # PC bins, area, occupancy, ...
plot(birds_e)

birds_g = Assemblage(phylocom_g, coords_g)
addsitestats!(birds_g, sitestats_g, :ID_geo)   # CHELSA bioclim, PC1-3, area, ...
plot(birds_g)


### Heavy step: GND + SOS for every node, both spaces, cached to disk ----------
# `node_analysis` computes GND and the per-cell SOS together (SOS is needed for
# GND anyway). This is the slow part - randomisations over the whole tree, and
# ~18k cells for the geographic scan - so cache it: re-running the script just
# reloads the results and jumps straight to the plotting below.
cachefile = "data/node_analysis.jld2"
if !isfile(cachefile)
    res_e = node_metrics(birds_e, tree; nsims = 200)
    res_g = node_metrics(birds_g, tree; nsims = 200)
    jldsave(cachefile; res_e, res_g)
end
res_e, res_g = load(cachefile, "res_e", "res_g")   # each a NodeMetrics (gnd/rms/spatial/ses/pval + sos)

### ---- Exploratory plotting (from the cached NodeMetrics; `_e` vs `_g`) ----- ###

# strongly divergent nodes in each space. by = :gnd keeps the original GND > 0.8
# selection; switch to the default (RMS-SOS > 1.5) or by = :pval to use the
# effect-size / null-calibrated scores instead.
divergent_e = divergent_nodes(res_e; by = :gnd, threshold = 0.8)
divergent_g = divergent_nodes(res_g; by = :gnd, threshold = 0.8)
divergent = divergent_e ∩ divergent_g

# GND of just the divergent nodes mapped onto the tree (plot_gnd marks every node
# in the Dict it is given, so pass the divergent subset rather than the full result)
plot_gnd(tree, Dict(n => res_e.gnd[n] for n in divergent_e))
plot_gnd(tree, Dict(n => res_g.gnd[n] for n in divergent_g))

# SOS of the most divergent node mapped onto each space (cached SOS, no recompute)
focal_e = argmax(n -> res_e.gnd[n], divergent_e)
plot(res_e.sos[focal_e], birds_e, fillcolor = :RdYlBu, clim = (-8, 8), title = "env SOS - $focal_e")
focal_g = argmax(n -> res_g.gnd[n], divergent_g)
plot(res_g.sos[focal_g], birds_g, fillcolor = :RdYlBu, clim = (-8, 8), title = "geo SOS - $focal_g")

# parent/SOS/children panel for that node (4th arg = cached SOS, no recompute)
plot_node(birds_e, tree, focal_e, res_e)
plot_node(birds_g, tree, focal_g, res_g)

# ordinate the divergent nodes by SOS-pattern similarity (cached SOS -> distances
# from Nodiv -> MDS; presentation stays here)
function sos_mds_plot(res, nodes, title)
    coords = predict(fit(MDS, sos_distances(res, nodes); distances = true, maxoutdim = 2))
    scatter(coords[1, :], coords[2, :], label = "",
            series_annotations = text.(nodes, 6, :bottom),
            xlabel = "MDS axis 1", ylabel = "MDS axis 2", title = title)
end
sos_mds_plot(res_e, divergent, "SOS-pattern similarity (environmental)")
sos_mds_plot(res_g, divergent, "SOS-pattern similarity (geographic)")

focal = "Node 10531"
plot_node(birds_e, tree, focal, res_e)
plot_node(birds_g, tree, focal, res_g)

same = divergent_e ∩ divergent_g

nodes = collect(keys(res_e.gnd))
dat = DataFrame(
    :logit_g => [logit(res_g.gnd[n]) for n in nodes],
    :logit_e => [logit(res_e.gnd[n]) for n in nodes]
)
dat = filter(row -> all(x -> !ismissing(x) && isfinite(x), row), dat)

scatter(dat.logit_g, dat.logit_e,
        xlabel = "geo GND", ylabel = "env GND", label = "")

mod = lm(@formula(logit_e ~ logit_g), dat)


sizes = Dict(node => noccupied(get_clade(birds_e, tree, node)) for node in nodes)
histogram(collect(values(sizes)))

plot(tree, treetype = :fan, marker_z = sizes, showtips = false, msw = 0)

scatter([sizes[n] for n in nodes], [res_e.gnd[n] for n in nodes],
        xlabel = "occupied env sites", ylabel = "env GND", label = "")

plot!([-2, 4], [-2, 4], c = :red, label = "")