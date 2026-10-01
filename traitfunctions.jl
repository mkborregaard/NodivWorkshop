# The functions of the trait analysis, shared by traitsscript.jl and traitsexplorer.jl: the
# PCA of the AVONET morphometrics, the trait probability densities (TPD) of clades and
# their overlap, and the trait explorer. Needs DataFrames, GLMakie, MultivariateStats,
# NodivMakie, Phylo, SpatialEcology, Statistics and StatsFuns loaded.

### ---- Traits: PCA of the AVONET morphometrics ---- ###

# Shapiro-Francia W': the squared correlation of the sorted values with the normal quantiles
function normality(x)
    n = length(x)
    q = norminvcdf.(((1:n) .- 0.375) ./ (n + 0.25))
    return cor(sort(x), q)^2
end

# PCA of the columns of `df`: a column is logged where that raises its W' by more than
# `tol`, then all are z-transformed. Returns the PCA, the logged columns and the
# z-transformed matrix (species x traits) it was fit on.
function trait_pca(df; tol=0.01)
    X = Matrix{Float64}(df)
    logged = [all(>(0), x) && normality(log.(x)) - normality(x) > tol for x in eachcol(X)]
    X[:, logged] .= log.(X[:, logged])
    Z = (X .- mean(X; dims=1)) ./ std(X; dims=1)
    return (; pca=fit(PCA, permutedims(Z); pratio=1), logged=names(df)[logged], z=Z)
end

# The PCA of the morphometrics of the species in `asm`: a DataFrame of their scores on the
# first `naxes` axes (columns pca1, pca2, ... and species, for `addtraits!`), the fit from
# `trait_pca`, and the proportion of the variance on each axis
function trait_pcs(asm; naxes=4)
    t = traits(asm)
    pcafit = trait_pca(t[:, Between(:Beak_Length_Culmen, :Mass)])
    explained = principalvars(pcafit.pca) ./ var(pcafit.pca)
    scores = predict(pcafit.pca, permutedims(pcafit.z))[1:naxes, :]
    pcs = DataFrame(permutedims(scores), ["pca$i" for i in 1:naxes])
    pcs.species = t.name
    return pcs, pcafit, explained
end

### ---- Trait probability densities and their overlap ---- ###

# Species => point in trait space, from the columns `x` and `y` of an assemblage's traits
function trait_points(asm, x, y)
    t = traits(asm)
    return Dict(zip(t.name, Point2d.(t[!, x], t[!, y])))
end

# The trait-space points of the species of each of `node`'s two child clades
function child_points(tree, node, pts)
    return [
        [pts[sp] for sp in nodespecies(tree, getnodename(tree, c)) if haskey(pts, sp)] for
        c in getchildren(tree, node)[1:2]
    ]
end

# A grid over the space of `pts` (species => point, or a vector of points): `n` points along
# each axis, spanning the range of the points extended by `extend` of it on each side
function trait_grid(pts; n=200, extend=0.15)
    function axis(v)
        lo, hi = extrema(v)
        r = extend * (hi - lo)
        return range(lo - r, hi + r; length=n)
    end
    return axis(first.(values(pts))), axis(last.(values(pts)))
end

# The density above which the TPD `d` (summing to 1) holds the fraction `p` of its
# probability
function tpd_level(d, p)
    v = sort(vec(d); rev=true)
    return v[min(searchsortedfirst(cumsum(v), p), length(v))]
end

# The trait probability density (TPD, Carmona et al. 2016) of 2-d points on `grid`: a
# Gaussian kernel density with a normal-reference bandwidth in each dimension, evaluated
# at the grid points, cut to the cells holding the fraction `alpha` of the probability
# (as Carmona's TPDs) and scaled to sum to 1; `nothing` for fewer than three distinct points
function tpd(pts, (xs, ys); alpha=0.95)
    length(unique(pts)) < 3 && return nothing
    n = length(pts)
    function kernel(v, grid)
        h = std(v) * n^(-1 / 6)
        return exp.(-0.5 .* ((grid' .- v) ./ h) .^ 2)
    end
    x, y = first.(pts), last.(pts)
    (std(x) > 0 && std(y) > 0) || return nothing
    d = kernel(x, xs)' * kernel(y, ys)
    d ./= sum(d)
    alpha < 1 && (d[d .< tpd_level(d, alpha)] .= 0)
    return d ./ sum(d)
end

# The overlap of two TPDs on the same grid, the summed minimum of the two in each cell
# (1 - Carmona's TPD dissimilarity); NaN if either is missing
tpd_overlap(d1, d2) = d1 === nothing || d2 === nothing ? NaN : sum(min.(d1, d2))

# Node => the TPD overlap of its two child clades in the trait space of `pts`
function trait_overlaps(tree, pts, nodes; alpha=0.95)
    grid = trait_grid(pts)
    return Dict(
        n => tpd_overlap((tpd(c, grid; alpha) for c in child_points(tree, n, pts))...) for
        n in nodes
    )
end

### ---- Clade densities in environmental space ---- ###

# The sites of `asm` as points in its space (for the environmental assemblage, the midpoints
# of its PC1-PC2 bins)
function site_points(asm)
    c = coordinates(asm)
    return Point2d.(c[:, 1], c[:, 2])
end

# Species => the indices of the sites of `asm` it occurs in
function species_sites(asm)
    occ = permutedims(occurrences(asm))
    return Dict(sp => findall(occ[:, j]) for (j, sp) in enumerate(speciesnames(asm)))
end

# For each of `node`'s two child clades, site index => the number of its species there (in
# `spsites`, species => site indices)
function child_site_counts(tree, node, spsites)
    return map(getchildren(tree, node)[1:2]) do c
        counts = Dict{Int,Int}()
        for sp in nodespecies(tree, getnodename(tree, c)), i in get(spsites, sp, Int[])
            counts[i] = get(counts, i, 0) + 1
        end
        return counts
    end
end

# The points of each of `node`'s two child clades: the union of the sites (in `sites`) of
# its species, each site counted once
function child_site_points(tree, node, sites, spsites)
    return [sites[sort!(collect(keys(c)))] for c in child_site_counts(tree, node, spsites)]
end

# Node => the overlap of the densities of its two child clades in the space of `asm`, each
# fit to the union of the sites of its species as the TPDs of `tpd`
function site_overlaps(tree, asm, nodes; alpha=0.95)
    sites = site_points(asm)
    spsites = species_sites(asm)
    grid = trait_grid(sites)
    return Dict(
        n => tpd_overlap(
            (tpd(c, grid; alpha) for c in child_site_points(tree, n, sites, spsites))...
        ) for n in nodes
    )
end

### ---- Divergence classes ---- ###

# The divergence classes, label => colour: which of geography, environment and traits the
# two child clades of a node diverge in. Indexed by 1 + geo + 2env + 4trait.
const DIVERGENCE_CLASSES = [
    "none" => RGBf(0.85, 0.85, 0.85),
    "geo" => RGBf(0 / 255, 114 / 255, 178 / 255),
    "env" => RGBf(0 / 255, 158 / 255, 115 / 255),
    "geo+env" => RGBf(86 / 255, 180 / 255, 233 / 255),
    "trait" => RGBf(240 / 255, 228 / 255, 66 / 255),
    "geo+trait" => RGBf(204 / 255, 121 / 255, 167 / 255),
    "env+trait" => RGBf(230 / 255, 159 / 255, 0 / 255),
    "all" => RGBf(0.1, 0.1, 0.1),
]

# Node => its divergence class (an index into DIVERGENCE_CLASSES): divergent in geography
# and in environment where its value in `values_g` and in `values_e` is above
# `threshold`, in traits where the trait overlap of its child clades in `overlap` is below
# `overlap_threshold`. Nodes missing any of the three are left out.
function divergence_classes(
    values_g, values_e, overlap; threshold=1.5, overlap_threshold=0.2
)
    nodes = [
        n for n in keys(overlap) if isfinite(overlap[n]) &&
        isfinite(get(values_g, n, NaN)) &&
        isfinite(get(values_e, n, NaN))
    ]
    return Dict(
        n =>
            1 +
            (values_g[n] > threshold) +
            2 * (values_e[n] > threshold) +
            4 * (overlap[n] < overlap_threshold) for n in nodes
    )
end

### ---- The trait explorer ---- ###

# `color` darkened by `f` (0 = black, 1 = unchanged)
function darken(color, f=0.6)
    c = to_color(color)
    return RGBAf(f * c.r, f * c.g, f * c.b, c.alpha)
end

# All points `allpts` in grey, with the points of the two child clades in `clades` (an
# Observable of the two point vectors) in the explorer's clade colours, the smaller clade on
# top, each with contours in a darker shade around the fractions `probs` of the probability
# of its density on `grid` (the outermost where `tpd` cuts it, by default). With
# `showpoints=false` the clades' points are left to the caller.
function clade_density_panel!(
    gp, allpts, grid, clades; probs=(0.95, 0.5, 0.25), showpoints=true, axis=(;)
)
    colors = clade_colors(:RdYlBu)
    ax = Axis(gp; xgridvisible=false, ygridvisible=false, axis...)
    scatter!(ax, allpts; color=:gray80, markersize=3, inspectable=false)
    for (k, color) in enumerate(colors)
        cladepts = lift(c -> c[k], clades)
        if showpoints
            sc = scatter!(ax, cladepts; color, markersize=5, inspectable=false)
            on(clades; update=true) do c
                translate!(sc, 0, 0, length(c[k]) <= length(c[3 - k]))
            end
        end
        # contours of the uncut TPD are smooth
        dens = lift(c -> something(tpd(c, grid; alpha=1), zeros(length.(grid))), cladepts)
        levels = lift(d -> iszero(d) ? [1.0] : sort([tpd_level(d, p) for p in probs]), dens)
        outline = contour!(
            ax, grid..., dens; levels, color=darken(color), linewidth=1, inspectable=false
        )
        translate!(outline, 0, 0, 2)
    end
    return ax
end

# All species in trait space, with the two child clades of `node` (an Observable) and the
# contours of their TPDs, as in `clade_density_panel!`
function trait_panel!(gp, asm, tree, node, x, y; probs=(0.95, 0.5, 0.25), axis=(;))
    pts = trait_points(asm, x, y)
    clades = lift(n -> child_points(tree, n, pts), node)
    return clade_density_panel!(
        gp, collect(values(pts)), trait_grid(pts), clades; probs, axis
    )
end

# The smallest spacing of the sites `pts` along x and along y (the size of the bins of the
# environmental space)
function bin_size(pts)
    spacing(v) = minimum(diff(sort!(unique(v))))
    return spacing(first.(pts)), spacing(last.(pts))
end

# All sites of `asm` in its space, with the sites of the two child clades of `node` (an
# Observable) as hollow circles shifted left and right within their bin, the area of each
# proportional to the number of the clade's species there (up to `maxsize` pixels across),
# and the contours of the clades' densities as in `clade_density_panel!`; the title gives
# the node's overlap in `overlap`
function site_panel!(
    gp, asm, tree, node, overlap; probs=(0.95, 0.5, 0.25), maxsize=12, axis=(;)
)
    sites = site_points(asm)
    spsites = species_sites(asm)
    dx, _ = bin_size(sites)
    counts = lift(n -> child_site_counts(tree, n, spsites), node)
    clades = lift(cs -> [sites[sort!(collect(keys(c)))] for c in cs], counts)
    title = lift(n -> "KDE overlap: $(round(get(overlap, n, NaN); digits=3))", node)
    ax = clade_density_panel!(
        gp,
        sites,
        trait_grid(sites),
        clades;
        probs,
        showpoints=false,
        axis=(; title, axis...),
    )
    marks = lift(counts) do cs
        cmax = maximum(c -> maximum(values(c); init=1), cs)
        return map(enumerate(cs)) do (k, c)
            s = sort!(collect(keys(c)))
            shift = Point2d((2k - 3) * dx / 5, 0)
            return (
                sites[s] .+ Ref(shift),
                [max(2, maxsize * sqrt(c[i] / cmax)) for i in s],
                [c[i] for i in s],
            )
        end
    end
    for (k, color) in enumerate(clade_colors(:RdYlBu))
        label(_, i, _) = "$(marks[][k][3][i]) species"
        scatter!(
            ax,
            lift(m -> m[k][1], marks);
            markersize=lift(m -> m[k][2], marks),
            color=:transparent,
            strokecolor=color,
            strokewidth=1,
            inspector_label=label,
        )
    end
    return ax
end

# Each node's log value in the environmental `values_e` against that in the geographic
# `values_g`, coloured by the trait overlap of its child clades in `overlap`, with the 1:1
# line; the node shown (an Observable) is ringed, and clicking a point shows that node if it
# is `selectable`, else says why not in `status`
function overlap_panel!(
    gp,
    node,
    overlap,
    values_g,
    values_e;
    selectable=n -> true,
    unselectable="",
    status=nothing,
    axis=(;),
)
    nodes = [
        n for n in keys(overlap) if
        isfinite(overlap[n]) && get(values_g, n, NaN) > 0 && get(values_e, n, NaN) > 0
    ]
    x = [log(values_g[n]) for n in nodes]
    y = [log(values_e[n]) for n in nodes]
    z = [overlap[n] for n in nodes]
    ax = Axis(gp[1, 1]; axis...)
    label(_, i, _) = "$(nodes[i])\ngeo: $(round(x[i]; digits=3))\nenv: $(round(y[i]; digits=3))\noverlap: $(round(z[i]; digits=3))"
    sc = scatter!(ax, x, y; color=z, markersize=5, inspector_label=label)
    Colorbar(gp[1, 2], sc; label="TPD overlap")
    ablines!(ax, 0, 1; color=:red, inspectable=false)
    index = Dict(zip(nodes, eachindex(nodes)))
    ring = lift(n -> haskey(index, n) ? [Point2d(x[index[n]], y[index[n]])] : Point2d[], node)
    scatter!(
        ax,
        ring;
        color=:transparent,
        strokecolor=:black,
        strokewidth=2,
        markersize=14,
        inspectable=false,
    )
    scene = ax.scene
    on(events(scene).mousebutton; priority=2) do ev
        (ev.button == Mouse.left && ev.action == Mouse.press) || return Consume(false)
        Makie.is_mouseinside(scene) || return Consume(false)
        p, i = pick(scene, events(scene).mouseposition[], 10)
        p === sc || return Consume(false)
        if selectable(nodes[i])
            node[] = nodes[i]
        elseif status !== nothing
            status[] = "$(nodes[i]): $unselectable"
        end
        return Consume(true)
    end
    return ax
end

# A slider under the tree of the explorer tree `tr` that marks the nodes `marked(t)` (node
# name => value) for the threshold `t` it is set to
function threshold_slider!(tr, marked; range=0:0.1:3, startvalue=2, label="threshold")
    sg = SliderGrid(tr.layout[4, 1], (; label, range, startvalue); tellwidth=false)
    on(sg.sliders[1].value) do t
        m = marked(t)
        Makie.update!(tr.treeplot; nodecolor=m, shownodes=collect(keys(m)))
    end
    return sg
end

# A node explorer of the two spaces and trait space: the tree, the SOS of the node shown in
# geographic and environmental space, its two child clades on the PCA axes `pcs`, and every
# node's log env `metric` against its log geo `metric`, coloured by the trait overlap of its
# child clades in `overlap` (computed on the same axes), where clicking a node shows it. If
# `env_overlap` (node => the overlap of the child clades in environmental space, from
# `site_overlaps`) is given, the last panel instead shows the two child clades in
# environmental space. Each space is passed as an (assemblage, NodeMetrics) pair. The
# tree marks `marked`, a Dict of node name => value, or, if `marked` is a function of a
# threshold returning one, the nodes for the threshold set with a slider of `thresholds`
# under the tree, starting at `threshold`. With `classes`, a vector of label => colour
# (e.g. DIVERGENCE_CLASSES), the values of `marked` are indices into it: the nodes are
# coloured by class, with the labels of the classes among them on the colour bar under the
# title `classlabel`.
function trait_explorer(
    tree,
    marked,
    nodevalues,
    (birds_g, res_g),
    (birds_e, res_e),
    explained,
    overlap;
    metric,
    pcs=(1, 2),
    threshold=2,
    thresholds=0:0.1:3,
    env_overlap=nothing,
    classes=nothing,
    classlabel="divergence class",
    images=nothing,
    imageoptions=(;),
)
    fig = Figure(; size=(1600, 850))
    marked0 = marked isa AbstractDict ? marked : marked(threshold)
    # the start node: the highest marked value, ties (classes) broken by `nodevalues`
    node = Observable(argmax(n -> (marked0[n], get(nodevalues, n, -Inf)), keys(marked0)))
    selectable(n) = has_sos(tree, res_g.sos, n) && has_sos(tree, res_e.sos, n)
    unselectable = "no SOS in both spaces"
    # with classes, the nodes are coloured by their position among the classes shown
    if classes === nothing
        shown, nodecolors, classcolors = nothing, marked0, (;)
    else
        shown = sort!(unique(values(marked0)))
        position = Dict(c => i for (i, c) in enumerate(shown))
        nodecolors = Dict(n => position[c] for (n, c) in marked0)
        classcolors = (;
            colormap=cgrad(last.(classes[shown]); categorical=true),
            colorrange=(0.5, length(shown) + 0.5),
        )
    end
    tr = explorer_tree!(
        fig[1, 1],
        tree,
        node,
        nodecolors;
        values=nodevalues,
        label="geo $metric",
        selectable,
        unselectable,
        images,
        imageoptions,
        rangesize=birds_g,
        classcolors...,
    )
    if classes !== nothing
        cb = only(contents(tr.layout[3, 1]))
        cb.ticks = (eachindex(shown), first.(classes[shown]))
        cb.label = classlabel
        cb.width = Relative(0.9)
        hover = tr.treeplot.hoverlabel[]
        tr.treeplot.hoverlabel = function (n)
            haskey(marked0, n) || return hover(n)
            return hover(n) * "\n" * first(classes[marked0[n]])
        end
    end
    if !(marked isa AbstractDict)
        threshold_slider!(
            tr, marked; range=thresholds, startvalue=threshold, label="$metric threshold"
        )
    end
    panels = fig[1, 2] = GridLayout()
    sos_map!(panels[1, 1], birds_g, node, res_g; title="Geographic SOS")
    sos_map!(panels[1, 2], birds_e, node, res_e; title="Environmental SOS")
    pc_label(i) = "pca$i ($(round(100explained[i]; digits = 1))%)"
    trait_panel!(
        panels[2, 1],
        birds_g,
        tree,
        node,
        (Symbol("pca$i") for i in pcs)...;
        axis=(; xlabel=pc_label(pcs[1]), ylabel=pc_label(pcs[2])),
    )
    if env_overlap !== nothing
        site_panel!(
            panels[2, 2],
            birds_e,
            tree,
            node,
            env_overlap;
            axis=(; xlabel="env PC1", ylabel="env PC2"),
        )
    else
        overlap_panel!(
            panels[2, 2],
            node,
            overlap,
            getfield(res_g, metric),
            getfield(res_e, metric);
            selectable,
            unselectable,
            status=tr.status,
            axis=(;
                xlabel="log geo $metric",
                ylabel="log env $metric",
                title="Coloured by TPD overlap (pca$(pcs[1])-$(pcs[2]))",
            ),
        )
    end
    colsize!(fig.layout, 1, Relative(0.45))
    DataInspector(fig)
    return fig, tr
end
