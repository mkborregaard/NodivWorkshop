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

# A grid over the trait space of `pts` (species => point): `n` points along each axis,
# spanning the range of the points extended by `extend` of it on each side
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

### ---- The trait explorer ---- ###

# `color` darkened by `f` (0 = black, 1 = unchanged)
function darken(color, f=0.6)
    c = to_color(color)
    return RGBAf(f * c.r, f * c.g, f * c.b, c.alpha)
end

# All species in trait space in grey, with the two child clades of `node` (an Observable) in
# the explorer's clade colours, the smaller clade on top, each with contours in a darker
# shade around the fractions `probs` of the probability of its TPD (the outermost where
# `tpd` cuts it, by default)
function trait_panel!(gp, asm, tree, node, x, y; probs=(0.95, 0.5, 0.25), axis=(;))
    pts = trait_points(asm, x, y)
    grid = trait_grid(pts)
    colors = clade_colors(:RdYlBu)
    ax = Axis(gp; xgridvisible=false, ygridvisible=false, axis...)
    scatter!(ax, collect(values(pts)); color=:gray80, markersize=3, inspectable=false)
    clades = lift(n -> child_points(tree, n, pts), node)
    for (k, color) in enumerate(colors)
        cladepts = lift(c -> c[k], clades)
        sc = scatter!(ax, cladepts; color, markersize=5, inspectable=false)
        on(c -> translate!(sc, 0, 0, length(c[k]) <= length(c[3 - k])), clades; update=true)
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

# Each node's log value in `values` against its trait overlap in `overlap`, with the
# least-squares line; the node shown (an Observable) is ringed, and clicking a point shows
# that node if it is `selectable`, else says why not in `status`
function overlap_panel!(
    gp, node, overlap, values; selectable=n -> true, unselectable="", status=nothing, axis=(;)
)
    nodes = [
        n for n in keys(overlap) if isfinite(overlap[n]) && get(values, n, NaN) > 0
    ]
    x = [overlap[n] for n in nodes]
    y = [log(values[n]) for n in nodes]
    ax = Axis(gp; axis...)
    label(_, i, _) = "$(nodes[i])\noverlap = $(round(x[i]; digits=3))"
    sc = scatter!(ax, x, y; color=:gray50, markersize=5, inspector_label=label)
    ablines!(ax, ([ones(length(x)) x] \ y)...; color=:red, inspectable=false)
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

# A node explorer of the two spaces and trait space: the tree, the SOS of the node shown in
# geographic and environmental space, its two child clades on PCA axes 1-2, and every
# node's log `nodevalues` against the trait overlap of its child clades in `overlap`, where
# clicking a node shows it. Each space is passed as an (assemblage, NodeMetrics) pair.
function trait_explorer(
    tree,
    marked,
    nodevalues,
    (birds_g, res_g),
    (birds_e, res_e),
    explained,
    overlap;
    metric,
    images=nothing,
    imageoptions=(;),
)
    fig = Figure(; size=(1600, 850))
    node = Observable(argmax(n -> marked[n], keys(marked)))
    selectable(n) = has_sos(tree, res_g.sos, n) && has_sos(tree, res_e.sos, n)
    unselectable = "no SOS in both spaces"
    tr = explorer_tree!(
        fig[1, 1],
        tree,
        node,
        marked;
        values=nodevalues,
        label="geo $metric",
        selectable,
        unselectable,
        images,
        imageoptions,
        rangesize=birds_g,
    )
    panels = fig[1, 2] = GridLayout()
    sos_map!(panels[1, 1], birds_g, node, res_g; title="Geographic SOS")
    sos_map!(panels[1, 2], birds_e, node, res_e; title="Environmental SOS")
    pc_label(i) = "pca$i ($(round(100explained[i]; digits = 1))%)"
    trait_panel!(
        panels[2, 1],
        birds_g,
        tree,
        node,
        :pca1,
        :pca2;
        axis=(; xlabel=pc_label(1), ylabel=pc_label(2)),
    )
    overlap_panel!(
        panels[2, 2],
        node,
        overlap,
        nodevalues;
        selectable,
        unselectable,
        status=tr.status,
        axis=(; xlabel="TPD overlap (pca1-2)", ylabel="log geo $metric"),
    )
    colsize!(fig.layout, 1, Relative(0.45))
    DataInspector(fig)
    return fig, tr
end
