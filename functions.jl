"""
    strsite!(df)

Convert the `site` column of `df` to strings, so site ids read by CSV stay strings.
Returns `df`.
"""
function strsite!(df)
    df.site = string.(df.site)
    return df
end

"""
    mapfigure(args...; title="", label="", kw...)

A map of one value per site (`sitemap(args...)`), with a colour bar labelled `label` beside
it. Other keywords go to `sitemap`. Returns the figure.
"""
function mapfigure(args...; title="", label="", kw...)
    fig, ax, p = sitemap(args...; axis=(; title), kw...)
    Colorbar(fig[1, 2], p; label)
    return fig
end

"""
    cluster_idmap(groups)

The clusters worth showing: the cut clusters of `groups` with more than one node,
relabelled 1..m in the order of their `cutree` ids, as a `Dict` from cluster id to label.
The heatmap and the tree both number and colour the clusters by this, so the two views
cross-reference directly.
"""
function cluster_idmap(groups)
    ids = sort(collect(keys(sos_cluster_sizes(groups).members)))
    return Dict(c => i for (i, c) in enumerate(ids))
end

"""
    clustercolors(m)

`m` colours for the labelled clusters, cycling through the `:tab20` colormap.
"""
function clustercolors(m)
    c = Makie.to_colormap(:tab20)
    return [c[mod1(i, length(c))] for i in 1:m]
end

"""
    sos_cluster_sizes(groups)

Cluster-size table for a `groups` Dict (node => cluster): how many nodes fall in each cut
cluster, and which clusters are non-trivial (size > 1). A near-flat table of singletons is
the "largely idiosyncratic" null result; a few multi-node clusters are the co-patterned
exceptions.

Returns a NamedTuple with `nclusters`, `nsingletons` and `members`, a `Dict` from each
non-trivial cluster to its sorted node names.
"""
function sos_cluster_sizes(groups)
    counts = Dict{Int,Int}()
    for c in values(groups)
        counts[c] = get(counts, c, 0) + 1
    end
    nontrivial = [c for (c, k) in counts if k > 1]
    members = Dict(c => sort([n for (n, g) in groups if g == c]) for c in nontrivial)
    return (; nclusters=length(counts), nsingletons=count(==(1), values(counts)), members)
end

"""
    plot_cluster_tree(tree, groups, title; markersize=12, kw...)

Map the heatmap clusters onto the phylogeny: a marker at each node of a non-trivial cluster
(size > 1), one colour per cluster, and nothing elsewhere, so the co-patterned groups stand
out against the tree; idiosyncratic singletons are left unmarked. `groups` is the Dict from
`sos_clusters`. The cluster numbers and colours are those outlined on the heatmap
(`cluster_idmap`). Other keywords go to `treeplot`. Returns the figure.
"""
function plot_cluster_tree(tree, groups, title; markersize=12, kw...)
    idmap = cluster_idmap(groups)
    shown = Dict(n => idmap[c] for (n, c) in groups if haskey(idmap, c))
    fig, ax, p = treeplot(
        tree;
        treetype=:fan,
        showtips=false,
        nodegroup=shown,
        groupcolors=clustercolors(length(idmap)),
        markersize,
        strokewidth=0.5,
        strokecolor=:gray30,
        axis=(; title),
        figure=(; size=(900, 800)),
        kw...,
    )
    Legend(fig[1, 2], ax, "cluster")
    return fig
end

"""
    plot_node_pdf(assemblage, tree, nodes, res, outfile)

Write one 4-panel node panel (parent / SOS / child 1 / child 2) per node of `nodes` to
`outfile`, a single multi-page PDF with one node per page. The SOS panel comes from the
cached `res`, with no recompute. Returns `outfile`.

The panel is built once and switched from node to node; each page is saved with CairoMakie
(GLMakie cannot write PDF) and the pages merged with `pdfunite`, from poppler, which must
be on the PATH.
"""
function plot_node_pdf(assemblage, tree, nodes, res, outfile)
    if Sys.which("pdfunite") === nothing
        error(
            "plot_node_pdf needs `pdfunite` (poppler) on PATH, e.g. `brew install poppler`"
        )
    end
    tmp = mktempdir()
    pages = String[]
    fig, panel = nodepanel(assemblage, tree, first(nodes), res)
    for (i, node) in enumerate(nodes)
        panel.node[] = node
        page = joinpath(tmp, string(lpad(i, 3, '0'), ".pdf"))
        save(page, fig; backend=CairoMakie)
        push!(pages, page)
    end
    run(`pdfunite $pages $outfile`)
    rm(tmp; recursive=true)
    @info "wrote node-panel PDF" outfile npages = length(pages)
    return outfile
end
