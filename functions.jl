"""
    string_sites!(df)

Convert the `site` column of `df` to strings, so site ids read by CSV stay strings.
Returns `df`.
"""
function string_sites!(df)
    df.site = string.(df.site)
    return df
end

# The trait rows (from `addtraits!`) of the species in `asm` that descend from `node`
function _clade_traits(asm, tree, node)
    species = Set(nodespecies(tree, node))
    return filter(r -> r.name in species, traits(asm))
end

"""
    genera(asm, tree, node)

The sorted genera of the species in `asm` that descend from `node` in `tree`. The genus is
the first part of the species name, e.g. "Accipiter" of "Accipiter_badius".
"""
function genera(asm, tree, node)
    species = _clade_traits(asm, tree, node).name
    return sort!(unique(String.(first.(split.(species, "_")))))
end

"""
    families(asm, tree, node)

The sorted families (the AVONET `Family1` trait) of the species in `asm` that descend from
`node` in `tree`.
"""
function families(asm, tree, node)
    return sort!(unique(String.(_clade_traits(asm, tree, node).Family1)))
end

# `taxa` as lines of `perline`, cut to the first `n`
function _taxa_list(taxa, n; perline=4)
    shown = first(taxa, n)
    lines = [
        join(shown[i:min(i + perline - 1, end)], ", ") for i in 1:perline:length(shown)
    ]
    length(taxa) > n && push!(lines, "… and $(length(taxa) - n) more")
    return join(lines, "\n")
end

"""
    taxa_text(asm, tree, node; nfamilies=8, ngenera=6)

The taxa of the species in `asm` below `node`, for a hover label: the first `nfamilies`
families if there are several, else the family and its first `ngenera` genera.
"""
function taxa_text(asm, tree, node; nfamilies=8, ngenera=6)
    fams = families(asm, tree, node)
    isempty(fams) && return ""
    length(fams) > 1 && return "$(length(fams)) families:\n" * _taxa_list(fams, nfamilies)
    return only(fams) * "\n" * _taxa_list(genera(asm, tree, node), ngenera)
end

"""
    taxa_hover!(explorer, asm, tree)

Add [`taxa_text`](@ref) to the hover labels of the tree of `explorer` (from
`node_explorer` or `explorer_tree!`), and of its ordination if it has one. Returns
`explorer`.
"""
function taxa_hover!(explorer, asm, tree)
    base = explorer.treeplot.hoverlabel[]
    cache = Dict{String,String}()
    label(n) = get!(() -> string(base(n), "\n", taxa_text(asm, tree, n)), cache, n)
    explorer.treeplot.hoverlabel = label
    if hasproperty(explorer, :ordination) && explorer.ordination !== nothing
        explorer.ordination.hoverlabel = label
    end
    return explorer
end

"""
    taxa_image_hover!(explorer, asm)

Add the family and order under the species name in the hover labels of the species images
of `explorer` (from `node_explorer` or `explorer_tree!`): those around the tree, and those
in the corner of the child-clade maps. Returns `explorer`.

NodivMakie has no option for the image labels, so this re-attaches its internal
`_image_hover!` with the new text (a replacement label must be of the same type as the old
one).
"""
function taxa_image_hover!(explorer, asm)
    t = traits(asm)
    taxon = Dict(zip(t.name, zip(t.Family1, t.Order1)))
    function relabel(text)
        lines = split(text, "\n")
        sp = replace(first(lines), " " => "_")
        haskey(taxon, sp) || return text
        fam, ord = taxon[sp]
        return join([first(lines); "$fam, $ord"; lines[2:end]], "\n")
    end
    if explorer.images !== nothing
        scene = explorer.images.axis.scene
        for p in explorer.images.plots
            p isa Image || continue
            NodivMakie._image_hover!(p, scene, relabel(p.inspector_label[](p, 1, nothing)))
        end
    end
    hasproperty(explorer, :panel) || return explorer
    # the child-clade images change with the node shown
    for ax in explorer.panel.axes[3:4]
        ax === nothing && continue
        for p in ax.scene.plots
            p isa Image || continue
            old = p.inspector_label[]
            text = lift(_ -> relabel(old(p, 1, nothing)), explorer.panel.node)
            NodivMakie._image_hover!(p, ax.scene, text)
        end
    end
    return explorer
end
