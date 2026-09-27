# Tests of the functions in functions.jl, on a small synthetic tree and assemblage. The
# hover labels of the species images rely on NodivMakie internals, so these tests catch a
# NodivMakie update that breaks them.
#
#     julia --project=. test/runtests.jl

using CairoMakie
using DataFrames
using Nodiv
using NodivMakie
using Phylo
using SpatialEcology
using Test

include(joinpath(@__DIR__, "..", "functions.jl"))

const NEWICK = "((Aa_x:1,Aa_y:1)n1:1,(Bb_x:1.5,(Cc_x:1,Cc_y:1)n3:0.5)n2:1)root;"
const SPECIES = ["Aa_x", "Aa_y", "Bb_x", "Cc_x", "Cc_y"]

tree = parsenewick(NEWICK)

# Five species on a 4 x 3 grid; Aa in Fam1 (Ord1), Bb and Cc in Fam2 (Ord2)
grid = [(x, y) for y in 1:3 for x in 1:4]
occ = zeros(Int, 5, 12)
occ[1, 1:6] .= 1
occ[2, 4:9] .= 1
occ[3, 7:12] .= 1
occ[4, [1, 5, 9]] .= 1
occ[5, 10:12] .= 1
sites = ["s$i" for i in 1:12]
asm = Assemblage(occ, Float64[first.(grid) last.(grid)], sites, SPECIES)
avonet = DataFrame(;
    species=SPECIES,
    Family1=["Fam1", "Fam1", "Fam2", "Fam2", "Fam2"],
    Order1=["Ord1", "Ord1", "Ord2", "Ord2", "Ord2"],
)
addtraits!(asm, avonet, :species)

internal = ["root", "n1", "n2", "n3"]
sos = Dict(n => collect(range(-8, 8; length=12)) .* k for (k, n) in enumerate(internal))
res = NodeAnalysis(
    internal, Dict("root" => 0.2, "n1" => 0.4, "n2" => 0.9, "n3" => 0.5), sos
)

# A synthetic image of every species
imagedir = mktempdir()
for sp in SPECIES
    img = fill(RGBAf(0.2, 0.4, 0.6, 1), 40, 40)
    NodivMakie.FileIO.save(joinpath(imagedir, "$sp.png"), img)
end

# The labels of the species images of `plots`, as the DataInspector shows them
image_labels(plots) = [p.inspector_label[](p, 1, nothing) for p in plots if p isa Image]
function clade_image_labels(ex)
    return reduce(vcat, [image_labels(ax.scene.plots) for ax in ex.panel.axes[3:4]])
end

@testset "NodivWorkshop functions" begin
    @testset "genera and families" begin
        @test genera(asm, tree, "root") == ["Aa", "Bb", "Cc"]
        @test genera(asm, tree, "n2") == ["Bb", "Cc"]
        @test genera(asm, tree, "Cc_y") == ["Cc"]
        @test families(asm, tree, "root") == ["Fam1", "Fam2"]
        @test families(asm, tree, "n2") == ["Fam2"]
    end

    @testset "taxa_text" begin
        @test _taxa_list(string.('a':'c'), 8) == "a, b, c"
        @test _taxa_list(string.('a':'j'), 8) == "a, b, c, d\ne, f, g, h\n… and 2 more"
        @test taxa_text(asm, tree, "root") == "2 families:\nFam1, Fam2"
        @test taxa_text(asm, tree, "n2") == "Fam2\nBb, Cc"
        @test taxa_text(asm, tree, "n2"; ngenera=1) == "Fam2\nBb\n… and 1 more"
        @test taxa_text(asm, tree, "n3") == "Fam2\nCc"
    end

    @testset "node_explorer hover labels" begin
        fig, ex = node_explorer(asm, tree, res; nodes=:all, images=imagedir)
        @test ex.ordination !== nothing
        before = ex.treeplot.hoverlabel[]("n2")
        treelabels = image_labels(ex.images.plots)
        @test !isempty(treelabels)

        taxa_hover!(ex, asm, tree)
        taxa_image_hover!(ex, asm)
        @test ex.treeplot.hoverlabel[]("n2") == before * "\nFam2\nBb, Cc"
        @test endswith(ex.treeplot.hoverlabel[]("root"), "\n2 families:\nFam1, Fam2")
        @test ex.ordination.hoverlabel[]("root") == ex.treeplot.hoverlabel[]("root")

        # tree images: the family and order under the species name, the rest kept
        for (old, new) in zip(treelabels, image_labels(ex.images.plots))
            lines, newlines = split(old, "\n"), split(new, "\n")
            sp = replace(first(lines), " " => "_")
            fam, ord = only(eachrow(avonet[avonet.species .== sp, :]))[[:Family1, :Order1]]
            @test newlines == [first(lines); "$fam, $ord"; lines[2:end]]
        end
        # the hover callbacks resolve (a replacement of another type would not)
        for p in ex.images.plots
            p isa Image && @test p.inspector_hover[] isa Function
        end

        # child-clade images follow the node shown
        ex.panel.node[] = "n1"
        @test sort(clade_image_labels(ex)) == ["Aa x\nFam1, Ord1", "Aa y\nFam1, Ord1"]
        ex.panel.node[] = "n3"
        @test sort(clade_image_labels(ex)) == ["Cc x\nFam2, Ord2", "Cc y\nFam2, Ord2"]
        @test size(Makie.colorbuffer(fig)) != (0, 0)
    end

    @testset "color_by_clusters!" begin
        fig, ex = node_explorer(asm, tree, res; nodes=:all)
        nodes = ex.ordination.ordination[].nodes
        clusters = sos_clusters(sos_distances(res, nodes), nodes; simcut=0.5)
        color_by_clusters!(ex, clusters)
        colors = cluster_colors(length(clusters.labels))
        for (n, c) in clusters.groups
            expected = haskey(clusters.labels, c) ? colors[clusters.labels[c]] : :gray70
            @test ex.ordination.nodecolor[][n] == to_color(expected)
        end
        i = findfirst(==("n2"), nodes)
        @test ex.ordination.point_colors[][i] == ex.ordination.nodecolor[]["n2"]
        @test size(Makie.colorbuffer(fig)) != (0, 0)
    end

    @testset "explorer_tree! hover labels" begin
        # the tree of an explorer with other panels, as in trait_explorer
        fig = Figure()
        et = explorer_tree!(
            fig[1, 1],
            tree,
            Observable("root"),
            Dict("n1" => 0.4, "n2" => 0.9);
            images=imagedir,
            rangesize=asm,
        )
        taxa_hover!(et, asm, tree)
        taxa_image_hover!(et, asm)
        @test endswith(et.treeplot.hoverlabel[]("n3"), "\nFam2\nCc")
        labels = image_labels(et.images.plots)
        @test !isempty(labels) && all(l -> occursin(r"\nFam\d, Ord\d", l), labels)
    end
end
