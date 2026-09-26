# Design spec: grouping divergent nodes by SOS-pattern similarity

**Context.** Node-based analysis (Borregaard et al. 2014) applied to New World birds in two spaces — geographic (`birds_g`, `res_g`) and environmental (`birds_e`, `res_e`). After selecting strongly divergent nodes (`divergent`, the `divergent_e ∩ divergent_g` intersection, n ≈ 100), the goal is to **identify groups of nodes that exhibit very similar SOS patterns**. The current approach computes a correlation-based distance matrix (`sos_distances` in `Nodiv`) and feeds it to classical MDS (`fit(MDS, …; distances = true, maxoutdim = 2)`), plotted by `sos_mds_plot` in `script.jl`. This document specifies the corrected approach and the reasoning behind it. It is a design spec, not an implementation; the implementer should read the current body of `sos_distances` first, since the parameter names and the `res.sos[node]` access pattern below are inferred, not confirmed.

## The goal, stated precisely

Find clusters of nodes whose per-cell SOS maps are mutually similar, against a background of nodes whose patterns are mutually distinct. "Similar pattern" must treat the two daughter labels symmetrically (see sign, below), and must not reward pairs that merely co-occupy space without sharing the same over/under-representation structure. A valid outcome includes finding *few or no* tight groups — that the divergent nodes are mostly idiosyncratic is a legitimate scientific result, not a failure to be engineered away.

## What the current MDS picture is actually telling us

The circular/ring arrangement with little clustering is not a weak or failed embedding. It is an accurate report. Classical MDS embeds a matrix of near-uniform, near-maximal distances by spreading the points evenly on a circle, because that is the 2-D configuration in which all pairwise distances are roughly equal. The ring therefore means: *almost all of these nodes are mutually near-orthogonal in pattern, and there is little low-dimensional structure to spatialize.* This follows directly from the empirical distance distribution (next section). A different embedder will not "recover" groups that the distances do not support; it will only redistribute the same near-equidistance differently.

If a demonstration is wanted, fit MDS once at `maxoutdim ≈ min(10, n−1)` and inspect the eigenvalues. If axes 3+ carry weight comparable to axes 1–2, the 2-D scatter is a projection artefact, which confirms the ring reading. This is a one-time diagnostic, not the analysis.

## Agreed design principles

**1. Sign of SOS is arbitrary per node, so use 1 − |r|, not 1 − r.** Which daughter is "clade 1" is incidental to the node, so a mirror-image SOS map represents the *same* divergence geography with the labels swapped. Folding by absolute value makes positive and negative correlations equally indicative of pattern similarity, which is what the question requires. (In this dataset there are essentially no strong negative correlations, so in practice |r| and r nearly coincide; |r| remains the principled choice regardless.)

**2. Non-overlap is real difference, not missing information.** Two clades occupying disjoint regions are maximally divergent in *where* their over/under-representation falls. Encoding a disjoint pair as distance 1.0 is the honest encoding, and the resulting pile-up of distances near 1.0 is true structure, not an artefact. Methods that hide this structure (see "what to avoid") would misrepresent the data. The implication is the opposite of smoothing it away: the analysis must preserve the distinction between "genuinely similar" and "everything else," rather than compressing the background into a manufactured gradient.

**3. The empirical distance range is compressed toward the maximum.** Pattern correlations run from roughly random up to strongly similar, with essentially nothing more anticorrelated than chance. On the 1 − |r| distance this puts the few similar pairs near 0, unrelated-but-overlapping pairs near 1 − |r_random|, and disjoint pairs at exactly 1 — i.e. a thin tail of close pairs against a dense band near the maximum. There is no genuine "far pole" to pull a low-dimensional structure out, which is exactly why the embedding rings. The right question is therefore not "which 2-D method recovers the groups" but "are there any tight groups at all, or only a handful of similar pairs against a sea of orthogonality" — a question best answered on the similarity matrix directly, not on any spatialization of it.

**4. No ancestor–descendant collapse.** SOS at a node compares *that node's two daughters'* relative occupancy under a null; a node and its parent are built from different partitions of different species sets, so there is no design-level reason for them to share a pattern. If nested nodes do cluster, that is a contingent finding (a divergence concentrated in one descendant lineage and persisting down the tree), to be inspected on the phylogeny — not a pseudoreplication artefact to remove a priori. Do not exclude or down-weight ancestor–descendant pairs.

## Distance definition

Distance between nodes *k* and *l*:

```
D(k,l) = 1 − |r|, computed over cells where both nodes are occupied,
         provided the shared occupied support is at least `minoverlap` cells;
         otherwise D(k,l) = 1.
```

Specifics:

- **Support = occupied cells.** Define occupancy per node explicitly (clade present), not as `SOS != 0`. Under the paper's definition SOS ≈ 0 at occupied cells where both daughters are equally represented, so an `SOS != 0` mask would silently drop equal-representation cells in addition to absent ones — a third conflation on top of the two named above. Confirm which mask the current `sos_distances` uses; this is the most likely latent bug.
- **Minimum-overlap guard.** Below `minoverlap` shared occupied cells, set distance to 1 rather than trusting a correlation estimated on a handful of cells. This handles the spurious-high-correlation-on-tiny-overlap case by rule. A floor on the order of 5–10 cells is reasonable for the geographic scan; choose by the precision you are willing to defend, and note that the environmental scan has far fewer bins so the floor must be set separately per space.
- **|Spearman| over |Pearson| if SOS is heavy-tailed,** which it tends to be at strongly divergent nodes where a few cells carry extreme values. Worth checking the marginal SOS distributions and choosing accordingly.

**One genuine judgment call — whether to weight by overlap extent.** Correlation conditions *on* the shared support and is then blind to how large it is: two nodes sharing three cells with perfect correlation get the same distance as two sharing their whole range. Given principle 2 (partial overlap = partial comparability), there is a case for multiplying by an overlap term, `D = 1 − O·|r|` with `O` a Jaccard/Sørensen/Simpson index of occupied-cell overlap, so that partial-overlap-but-high-correlation does not masquerade as full similarity. The cost is that this re-couples overlap and pattern into one number. Recommendation: compute both the overlap term and |r| and keep them inspectable; default to the overlap-weighted distance for the clustering, but retain the unweighted version so the two can be compared. The downstream method (below) is robust to this choice, which is part of why it is preferred.

## Recommended method

**Primary: hierarchical clustering on the distance matrix, read as a clustered heatmap with dendrogram.** At n ≈ 100 this is fully legible and it sidesteps spatialization entirely: a heatmap never has to reconcile the disjoint 1.0s against the rest, it simply shows them as the uniform background they are. Reorder rows/columns by the dendrogram and the dense blocks (real groups of co-patterned nodes), if any, appear on the diagonal while the orthogonal background stays uniform. This directly answers the question in principle 3.

- **Linkage: average or complete; not Ward.** Ward assumes roughly spherical Euclidean clusters and will impose block structure that is not there. Complete linkage is the conservative default — it groups only all-pairs-similar nodes and will not chain marginal pairs into spurious groups, which matters when the similar tail is thin.
- **Read groups by cutting the dendrogram at a justifiable height,** e.g. the height corresponding to |r| ≈ 0.7, rather than by eyeballing a scatter. State the cut threshold and its rationale in the methods.
- Map the resulting groups back onto the phylogeny and onto the SOS maps for interpretation, including any nested nodes that co-cluster.

**Secondary / confirmatory: thresholded similarity graph + community detection (Leiden or Louvain).** Build edges only between pairs with |r| above a high threshold (and overlap above the floor); disjoint and orthogonal pairs never become edges and so drop out of the structure rather than repelling anything, leaving genuinely co-patterned nodes as connected components/communities. Run this *after* the heatmap: the heatmap tells you whether there is enough block structure to make community detection worthwhile. If there is, the graph gives a cleaner group assignment and scales if the node set later grows.

**Alternative single-distance route: HDBSCAN on the (overlap-weighted) distance matrix.** Density-based, finds dense groups and labels the remainder as noise — which is the honest treatment of lonely nodes and matches principle 2/3. Use if a single-distance pipeline is preferred over the graph.

## What to avoid, and why

- **UMAP / t-SNE at this n.** Two distinct problems. First, at n ≈ 100 they are unreliable and manufacture clusters from noise. Second, they build *relative* neighborhoods — every node gets its k nearest regardless of absolute distance — so a node whose nearest neighbors sit at 0.9 is handed a fake neighborhood, exactly the misrepresentation principle 2 warns against. The goal needs an *absolute* notion of "similar enough," which thresholded graphs and HDBSCAN provide and neighbor-embeddings do not. Reserve UMAP for a future regime with hundreds of nodes, and even then only as a viewer of communities defined on the true graph, never as the definition of the groups.
- **Letting the 2-D MDS scatter define groups.** Keep MDS only as the eigenvalue diagnostic above, if at all.
- **Ward linkage** (imposes spherical structure) and **`SOS != 0` as the overlap mask** (drops equal-representation cells), per above.

## Two-space interpretation

The environmental and geographic distances are not comparable in magnitude — environmental "occupancy" is over tens of climate bins, geographic over ~18k cells — so do not read the two ordinations/heatmaps as if their axes or distances share units. Treat them as two different questions: the geographic analysis asks whether two nodes mark *the same specific break*; the environmental analysis asks whether they mark *the same kind of transition, abstracted from location*. That asymmetry is information. A pair similar in environmental but not geographic space marks the same climatic transition in different places — itself a substantive pattern.

The environmental MDS also rings, which most likely means the same thing there: even though clades share climate bins, their SOS *patterns over* those bins are mostly mutually orthogonal, because co-occupancy of a bin does not imply the same daughter dominates it. Confirm on the environmental heatmap; expect "mostly uniform background + a few hot pairs" there too.

## Reporting

If the heatmaps show mostly uniform background with a few off-diagonal hot pairs and no substantial blocks, that null result — *the divergent nodes are largely idiosyncratic in SOS pattern, with a small number of co-patterned exceptions (name them)* — is the scientific finding and should be reported as such, with the heatmap as the evidence. Do not escalate to more aggressive embedding or clustering in order to produce groups the distances do not support.

## Implementation touch-points

- `Nodiv`: `sos_distances(res, nodes; …)` — the distance definition lives here. Read its current body first; verify the overlap mask and the sign handling, add the minimum-overlap guard and |r| folding, and expose the overlap-weighting option. Keep the existing call signature so `script.jl` is unaffected.
- `script.jl`: replace `sos_mds_plot` with a clustered-heatmap-plus-dendrogram routine as the primary view; keep an MDS eigenvalue diagnostic available; add the graph/community-detection path as a secondary. Run separately for `res_e` and `res_g` on the `divergent` node set, and do not cross-compare magnitudes between the two.
