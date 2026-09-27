# Design: grouping divergent nodes by SOS-pattern similarity

**Context.** Node-based analysis (Borregaard et al. 2014) applied to New World birds in two spaces — geographic (`birds_g`, `res_g`) and environmental (`birds_e`, `res_e`). Each space has its own set of strongly divergent nodes (RMS-SOS > 2): `divergent_e` (92 nodes) and `divergent_g` (38 nodes). The earlier intersection `divergent = divergent_e ∩ divergent_g` is no longer used. The goal is to **identify groups of nodes that exhibit very similar SOS patterns**, separately in each space. This document records the method decided on 2026-09-27 and the reasoning behind it, including the alternatives that were tried and rejected.

## The decision

- **Distance:** `D = 1 − |r|` from `sos_distances` in Nodiv (unweighted, Pearson, `minoverlap = 3`).
- **Grouping:** complete-linkage hierarchical clustering (`sos_clusters`), cut at an absolute similarity **|r| ≥ 0.6** (`SIMCUT` in `script.jl`). Every divergent node is in exactly one cluster; a cluster of one node is allowed. The cut is a choice the user owns and states in the methods. The data do not choose it.
- **Output:** the printed cluster membership, and the clusters on the phylogeny (`cluster_tree`). The dendrogram-ordered |r| heatmap (`sos_cluster_heatmap`) is optional.
- **Explorer:** the ordination panel of the node explorers (2-D classical MDS of the same `D`) stays, coloured by the clusters (`color_by_clusters!`). It is a projection for browsing, so distances in it are approximate. The clusters come from the full `D`, not from the 2-D picture.
- **Dropped:** the MDS eigenvalue plot, the thresholded similarity-graph communities (`sos_similarity_communities`), and the bootstrap-split clusters (`supported_clusters`).

What complete linkage at a cut means is easy to state: *every pair of nodes in a group has |r| ≥ 0.6.* The number of clusters then answers "how many groups of nodes share an SOS pattern at similarity ≥ 0.6". It does not answer "how many distinct patterns there are". The data do not support a method-independent answer to the latter (see below).

## Principles (unchanged)

**1. Sign of SOS is arbitrary per node, so use 1 − |r|, not 1 − r.** Which daughter is "clade 1" is incidental to the node. A mirror-image SOS map is the *same* divergence geography with the labels swapped.

**2. Non-overlap is real difference, not missing information.** Two clades in disjoint regions are maximally different in where their over- and under-representation falls, so a disjoint pair gets distance 1. The same holds for pairs sharing fewer than `minoverlap` cells with both SOS defined.

**3. No ancestor–descendant exclusion.** SOS at a node compares that node's two daughters. A node and its parent are built from different partitions of different species sets, so there is no design-level reason for them to share a pattern. The data bear this out: nested (ancestor–descendant) pairs are 22% of the similar pairs and 22% of all pairs in the environmental space. So co-clustered nodes mostly mark shared geography, not shared ancestry. Do not exclude or down-weight nested pairs.

**4. A result of few or no groups is legitimate.** If the divergent nodes are largely idiosyncratic in SOS pattern, that is the finding, not a failure to be engineered away.

## Distance, as implemented

`sos_distances(res, nodes; minoverlap=3)` correlates two SOS maps over the cells where both are finite: where both nodes' daughter clades are present and the null model varies. Occupied cells where the null model cannot vary have no SOS and do not count. A pair sharing fewer than `minoverlap` such cells gets distance 1, and a constant map gets distance 1. The options exist but are not used: `method = :spearman`, and `overlapweight = true` (`D = 1 − O·|r|`, with `O` the Sørensen index of the occupied cells).

**The spaces differ in size, and |r| is not comparable between them.** The environmental space has 489 climate bins, and the geographic space 17,542 cells. (Earlier versions of this document and of `script.jl` said "tens of PC bins"; that was wrong.)
- **Environmental:** pairs of divergent nodes share a median of 310 bins with SOS defined (5th percentile 71). Co-clustered pairs share a median of about 350. At 300 bins, r = 0.7 has a 95% CI of about 0.64–0.75, so the |r| values are not sampling noise.
- **Distributions:** environmental median |r| is 0.39 (90th percentile 0.74), with no pairs at distance 1. Geographic median |r| is 0.31 (90th percentile 0.81), and 11% of pairs are at distance 1. The environmental SOS maps share broad climate structure, so they are *not* mostly orthogonal. The expectation in earlier versions of this document, a near-uniform background with a few hot pairs, holds at best for the geographic space.

**`minoverlap`:** Nodiv suggests a floor of 5–10 cells for a geographic grid. At present no co-clustered pair in either space shares fewer than 10 cells, so 3 versus 5–10 changes nothing. It is a per-space setting to revisit if the node sets change. Range overlap does not drive the clustering either: co-clustered pairs share nearly all of the smaller node's SOS cells.

## Sensitivity to the cut

Complete linkage; total clusters = groups of two or more nodes + single nodes.

| cut |r| ≥ | environmental (92 nodes) | geographic (38 nodes) |
|---|---|---|
| 0.6 | 27 = 22 + 5 | 17 = 10 + 7 |
| 0.7 | 35 = 25 + 10 | 18 = 11 + 7 |
| 0.8 | 53 = 22 + 31 | 22 = 10 + 12 |

The geographic answer is stable across the range: about ten groups and a handful of single nodes. The environmental answer depends strongly on the cut, which fits a continuum of variants rather than discrete patterns (next section).

## Why the number of clusters is not chosen from the data

**Silhouette peak.** Partitions for k = 2..40 were scored by mean silhouette (single nodes scored 0) and by bootstrap stability (ARI of `cutree(k)` between the full data and 100 cell resamples). Average linkage beat complete linkage on both scores, in both spaces.
- **Geographic:** a clear optimum at k = 9–11 (silhouette 0.37–0.40, ARI 0.85–0.90).
- **Environmental:** weak structure at every k (silhouette 0.20–0.29; peak k = 8 with ARI 0.70). A silhouette below about 0.25 is conventionally read as no substantial structure.

The silhouette peak always returns some k ≥ 2 and penalises single nodes, so it cannot return the "largely idiosyncratic" result. That makes it relative in the same way as the per-space percentile cut on |r|, which was rejected for imposing the result. A null reference would be needed, but none of the available ones is valid:
- **Shuffling `D` among the pairs** breaks the transitivity that real maps have, and it gives *higher* null silhouettes than observed (environmental 0.29 vs 0.38, geographic 0.40 vs 0.51). It is uninformative, not evidence against structure.
- **Permuting the cells** doesn't work either. A joint permutation of all maps leaves `D` unchanged; independent permutations turn every map into noise with |r| ≈ 0.

**Weighted modularity (check only, not in the repo).** Complete graph with edge weight |r|, standard weighted modularity with the degree null model, and greedy agglomeration (Clauset, Newman & Moore 2004). The null is the same greedy maximisation on 50 graphs with the weights shuffled among the node pairs, which keeps the distribution of |r| but not each node's strength.
- **Environmental:** two communities (55 and 37 nodes), Q = 0.093 against a null mean of 0.040 (maximum of 50 shuffles 0.045). **This weak two-way split is the only structure in either space that beat a null.** It probably corresponds to the two broad moderate-|r| blocks visible in the environmental heatmap (not checked node by node).
- **Geographic:** three communities, Q = 0.100 against a null mean of 0.085 (maximum 0.097). No evidence of structure beyond the null.

**The methods disagree.** The modularity and silhouette-peak partitions agree poorly (ARI 0.18 environmental, 0.27 geographic). At the same k (2 environmental, 3 geographic), average linkage and modularity give an ARI of about −0.02. The number of patterns thus depends entirely on the method. An absolute, stated cut is the honest alternative.

**Complete rather than average linkage** despite average linkage's better scores: complete linkage keeps the guarantee that every pair in a group is above the cut. That guarantee is what makes an absolute cut interpretable. Complete linkage also does not chain marginal pairs into one group.

## Rejected: bootstrap splitting (`supported_clusters`)

The complete-linkage clusters were resampled over the cells, 1000 times with replacement. Any cluster that did not re-form in at least 95% of resamples was split down the dendrogram. This was dropped for two reasons.
- **It can only split, so it only raises the count.** In the environmental space at a cut of 0.7 it ended with 69 clusters: 20 groups and 49 single nodes.
- **The support is only an upper bound.** Resampling cells independently ignores the spatial (and climatic) autocorrelation of neighbouring cells.

Bootstrap support could still be reported as information about individual clusters, but it is not a rule for forming them.

## What to avoid, and why

- **UMAP / t-SNE at this n.** They are unreliable at n ≈ 100 and manufacture clusters from noise. They also build *relative* neighbourhoods, whereas the question needs an absolute notion of "similar enough".
- **Letting the 2-D MDS scatter define groups.** It is a view for browsing only.
- **Ward linkage.** It imposes spherical, Euclidean block structure.
- **A per-space percentile cut, or any rule that fixes the proportion of pairs that count as similar.** It makes the answer relative to each space and imposes the result.

## Two-space interpretation

Treat the two spaces as two questions and do not compare magnitudes between them.
- **Geographic:** do two nodes mark *the same specific break*?
- **Environmental:** do they mark *the same kind of climatic transition*, abstracted from location?

A pair similar in environmental but not geographic space marks the same climatic transition in different places, which is itself a substantive pattern.

## Follow-up

A data-chosen number of patterns would need a valid null: surrogate SOS maps that preserve their autocorrelation, as in spin or shift tests. These don't fit easily on a non-toroidal geographic grid or on the environmental bins. Block resampling of cells would likewise give more honest cluster support than independent resampling. Neither is needed for the current method.
