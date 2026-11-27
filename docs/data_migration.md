# Migrating to `Bird_data_PCA30.rds`

Status: **executed 2026-11-24** on the decisions below — `McTavish_2025` with its
polytomies resolved stochastically, the supplied trait PCA in place of the one the trait
scripts used to fit, and the raw geographic layers carried over. Record of what the new
raw data contains, how it differs from what the analysis used before, and what was done.

## The two raw files

| | old `data_birds_matched_simplified.rds` | new `Bird_data_PCA30.rds` |
|---|---|---|
| geo grid | `g_space$grid_sf`, 17,648 cells, 49 cols | `g_space$g_space_sf`, 17,648 cells, 8 cols |
| geo PAM | `g_space$presence`, 2,993,566 rows | `g_space$g_space_PAM`, 3,003,599 rows |
| env grid | `e_space$grid_sf`, 900 cells | `e_space$e_space_sf`, 900 cells |
| env PAM | `e_space$presence`, 475,555 rows | `e_space$e_space_PAM`, 476,782 rows |
| phylogeny | one tree, 10,324 tips | `Jetz_2012`, `McTavish_2025`, `Claramund_2025`, 10,354 tips each |
| traits | `traits`, 10,316 x 37 | `traits$avonet_raw_traits` 10,354 x 38, plus `PCA_scores` (10,354 x 9) and `PCA_loadings` |
| PCA metadata | — | `env_PCA_loadings` (17 vars x 17 PCs) |

Field names changed throughout (`grid_sf` -> `g_space_sf`, `presence` -> `g_space_PAM`,
and so on), so `preprocess.jl` needs updating regardless of which parts are adopted.

## What actually differs

**Geography: identical.** Same 17,648 cells, same `ID_geo` set *and order*, centroids
agreeing to 5e-13 degrees, identical vertex count (492,473), `area_km` to 8.9e-16.
`PC1-PC3` agree at r = 0.999999 (max diff 0.014).

The new grid **drops 41 columns**: `area_m` and every raw CHELSA / MODIS / GLWD /
elevation layer. Nothing in the workshop code reads them (only a comment in
`create_objects.jl:53` mentions them), but they exist *only* in the old file. Since cell
order is identical they can be carried across by a merge on `ID_geo` if wanted.

**Environment: one cell explains the whole change.** Cell `ID_geo` 28903 had
`PC1` = 10.009, the single largest of 17,509 values, and is `NA` in the new file. That
drops the PC1 range to 7.394, so the 30-bin width goes 0.5478 -> 0.4607 and **17,223 of
17,648 cells land in a different `ID_env`** (only 794 keep `bin_PC1`). The 30x30 scheme
and the `e_1_1 ... e_30_30` labels are unchanged.

The new file also fixes `occupied` (489 bins, matching the data; the old `occupid` said
407, which was simply wrong), fills `area_km` for all 900 bins (old had 493 `NA`s), adds
`n_cells`, and drops the bin edges `xmin/xmax/ymin/ymax`.

**The old env space was self-inconsistent.** Deriving the env PAM from each file's own
geo->env mapping: the new file matches its stated PAM exactly (Jaccard 1.000); the old
file does not (Jaccard 0.413). The env presences are *identical* between the files for
shared species, so the presences already corresponded to the new binning and it was the
old grid's `PC`/`bin`/`ID_env` columns that were stale. Consequence: the current
`coords_e.csv` is wrong relative to the current `phylocom_e.csv`, displaced by up to
2.57 PC1 units. Both midpoint sets are monotone in bin index, so existing env maps are a
stretched version of the truth rather than scrambled.

**Taxonomy, not new sampling, drives the species differences.** 291 new names, 256
retired (`Amazilia` -> `Saucerottia`/`Chrysuronia`/`Leucolia`, `Alopecoenas` ->
`Pampusana`, `Arborophila` -> `Tropicoperdix`). Among the 10,063 shared species the geo
PAM gained 4,359 records and lost none. The new `AVONET_species` column carries the old
names and bridges the change (10,295 of 10,354 rows map onto the old tree's tips).

**Known warts in the new file.** 218 `ID_geo` values in the geo PAM do not exist in its
grid; they carry 4,417 of 3,003,599 records (0.15%) and dropping them loses **no**
species, leaving 17,542 assemblage sites — the current count. 140 cells have no
`ID_env`, 104 of them occupied. `mean_geo_dist_km` is `NA` for 451 of 900 bins.

**Latent bug in the current clean data**, unrelated to the new file: `bin_PC1`/`bin_PC2`
in `sitestats_g.csv` hold `-2147483648` for the 139 cells without PC values (R's integer
`NA` leaking through RCall as a Float64 sentinel) instead of `missing`.

## The tree

The tree currently in use is **Jetz 2012**, not McTavish as assumed: RF 484 (normalized
0.024) against `Jetz_2012`, versus 13,000 (0.649) for `McTavish_2025` and 12,285 (0.667)
for `Claramund_2025`, with the root depth matching to all eight significant figures. The
residual difference is resolution — the current tree has 26 polytomies against Jetz's
168 — so it is a better-resolved product of the same source, carrying a 30.47 Myr
tip-depth spread and 1,778 zero-length branches.

Decision metric: `Nodiv._analysis_prepass` (`src/functions.jl:281`) only analyses a node
with **exactly two children**, each with >= 3 species. Polytomies are not analysed at all.

| on the 10,354 species in both PAMs | current (Jetz-like) | Jetz_2012 | McTavish_2025 | Claramund_2025 |
|---|---|---|---|---|
| tips kept | 10,063 | 10,354 | 10,354 | 10,354 |
| internal nodes | 10,027 | 10,062 | 10,290 | 8,669 |
| polytomies (max children) | 26 (6) | 168 (12) | 49 (9) | 551 (39) |
| **analysable nodes** | 1,907 | 1,980 | **1,909** | **1,526** |
| after random resolution | 1,930 | 2,000 | 1,914 (+5) | 1,826 (+300) |
| root depth | 120.110 | 120.110 | 116.880 | 107.183 |
| tip-depth spread | 30.47 | 5e-8 | 1e-5 | 3e-8 |
| zero-length branches | 1,778 | 0 | 0 | 0 |
| monophyletic genera | 1,004/1,349 | 993/1,393 | 1,186/1,393 | 1,186/1,393 |
| monophyletic families | 154/190 | 154/190 | 180/190 | 181/190 |
| monophyletic orders | 29/29 | 29/29 | 29/29 | 28/29 |

McTavish and Claramunt differ in 32% of splits (RF 6,037); they share 6,461 clades,
62.8% of McTavish's and 74.5% of Claramunt's.

### Resolving polytomies stochastically

Polytomies can be resolved with `multi2di`, which makes them analysable at the cost of
inventing topology. Measuring that cost over 5 seeds, counting analysable nodes whose
clade is absent from the source tree:

| | McTavish_2025 | Claramund_2025 |
|---|---|---|
| nodes resolution must invent | 63 (**0.6%** of the tree) | 1,684 (**16.3%** of the tree) |
| analysable, unresolved | 1,909 | 1,526 |
| analysable, resolved | 1,914 | 1,824 |
| of those, arbitrary | 2.8 (**0.1%**) | 168.4 (**9.2%**) |
| source-supported analysable nodes | **1,911** | **1,655** |

For McTavish the operation is free: 49 polytomies, +5 nodes, one node in a thousand
arbitrary. For Claramunt it is the whole difference: a sixth of the tree invented and
roughly one in eleven analysed clades arbitrary, each then carrying GND/SOS values and
feeding the SOS-pattern grouping as though it were real.

**Recommendation: McTavish_2025, with its 49 polytomies resolved stochastically under a
fixed seed.** It yields 1,911 source-supported analysable nodes against Claramunt's
1,655 — and that gap survives resolution rather than being closed by it. The two are tied
on monophyletic genera; Claramunt's one extra monophyletic family is offset by one lost
order. Claramunt is already ultrametric, but McTavish's 1e-5 deviation is 8.6e-8 relative
— a rounding artifact, fixed by extending tip branches. Write the resolved tree to
`clean/tree.nwk`, which makes the result reproducible as an artifact regardless of seed.

Either choice gains 291 species over the current analysis and removes the 30 Myr
tip-depth spread, and raises monophyletic genera from 74% to 85% and families from 81% to
95%, so many more clades carry taxon names in the figures.

## What was done

`preprocess.jl` was rewritten for the new layout and run; `create_objects.jl` rebuilt
`objects.jld2`; the node analysis was recomputed from scratch. Reported by the run:

| | |
|---|---|
| tree | `McTavish_2025`, 10,354 tips, 49 polytomies resolved, 63 nodes invented |
| ultrametric | tip-depth spread 1.0e-5 -> 2.8e-14, root depth 116.880 |
| geographic presences | 2,999,182 kept; 4,417 dropped on 218 sites absent from the grid |
| environmental presences | 476,782 |
| taxon-named internal nodes | 1,342 (was ~1,187) |
| raw layers carried over | 41 columns, matched by `ID_geo` |

Verified afterwards: the sentinel is gone and `bin_PC1`/`bin_PC2` are missing exactly
where `PC1` is; `sitestats_g.csv` has its 51 columns with all 34 CHELSA layers;
`coords_e.csv` now equals `sitestats_e.csv`'s bin midpoints, so the coordinate mismatch
described above is fixed; the env PAM derived through `ID_env` equals the stated one
(Jaccard 1.000); species sets agree across `phylocom_e/g.csv`, `traits.csv`,
`trait_pca_scores.csv` and `tree.nwk` at 10,354; the tree is binary, ultrametric and has
10,353 uniquely named internal nodes; the PCA proportions sum to 1 and match the scores.
Re-running `preprocess.jl` reproduced its log byte-identically, so the seeded tree is
deterministic.

Kept in case the earlier state is wanted: `data/clean_pre_migration/` and
`data/node_analysis_jetz2012.jld2`.

Unrelated pre-existing breakage, found while verifying: 5 failures in
`test/runtests.jl`, all because commit 2f9b44a ("Show the order(s) of a clade in the node
hover labels") made `taxa_text` emit an order line while the tests, last updated in the
older e20fdac, still expected none. The behaviour matches the function's docstring, so
the tests were brought up to date rather than the behaviour reverted, and a case for the
`norders` cut-off was added. The suite passes, 38 tests.

## The species images

`bow_images/` is private and git-ignores its own contents, so it was outside the
migration — but the new tree invalidated its name crosswalk, and the folder had gone
missing from the repo. It was restored from `My Drive/Nodiv/bow_images`: 11,055 JPEGs in
`320/`, named by Birds of the World asset id, plus `manifest.csv` and `crosswalk.jl`.
`workshop_species/`, the folder the scripts actually look for, is not stored anywhere —
`crosswalk.jl` derives it as symlinks into `320/` named by the tree's species.

Rebuilding it needed a fix first: `crosswalk.jl`'s own Newick reader stripped whitespace
from each label but not quotes, and Phylo.jl's writer quotes every label, so each name
arrived as `"Struthio_camelus"` and nothing matched. `readlabel` now strips both. The
pre-migration `tree.nwk` is quoted too, so this was never a regression; the September run
must have read a differently written tree.

Rebuilt against `McTavish_2025`:

| | |
|---|---|
| exact name match | 9,553 |
| matched via context | 424 |
| matched via genus | 4 |
| linked | 9,981 of 10,354 (96.4%) |
| no BOW match | 373 |

Reusing the Jetz-era crosswalk would have matched only 7,381, since McTavish and Birds of
the World are both eBird/Clements-based while BirdTree is not — the fuzzy matching now
only has to cover the residue. The 373 are species Birds of the World lumps differently,
concentrated in Thraupidae, Columbidae and Picidae. A further 8 links point at assets that
failed in the original download (`failed_assets.txt`), so 9,973 species have an image.

Verified by rebuilding the exact-match set independently from `manifest.csv` and
comparing: all 9,553 present, no asset-id mismatch, no link naming a species absent from
the tree. Every copied file matches Drive byte for byte; one arrived empty and was
re-copied. Drive keeps its earlier outputs as `crosswalk_jetz2012.csv` and
`crosswalk_review_jetz2012.txt`.

## Steps taken

1. **Decide**: tree (McTavish recommended); whether to carry the 41 raw geo columns over.
2. **Rewrite `preprocess.jl`** for the new layout: new field names; drop the `mn_g_d_` /
   `occupid` renames and the `area_m` parse (all gone upstream); no tip-label
   underscoring or `node.label` stripping needed; select the tree by name; make it
   ultrametric by extension, asserting the spread before and after; resolve polytomies
   with `multi2di` under a fixed seed; drop the 4,417 PAM records on unmapped sites;
   write `bin_PC1`/`bin_PC2` as `missing`, not the sentinel.
   Keep the Behrmann col/row derivation — the geometry is identical.
3. **Regenerate `data/clean/`**: the existing eight files plus `env_pca_loadings.csv`,
   `trait_pca_loadings.csv` and `trait_pca_scores.csv`. Expect 10,354 species, 17,542 geo
   sites, 489 occupied env bins.
4. **Rebuild `objects.jld2`** via `create_objects.jl` (expected to need only its comment
   updated, since columns are keyed by `:ID_geo` / `:ID_env`).
5. **Recompute `node_metrics` for both spaces** and rewrite `node_analysis.jld2`. The
   user authorised this one-time long run on 2026-11-24, overriding the CLAUDE.md
   do-not-recompute rule.
6. **Repoint `script.jl`**: `focal = "Node 16367"` at lines 131 and 288 refers to a node
   that will not exist; all auto-generated node names change with the tree.
7. **Verify**: derived env PAM equals the stated one; site and species counts as above;
   `Pkg.test()`.

Note for step 3: the workshop's own `trait_pca` (`traitfunctions.jl:31`) spans
`Beak_Length_Culmen` to `Mass`, 11 columns including `Kipps_Distance` and the
hand-wing index, which are derived from the others. The new file's `PCA_loadings` uses
nine non-redundant traits instead. Adopting the supplied scores is a separate scientific
decision; shipping them to `clean/` as data changes nothing on its own. The trait column
`Hand-Wing_Index` is also renamed `Hand_Wing_Index`.
