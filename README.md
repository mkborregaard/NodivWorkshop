# NodivWorkshop

An analysis in progress with [Nodiv.jl](https://github.com/mkborregaard/Nodiv.jl): the
node-based SOS/GND analysis of birds in environmental and geographic space.
`preprocess.jl` prepares the inputs. `script.jl` is worked through interactively, a line or
block at a time and not necessarily in order; it is not meant to be run top to bottom.
`explorer.jl` opens the interactive node explorers straight from a terminal.
Helper functions that are finished live in `functions.jl`, which the script `include`s;
functions still being developed stay in the script. Once finished, general analysis code
moves to Nodiv and general plotting code to NodivMakie; only what is specific to this
data (e.g. the traits and the taxon names of nodes) stays here. All plotting is done with
[Makie](https://docs.makie.org), through
[NodivMakie.jl](https://github.com/mkborregaard/NodivMakie.jl) for the trees, maps, node
panels and cluster figures, with GLMakie for interactive windows.
`docs/sos_pattern_grouping_design.md` is the design spec for grouping the divergent nodes
by SOS-pattern similarity.

## Data lives in Google Drive, not in git

The large inputs and the output figures are **not** tracked in this repo. They live in a
shared Google Drive folder that contains two subfolders:

```
<Google Drive>/…/EnvSpace_Workshop/Nodiv project data/NodivWorkshop/
├── data/      # raw RDS, cleaned inputs + the cached node_analysis.jld2 (res_e / res_g)
└── figures/   # outputs
```

In a working clone, `data/` and `figures/` are **symlinks** into that Drive folder
(both are gitignored). The expensive `node_metrics(...; nsims=…)` run is cached to
`data/node_analysis.jld2`. Don't recompute it: the script loads it from there.

`preprocess.jl` reads the matched raw data, `data/data_birds_matched_simplified.rds`
(phylogeny, presences and grids for both spaces, and AVONET traits), through RCall, so it
needs R with the `sf` and `ape` packages. It writes the cleaned CSVs and the pruned tree
that `script.jl` loads. The geographic grid is a Behrmann equal-area grid, so its cells are
placed on an exact regular grid in Behrmann coordinates (km).

## Setup on a new machine

1. Clone this repo and make sure the shared Drive folder is synced locally (you need
   access to it).
2. From the repo, run:

   ```bash
   julia setup.jl "/path/to/.../Nodiv project data/NodivWorkshop"
   ```

   That path is the Drive folder holding `data/` and `figures/`. On macOS you can usually
   omit it: `setup.jl` auto-discovers the Google Drive folder. You can also set
   `NODIVWORKSHOP_DATA` instead of passing the argument.

   `setup.jl` symlinks `data/` and `figures/` into the repo and instantiates the Julia
   environment.

3. Run the analysis in a REPL (or step through it in VS Code), so the figures can be shown:

   ```bash
   julia --project=.
   ```

   and then evaluate `script.jl` line by line or block by block (e.g. in VS Code). The
   `using` block and `include("functions.jl")` at the top come first.

> On Windows, creating the symlinks needs Developer Mode or an elevated shell.

## Figures

Every figure is kept in a variable (`richness_g`, `metric_tree_e`, `heat_g`,
`tree_clusters_e`, …). Evaluate one to show it. Below each figure is a commented-out
`save("figures/…", fig)` line: uncomment it to write that figure. GLMakie writes raster
formats; for vector files use `save("figures/name.pdf", fig; backend = CairoMakie)`. The
(also commented-out) `node_panel_pdf` calls at the end write one multi-page PDF of node
panels per space, which needs `pdfunite` (poppler, `brew install poppler`).

The interactive entry point is the **node explorer**, one per space
(`explorer_fig_e`, `explorer_fig_g`). In a REPL they open in their own windows. Each shows
the fan tree with the divergent nodes marked, next to the node's SOS map, the richness of
its two child clades, and an ordination (MDS) of the divergent nodes by the similarity of
their SOS maps.
- Click a node marker, or the branch leading to a node, to show that node. Clicking a
  point in the ordination does the same; the node shown has a ring there.
- Hover over nodes, branches, ordination points and map cells for labels.
- The two explorers are linked: a node picked in one is shown in the other, if it has
  an SOS there.
- `explorer_g.panel.node[] = "Node 17672"` shows a node from code.

Set `NODIVWORKSHOP_WINDOWS=false` to keep the explorers from opening windows (e.g. when
running the script headless).

To use just the explorers, without working through `script.jl`, run `explorer.jl` from a
terminal:

```bash
julia explorer.jl
```

It needs the setup above and the cached `data/node_analysis.jld2`, which it loads (it
never recomputes the analysis). It opens the two linked explorers, environmental and
geographic, each in its own window, and ends when both windows are closed. Getting to the
windows takes a few minutes, mostly loading and compiling the packages.

Hovering over a node in an explorer shows the families below it (or, for a single family,
its genera), and hovering over a species image shows the species' family and order.

### Species images (optional, not in the repo)

The explorers can show Birds of the World illustrations: around the tree, one per clade,
and in the corner of the child-clade maps. The images are copyrighted, so they are **not
in this repo** and must never be committed. If you have them, put them (or a symlink to
them) in `bow_images/workshop_species/`, one file per species named by the tree's
species name (e.g. `Carduelis_hornemanni.jpg`). `bow_images/` carries its own
`.gitignore` that ignores everything in it. The script uses the images only if that
folder exists, and runs the same without them.

## Dependencies

Nodiv and the other dependencies come from the General registry. NodivMakie is not
registered: it is installed from GitHub via `[sources]` in `Project.toml`, and
`Manifest.toml` pins the commit. To move to the latest NodivMakie, run `Pkg.update("NodivMakie")`.

To test local, unreleased changes to Nodiv or NodivMakie, dev them into this
environment on your own machine:

```julia
julia --project=. -e 'using Pkg; Pkg.develop(path="/path/to/Nodiv")'
```

and `Pkg.free("Nodiv")` to switch back. (Keep the dev override out of the committed
`Project.toml`/`Manifest.toml`.)

After updating NodivMakie, run the tests of `functions.jl`. The explorers' hover labels
for the species images use NodivMakie internals, and the tests catch an update that
breaks them:

```bash
julia --project=. test/runtests.jl
```
