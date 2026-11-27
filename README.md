# NodivWorkshop

The node-based SOS/GND analysis of birds in environmental and geographic space, with
[Nodiv.jl](https://github.com/mkborregaard/Nodiv.jl) and
[NodivMakie.jl](https://github.com/mkborregaard/NodivMakie.jl).

## Setup

You need access to the shared Google Drive folder `Nodiv`. It holds the data, the figures
and the species images, none of which are in this repo.

1. **Google Drive.** Install [Google Drive for desktop](https://www.google.com/drive/download/)
   and sign in. At [drive.google.com](https://drive.google.com), open *Shared with me*,
   right-click `Nodiv` and choose *Organize → Add shortcut → My Drive*. That makes the
   folder appear on your disk.
2. **Julia.** Install it with [juliaup](https://github.com/JuliaLang/juliaup):
   - macOS: `curl -fsSL https://install.julialang.org | sh` (accept the defaults)
   - Windows: `winget install julia -s msstore`

   Open a new terminal and check that `julia --version` works. If macOS says
   `command not found: julia`, add juliaup to your PATH and open a new terminal again:

   ```bash
   echo 'export PATH="$HOME/.juliaup/bin:$PATH"' | tee -a ~/.zprofile ~/.zshrc
   ```

3. **The repo.** Clone it. On macOS, if `git` is missing, the first `git` command offers
   to install it. On Windows, run `winget install Git.Git` first.

   ```bash
   git clone https://github.com/mkborregaard/NodivWorkshop.git
   cd NodivWorkshop
   ```

4. **Set up.** This links `data/`, `figures/` and `bow_images/` from the Drive folder into
   the repo and installs the Julia packages:

   ```bash
   julia setup.jl
   ```

   If it can't find the Drive folder, pass it: `julia setup.jl "/path/to/Nodiv/NodivWorkshop"`.
   On Windows, the links need Developer Mode (*Settings → System → For developers*). A
   message that RCall found no R installation is harmless: R is only needed to rerun
   `preprocess.jl`.
5. **Run the explorers.**

   ```bash
   julia explorer.jl
   ```

   The first run takes several minutes: Julia compiles the packages and Drive downloads
   about 400 MB of data. See [The explorers](#the-explorers) for what you can do in them.

## Updating

When the code, the package versions or the image list change, run this from the
`NodivWorkshop` folder:

```bash
git pull
```

```bash
julia setup.jl
```

`setup.jl` remembers the Drive folder and is safe to rerun. If you have edited
files in the repo yourself, `git pull` may refuse to run. Use `git stash` to set your
edits aside and `git stash pop` to get them back afterwards. If the packages need a newer
Julia, run `juliaup update`.

New data and figures in the Drive folder reach you through Drive. You don't need to do
anything.

**For the maintainer.** What each kind of change needs before people can get it:

| Change | Before `git pull` and `julia setup.jl` pick it up |
|---|---|
| Workshop scripts | Commit and push. |
| Nodiv or NodivMakie | Release the package. In the workshop, run `Pkg.update("Nodiv")`, or for NodivMakie `Pkg.add(url="https://github.com/mkborregaard/NodivMakie.jl", rev="v<new tag>")`. Run the tests, then commit `Project.toml` and `Manifest.toml`. |
| Species image list | Rerun `crosswalk.jl` in your local `bow_images/`, then copy `crosswalk.csv` and any new images in `320/` to `Nodiv/bow_images` in Drive. Nothing to commit. |
| `objects.jld2`, `node_analysis.jld2`, figures | Replace the files in the Drive folder. People get them without running anything. |

## What's here

| File | What it is |
|---|---|
| `explorer.jl` | Opens the two linked node explorers (environmental and geographic space) from a terminal. |
| `traitsexplorer.jl` | Opens the trait explorer from a terminal. |
| `script.jl` | The main analysis. Work through it interactively, a line or a block at a time. It is not meant to be run top to bottom. |
| `traitsscript.jl` | The analysis of the clades' traits: the PCA of the AVONET morphometrics, the trait explorer, and the overlap of the child clades in trait space. Worked through the same way. |
| `functions.jl`, `traitfunctions.jl` | Finished helper functions for the scripts. |
| `preprocess.jl` | Cleans the raw data into `data/clean/`. |
| `create_objects.jl` | Builds the tree and the assemblages from the cleaned data. |
| `setup.jl` | Setup and updates (see above). |
| `docs/sos_pattern_grouping_design.md` | The design spec for grouping divergent nodes by how similar their SOS patterns are. |
| `docs/data_migration.md` | Why the data and the tree are what they are. |

Helper functions that are still being developed stay in the scripts. Once they are
finished, general analysis code moves to Nodiv and general plotting code to NodivMakie.
Only what is specific to this data, such as the traits and the taxon names of nodes,
stays here.

## The explorers

`julia explorer.jl` opens one window per space. Each shows the fan tree with the
divergent nodes marked, the SOS map of the node shown, the richness of its two child
clades, and an ordination (MDS) of the divergent nodes by how similar their SOS maps are,
coloured by SOS cluster (grey for nodes on their own). The ordination is a 2-D projection
for browsing. The clusters come from the full distances.

- Click a node, the branch leading to it, or its point in the ordination to show it.
- Hover over nodes, branches, ordination points, map cells and species images for labels.
  A node's label lists the families below it, or the genera for a single family.
- The two windows are linked: a node picked in one is shown in the other, if it has an
  SOS there.
- The program ends when both windows are closed.

`julia traitsexplorer.jl` opens one window with the tree, the SOS maps, and the two child
clades of the node shown in trait space (PCA axes 3–4, beak shape and locomotor shape,
with contours around 95%, 50% and 25% of their trait densities; `TRAIT_PCS` in
`traitfunctions.jl` switches to axes 2–3, relative beak size and beak shape). It also has a scatter of every node's environmental rms
against its geographic rms, coloured by the overlap of those densities. Click a node in
the tree or the scatter to show it.

The species images are Birds of the World illustrations. They are copyrighted, so they
must never be committed: `bow_images/` is gitignored. Without them the explorers work the
same, only without the pictures.

## Working through `script.jl`

Open a REPL with `julia --project=.` (or use VS Code) and evaluate the `using` block and
`include("functions.jl")` at the top first. After that, go through the script in any order.

Every figure is kept in a variable (`richness_g`, `metric_tree_e`, `heat_g`,
`tree_clusters_e`, …). Evaluate the variable to show the figure. Below each figure is a
commented-out `save("figures/…", fig)` line. Uncomment it to write the figure to Drive.
GLMakie writes raster formats. For vector files, use
`save("figures/name.pdf", fig; backend = CairoMakie)`. The `node_panel_pdf` calls at the
end write one multi-page PDF of node panels per space. They need `pdfunite`
(`brew install poppler`).

The explorers are `explorer_fig_e` and `explorer_fig_g`. `explorer_g.panel.node[] = "Node 17672"`
shows a node from code. Set `NODIVWORKSHOP_WINDOWS=false` to keep them from opening
windows, for example when running headless.

## The data pipeline

The expensive step, `node_metrics(...; nsims=…)`, is cached in `data/node_analysis.jld2`
(`res_e`, `res_g`). The scripts load it from there. Don't recompute it.

1. `preprocess.jl` reads `data/Bird_data_PCA30.rds` through RCall, so it needs R with
   the `sf` and `ape` packages. The file holds the presences and grids for both spaces,
   three candidate phylogenies, the AVONET traits with a morphometric PCA, and the
   loadings of the environmental PCA. The script writes cleaned CSVs and the pruned tree
   to `data/clean/`. It uses the `McTavish_2025` phylogeny, resolves its few polytomies at
   random under a fixed seed, and makes it ultrametric. The geographic grid is a Behrmann
   equal-area grid, with the cells placed exactly on it (in km). The earlier
   `data_birds_matched_simplified.rds` is still read, if present, for the raw CHELSA /
   MODIS / GLWD / elevation layers, which the new file keeps only as their PCA.
2. `create_objects.jl` builds the tree and the two assemblages, `birds_e` and `birds_g`,
   with their site covariates, traits and the morphometric PCA axes `pca1 … pca9`. It
   caches them in `data/objects.jld2`, together with the variance on those axes as
   `trait_pca`. Rerun it whenever `preprocess.jl` has been rerun:

   ```bash
   julia --project=. create_objects.jl
   ```

## Developing Nodiv or NodivMakie

The workshop uses released versions: Nodiv from the General registry and NodivMakie from
GitHub at the tag in `[sources]` in `Project.toml`. To try unreleased changes, dev them
into the environment on your own machine:

```bash
julia --project=. -e 'using Pkg; Pkg.develop(path="/path/to/Nodiv")'
```

Run `Pkg.free("Nodiv")` to switch back. Never commit the dev override in `Project.toml`
or `Manifest.toml`.

After updating NodivMakie, run the tests. The hover labels on the species images use
NodivMakie internals, and the tests catch an update that breaks them:

```bash
julia --project=. test/runtests.jl
```
