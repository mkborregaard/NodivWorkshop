# NodivWorkshop

Integration workflow that drives [Nodiv.jl](https://github.com/mkborregaard/Nodiv.jl):
`script.jl` runs the node-based SOS/GND analysis end to end and produces the figures.
`preprocess.jl` prepares the inputs.

## Data lives in Google Drive, not in git

The large inputs and the output figures are **not** tracked in this repo. They live in a
shared Google Drive folder that contains two subfolders:

```
<Google Drive>/…/EnvSpace_Workshop/Nodiv project data/NodivWorkshop/
├── data/      # inputs + the cached node_analysis.jld2 (res_e / res_g)
└── figures/   # outputs
```

In a working clone, `data/` and `figures/` are **symlinks** into that Drive folder
(both are gitignored). The expensive `node_metrics(...; nsims=…)` run is cached to
`data/node_analysis.jld2` — don't recompute it; the script loads it from there.

## Setup on a new machine

1. Clone this repo and make sure the shared Drive folder is synced locally (you need
   access to it).
2. From the repo, run:

   ```bash
   julia setup.jl "/path/to/.../Nodiv project data/NodivWorkshop"
   ```

   That path is the Drive folder holding `data/` and `figures/`. On macOS you can usually
   omit it — `setup.jl` auto-discovers the Google Drive folder. You can also set
   `NODIVWORKSHOP_DATA` instead of passing the argument.

   `setup.jl` symlinks `data/` and `figures/` into the repo and instantiates the Julia
   environment (Nodiv and all dependencies come from the General registry).

3. Run the analysis:

   ```bash
   julia --project=. script.jl
   ```

> On Windows, creating the symlinks needs Developer Mode or an elevated shell.

## Depending on Nodiv

The workshop depends on the **registered** `Nodiv` (from the General registry) — no local
paths, so a clone is reproducible anywhere. To test local, unreleased changes to Nodiv,
dev it into this environment on your own machine:

```julia
julia --project=. -e 'using Pkg; Pkg.develop(path="/path/to/Nodiv")'
```

and `Pkg.free("Nodiv")` to switch back to the registered version. (Keep the dev override
out of committed `Project.toml`/`Manifest.toml`.)
