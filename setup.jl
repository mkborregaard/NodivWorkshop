#!/usr/bin/env julia
#
# setup.jl — one-time setup for a fresh clone or a new machine.
#
# The workshop keeps its large inputs and output figures OUT of git, in a shared
# Google Drive folder that contains `data/` and `figures/`. This script links those
# two folders into the repo and instantiates the Julia environment. Nodiv and the other
# dependencies come from the General registry, except NodivMakie (the plotting), which is
# not registered and is installed from GitHub (see `[sources]` in Project.toml). Nothing
# here is machine-specific.
#
# The optional species images (bow_images/, private and never in the repo) are not set up
# here; the script uses them only if bow_images/workshop_species exists. See README.md.
#
# Usage:
#   julia setup.jl /path/to/NodivWorkshop         # the Drive folder holding data/ & figures/
#   NODIVWORKSHOP_DATA=/path/... julia setup.jl    # ...or via an environment variable
#   julia setup.jl                                 # macOS only: auto-discover the Drive folder
#
# The explicit path (argument or env var) works on any OS. Auto-discovery is a macOS
# convenience only; on Linux/Windows pass the path yourself. On Windows, creating the
# symlinks needs Developer Mode or an elevated shell.

using Pkg

const REPO = @__DIR__

# Bounded, macOS-only probe of Google Drive for desktop. No deep walk: it only looks
# where the shared folder actually lives, so it never traverses the whole Drive.
function autodiscover()
    cs = joinpath(homedir(), "Library", "CloudStorage")
    isdir(cs) || return nothing
    tail = joinpath("EnvSpace_Workshop", "Nodiv project data", "NodivWorkshop")
    for gd in readdir(cs; join = true)
        startswith(basename(gd), "GoogleDrive-") || continue
        stb = joinpath(gd, ".shortcut-targets-by-id")   # shared-drive shortcut targets
        if isdir(stb)
            for id in readdir(stb; join = true)
                cand = joinpath(id, tail)
                isdir(joinpath(cand, "data")) && return cand
            end
        end
        cand = joinpath(gd, "My Drive", tail)  # or directly under My Drive
        isdir(joinpath(cand, "data")) && return cand
    end
    return nothing
end

function resolve_dataroot()
    !isempty(ARGS)                    && return ARGS[1]
    haskey(ENV, "NODIVWORKSHOP_DATA") && return ENV["NODIVWORKSHOP_DATA"]
    d = autodiscover()
    d === nothing && error("""
        Could not locate the shared data folder automatically.
        Pass it explicitly (the folder that contains `data/` and `figures/`):
            julia setup.jl "/path/to/.../Nodiv project data/NodivWorkshop"
        or set the NODIVWORKSHOP_DATA environment variable to that folder.""")
    return d
end

function link_dir(sub, root)
    tgt  = joinpath(root, sub)
    link = joinpath(REPO, sub)
    isdir(tgt) || mkpath(tgt)  # figures/ may not exist yet on a fresh share
    (islink(link) || ispath(link)) && rm(link; force = true, recursive = false)
    try
        symlink(tgt, link)
        println("  linked  $sub  ->  $tgt")
    catch err
        @warn "Could not create a symlink for `$sub`. On Windows this needs Developer Mode " *
              "or an elevated shell; otherwise create the link manually." exception = err
    end
end

root = resolve_dataroot()
isdir(joinpath(root, "data")) ||
    error("No `data/` under $root — is that the right folder?")
println("Data folder: $root")
for sub in ("data", "figures")
    link_dir(sub, root)
end

println("Instantiating the Julia environment (the General registry, and NodivMakie from GitHub)…")
Pkg.activate(REPO)
Pkg.instantiate()
println("Done. Open a REPL with  julia --project=.  and run  include(\"script.jl\")")
