#!/usr/bin/env julia
#
# setup.jl — set up a fresh clone, and bring an existing one up to date.
#
# The large inputs, the output figures and the species images are kept OUT of git, in the
# shared Google Drive folder `Nodiv/`:
#
#     Nodiv/
#     ├── NodivWorkshop/data/      inputs and caches (objects.jld2, node_analysis.jld2)
#     ├── NodivWorkshop/figures/   outputs
#     └── bow_images/              Birds of the World illustrations (private)
#
# This script links `data/` and `figures/` into the repo, builds `bow_images/` (a link to
# the images and one link per workshop species, named by the tree's names, from
# bow_images/crosswalk.csv), and installs the Julia packages pinned in Manifest.toml. It is
# safe to re-run, and re-running it is how a clone is updated after a `git pull`.
#
# Usage:
#   julia setup.jl                                 # find the Drive folder (macOS, Windows)
#   julia setup.jl "/path/to/Nodiv/NodivWorkshop"  # or pass it
#   NODIVWORKSHOP_DATA=/path/... julia setup.jl    # or set it in the environment
#
# After the first run the folder is remembered through the `data` link. On Windows,
# creating the links needs Developer Mode (Settings > System > For developers).

using Pkg

const REPO = @__DIR__
const TAIL = joinpath("Nodiv", "NodivWorkshop")

# Where Google Drive for desktop mounts: ~/Library/CloudStorage/GoogleDrive-* on macOS,
# a drive letter (G: by default) on Windows.
function drive_roots()
    if Sys.isapple()
        cs = joinpath(homedir(), "Library", "CloudStorage")
        isdir(cs) || return String[]
        return filter(d -> startswith(basename(d), "GoogleDrive-"), readdir(cs; join=true))
    elseif Sys.iswindows()
        return filter(isdir, ["$(l):\\" for l in 'D':'Z'])
    end
    return String[]
end

# Bounded probe: only where a shared folder added to My Drive can appear, never a deep walk.
function auto_discover()
    for gd in drive_roots()
        cand = joinpath(gd, "My Drive", TAIL)
        isdir(joinpath(cand, "data")) && return cand
        stb = joinpath(gd, ".shortcut-targets-by-id")
        isdir(stb) || continue
        for id in readdir(stb; join=true)
            cand = joinpath(id, TAIL)
            isdir(joinpath(cand, "data")) && return cand
        end
    end
    return nothing
end

function resolve_data_root()
    !isempty(ARGS) && return ARGS[1]
    haskey(ENV, "NODIVWORKSHOP_DATA") && return ENV["NODIVWORKSHOP_DATA"]
    link = joinpath(REPO, "data")
    islink(link) && return dirname(readlink(link))  # Set up before
    d = auto_discover()
    d === nothing && error("""
        Could not find the shared folder Nodiv/NodivWorkshop in Google Drive.
        Check that Google Drive for desktop is running and that the shared `Nodiv` folder
        is added to your My Drive, or pass the folder (the one containing `data/`):
            julia setup.jl "/path/to/Nodiv/NodivWorkshop\"""")
    return d
end

# Point `link` at `target`, replacing an earlier link but never a real folder.
function link_dir(target, link)
    if islink(link)
        rm(link)
    elseif ispath(link)
        println("  kept    $(relpath(link, REPO))  (a real folder, not a link)")
        return false
    end
    symlink(target, link; dir_target=true)
    println("  linked  $(relpath(link, REPO))  ->  $target")
    return true
end

function link_or_warn(target, link)
    try
        return link_dir(target, link)
    catch err
        msg = "Could not link $(relpath(link, REPO)). On Windows this needs Developer Mode."
        @warn msg exception = err
        return false
    end
end

# A small CSV reader (the packages are not installed yet): rows as Dicts by column name.
# A quoted field may hold commas.
function readcsv(path)
    fields(l) =
        [strip(m.captures[1], '"') for m in eachmatch(r"(\"[^\"]*\"|[^,]*)(?:,|$)", l)]
    lines = filter(!isempty, readlines(path))
    header = fields(first(lines))
    return [Dict(zip(header, fields(l))) for l in lines[2:end]]
end

# The species images, as `bow_images/workshop_species/<tree name>.jpg` links into the
# shared images, built from the shared crosswalk (written by bow_images/crosswalk.jl).
function link_images(src)
    dir = joinpath(REPO, "bow_images")
    if !isdir(src)
        println("  no species images at $src; the explorers run without them")
        return nothing
    end
    if isdir(joinpath(dir, "320")) && !islink(joinpath(dir, "320"))
        println("  kept    bow_images  (a local copy of the images)")
        return nothing
    end
    mkpath(dir)
    link_or_warn(joinpath(src, "320"), joinpath(dir, "320")) || return nothing
    key(s) = lowercase(replace(strip(s), r"[\s_\-]+" => "_"))  # As in crosswalk.jl
    manifest = readcsv(joinpath(src, "manifest.csv"))
    asset = Dict(key(r["sciName"]) => r["assetId"] for r in manifest)
    species = joinpath(dir, "workshop_species")
    rm(species; recursive=true, force=true)
    mkpath(species)
    n = 0
    for r in readcsv(joinpath(src, "crosswalk.csv"))
        id = get(asset, key(r["bow_name"]), nothing)
        id === nothing && continue
        name = replace(strip(r["workshop_name"]), r"\s+" => "_") * ".jpg"
        symlink(joinpath("..", "320", id * ".jpg"), joinpath(species, name))
        n += 1
    end
    println("  linked  bow_images/workshop_species  ($n species images)")
    return nothing
end

root = rstrip(normpath(resolve_data_root()), ['/', '\\'])
isdir(joinpath(root, "data")) || error("No `data/` in $root: is that the right folder?")
println("Shared folder: $root")
for sub in ("data", "figures")
    isdir(joinpath(root, sub)) || mkpath(joinpath(root, sub))
    link_or_warn(joinpath(root, sub), joinpath(REPO, sub))
end
link_images(joinpath(dirname(root), "bow_images"))

println("Installing the Julia packages pinned in Manifest.toml…")
Pkg.activate(REPO; io=devnull)
Pkg.instantiate()
println("Done. Open the explorers with  julia explorer.jl")
