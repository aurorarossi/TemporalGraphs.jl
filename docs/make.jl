using Documenter
using TemporalGraphs
using GNNGraphs  # load (and precompile) the extension before the doctests

DocMeta.setdocmeta!(TemporalGraphs, :DocTestSetup, :(using TemporalGraphs); recursive=true)

makedocs(;
    modules=[TemporalGraphs],
    sitename="TemporalGraphs.jl",
    authors="Aurora Rossi",
    format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true", edit_link="main",
                           size_threshold_warn=250 * 2^10, size_threshold=500 * 2^10),
    pages=[
        "Home" => "index.md",
        "Getting started" => "tutorial.md",
        "Temporal paths and distances" => "distances.md",
        "Centralities and statistics" => "centralities.md",
        "Motifs and reference models" => "motifs.md",
        "Connectivity, flows, spanners and isomorphisms" => "theory.md",
        "Graph representations" => "representations.md",
        "Datasets and machine learning" => "datasets.md",
        "Performance tips" => "performance.md",
        "API reference" => "api.md",
    ],
    checkdocs=:exports,
    doctest=true,
)

if get(ENV, "CI", "false") == "true"
    deploydocs(; repo="github.com/aurorarossi/TemporalGraphs.jl", devbranch="main")
end
