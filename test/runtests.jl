using TemporalGraphs
using Random
using Test
using Aqua
using JET
import Graphs
import GNNGraphs

include("instances.jl")
include("oracle.jl")

const DATA = joinpath(@__DIR__, "data")
paper_file() = load_ordered_edge_list(joinpath(DATA, "example_from_paper.tg"))

# Every file holds the tests of one area. Run a subset with, e.g.,
#     julia --project -e 'using Pkg; Pkg.test(test_args = ["motifs", "flows"])'
const GROUPS = ["core", "distances", "centralities", "temporal_betweenness", "statistics", "motifs",
                "reference_models", "connectivity", "flows", "trees", "isomorphism", "datasets", "drawing", "quality"]
const SELECTED = isempty(ARGS) ? GROUPS : filter(in(ARGS), GROUPS)
isempty(ARGS) || length(SELECTED) == length(ARGS) ||
    error("unknown test groups $(setdiff(ARGS, GROUPS)); available: $(join(GROUPS, ", "))")

@testset "TemporalGraphs.jl" begin
    for group in SELECTED
        include("$group.jl")
    end
end
