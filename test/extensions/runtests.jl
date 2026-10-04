# Tests of the MLDatasets.jl extension, in their own environment because MLDatasets is a
# heavy dependency. From the package directory (Julia ≥ 1.11, which reads [sources]):
#     julia --project=test/extensions -e 'using Pkg; Pkg.instantiate()'
#     julia --project=test/extensions test/extensions/runtests.jl
# On Julia 1.10 first run `Pkg.develop(path=".")` in that environment.
using TemporalGraphs
using Test
import MLDatasets

@testset "MLDatasets.jl extension" begin
    g = OrderedEdgeList(4, [(1, 2, 0, 1), (2, 3, 0, 1), (1, 3, 5, 1), (2, 1, 7, 2), (4, 4, 7, 1)])
    M = MLDatasets.TemporalSnapshotsGraph(g; resolution=5)
    @test M.num_snapshots == 2 && M.num_edges == [2, 3] && M.graph_data.times == [0, 5]
    @test Set((e.u, e.v, e.t) for e in edges(OrderedEdgeList(M, M.graph_data.times))) == Set((e.u, e.v, 5 * fld(e.t, 5)) for e in edges(g))
    x = MLDatasets.TemporalSnapshotsGraph(num_nodes=[3, 3, 3], edge_index=([1, 2, 3], [2, 3, 1], [1, 2, 3]))
    @test edges(OrderedEdgeList(x)) == [TemporalEdge(1, 2, 1, 1), TemporalEdge(2, 3, 2, 1), TemporalEdge(3, 1, 3, 1)]
    s = MLDatasets.Graph(num_nodes=3, edge_index=([1, 2], [2, 3]), edge_data=(timestamp=[10, 20],))
    @test edges(OrderedEdgeList(s)) == [TemporalEdge(1, 2, 10, 1), TemporalEdge(2, 3, 20, 1)]
    hg = MLDatasets.HeteroGraph(num_nodes=Dict("user" => 2, "movie" => 3),
                                edge_indices=Dict(("user", "rating", "movie") => ([1, 2], [3, 1])),
                                edge_data=Dict(("user", "rating", "movie") => Dict(:timestamp => [100, 50])), node_data=Dict())
    r = OrderedEdgeList(hg, ("user", "rating", "movie"))
    @test num_nodes(r) == 5 && edges(r) == [TemporalEdge(2, 3, 50, 1), TemporalEdge(1, 5, 100, 1)]
end
