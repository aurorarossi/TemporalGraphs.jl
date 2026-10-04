# Dataset catalog and the GNNGraphs.jl extension (MLDatasets.jl is tested in test/extensions).

@testset "datasets and GNNGraphs.jl" begin
    # the catalog, and loading from the cache without downloading
    cat = temporal_datasets()
    @test length(cat) >= 30 && allunique(d.name for d in cat) && all(d -> startswith(d.url, "http"), cat)
    @test occursin("CollegeMsg", sprint(show, MIME("text/plain"), first(cat)))
    @test_throws ArgumentError load_dataset("no-such-dataset")
    dir = mktempdir()
    mkpath(joinpath(dir, "CollegeMsg"))
    write(joinpath(dir, "CollegeMsg", "CollegeMsg.txt"), "10 20 100\n20 30 105\n% comment\n10 30 110\n")
    g = load_dataset("CollegeMsg"; dir=dir)
    @test num_nodes(g) == 3 && edges(g) == [TemporalEdge(1, 2, 100, 1), TemporalEdge(2, 3, 105, 1), TemporalEdge(1, 3, 110, 1)]
    @test original_ids(g) == [10, 20, 30]
    mkpath(joinpath(dir, "sp-hospital"))
    write(joinpath(dir, "sp-hospital", "hospital_lyon_contacts.dat"), "140\t1157\t1232\tMED\tADM\n160\t1157\t1191\tMED\tMED\n")
    h = load_dataset("sp-hospital"; dir=dir, transition_time=20)
    @test num_nodes(h) == 3 && num_edges(h) == 4 && all(e -> e.tt == 20, edges(h))  # undirected contacts
    @test num_edges(load_dataset("sp-hospital"; dir=dir, directed=true)) == 2
    mkpath(joinpath(dir, "soc-sign-bitcoin-otc"))
    write(joinpath(dir, "soc-sign-bitcoin-otc", "soc-sign-bitcoinotc.csv"), "6,2,4,1289241911.72836\n6,5,2,1289241941.53378\n")
    @test [e.t for e in edges(load_dataset("soc-sign-bitcoin-otc"; dir=dir))] == [1289241912, 1289241942]

    # GNNGraphs.jl
    g = OrderedEdgeList(4, [(1, 2, 0, 1), (2, 3, 0, 1), (1, 3, 5, 1), (2, 1, 7, 2), (4, 4, 7, 1)])
    G = GNNGraphs.GNNGraph(g)
    @test G.num_nodes == 4 && G.num_edges == 5 && G.edata.t == [0, 0, 5, 7, 7]
    @test edges(OrderedEdgeList(G)) == edges(g)
    @test_throws ArgumentError OrderedEdgeList(GNNGraphs.GNNGraph([1], [2]))
    TG = GNNGraphs.TemporalSnapshotsGNNGraph(g; resolution=5)
    @test TG.num_snapshots == 2 && TG.tgdata.times == [0, 5] && TG.num_nodes == [4, 4]
    @test Set((e.u, e.v, e.t) for e in edges(OrderedEdgeList(TG))) == Set((e.u, e.v, 5 * fld(e.t, 5)) for e in edges(g))
    A = augmented_event_gnngraph(g; β=2)
    B = augmented_event_graph(g; β=2)
    @test A.num_nodes == Graphs.nv(B.graph) && A.num_edges == Graphs.ne(B.graph) && A.ndata.label == B.labels
    @test vec(sum(A.ndata.x; dims=1)) == ones(Float32, A.num_nodes)

    # MLDatasets.jl is tested in test/extensions (it is a heavy dependency)
end
