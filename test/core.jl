# Temporal edges, input and output, transformations, snapshots and integration with other packages.

@testset "basic types" begin
    e1 = TemporalEdge(1, 2, 1, 1)
    e2 = TemporalEdge(1, 2, 2, 1)
    @test e1 < e2 && e1 != e2 && e1 == e1
    es = Set([TemporalEdge(1, 2, 1, 1), TemporalEdge(1, 2, 1, 1), TemporalEdge(1, 2, 2, 1),
              TemporalEdge(1, 3, 2, 1), TemporalEdge(1, 3, 2, 2), TemporalEdge(2, 3, 2, 2)])
    @test length(es) == 5
    v = sort([TemporalEdge(1, 2, 4, 1), TemporalEdge(1, 2, 98, 1), TemporalEdge(1, 2, 2, 1),
              TemporalEdge(1, 3, 1, 1), TemporalEdge(1, 3, 9, 2), TemporalEdge(2, 3, 9, 2)])
    @test issorted([e.t for e in v])
    @test TemporalEdge(1, 2, 3) == TemporalEdge(1, 2, 3, 1)
    @test sprint(show, TemporalEdge(1, 2, 3, 4)) == "(1 2 3 4)"
    @test_throws ArgumentError OrderedEdgeList(2, [TemporalEdge(1, 3, 1, 1)])
    g = OrderedEdgeList(3, [TemporalEdge(1, 2, 5, 1), TemporalEdge(2, 3, 1, 1)])
    @test issorted(edges(g)) && time_interval(g) == (1, 6)
    # input sorted by time is kept, but transition time 0 comes first at every time stamp
    g = OrderedEdgeList(7, [(6, 1, 1, 1), (1, 7, 1, 3), (7, 6, 1, 0), (2, 3, 2, 1), (3, 4, 2, 0)])
    @test edges(g) == [TemporalEdge(7, 6, 1, 0), TemporalEdge(6, 1, 1, 1), TemporalEdge(1, 7, 1, 3),
                       TemporalEdge(3, 4, 2, 0), TemporalEdge(2, 3, 2, 1)]
    @test earliest_arrival_times(g, 7)[1] == 2
    @test scale_timestamps(g, 3) == scale_timestamps(g, 3.0)
end

@testset "input and output" begin
    g = paper_file()
    @test num_nodes(g) == 4 && num_edges(g) == 7
    @test edges(g)[1] == TemporalEdge(1, 3, 1, 5)
    @test original_id(g, 1) == 0 && original_id(g, 3) == 3
    @test node_map(g)[3] == 3
    @test time_interval(g) == (1, 12)
    e = edges(g)[1]
    @test (original_id(g, e.u), original_id(g, e.v)) == (0, 3)

    s = get_statistics(g)
    @test repr(MIME"text/plain"(), s) == """
        number of nodes: 4
        number of edges: 7
        number of static edges: 5
        number of time stamps: 6
        number of transition times: 4
        min. time stamp: 1
        max. time stamp: 8
        min. transition time: 1
        max. transition time: 5
        min. temporal in-degree: 0
        max. temporal in-degree: 3
        min. temporal out-degree: 1
        max. temporal out-degree: 3"""

    mktempdir() do dir
        path = joinpath(dir, "g.tg")
        save_ordered_edge_list(path, g)
        h = load_ordered_edge_list(path)
        orig(x) = Set((original_id(x, e.u), original_id(x, e.v), e.t, e.tt) for e in edges(x))
        @test num_nodes(h) == num_nodes(g) && orig(h) == orig(g)
        write(path, "# comment\n10\t20\t3\n20,30,4,2\r\n\n5 6\n-1 10 7 1 extra 9.5\n")
        h = load_ordered_edge_list(path)
        @test num_nodes(h) == 4 && num_edges(h) == 3
        @test original_ids(h) == [10, 20, 30, -1]
        @test edges(h) == [TemporalEdge(1, 2, 3, 1), TemporalEdge(2, 3, 4, 2), TemporalEdge(4, 1, 7, 1)]
        u = load_ordered_edge_list(path; directed=false)
        @test num_edges(u) == 6
        write(path, "1 2 x\n")
        @test_throws ArgumentError load_ordered_edge_list(path)
    end
    @test num_edges(load_incident_lists(joinpath(DATA, "example_from_paper.tg"))) == 7
    @test num_trs_nodes(load_trs_graph(joinpath(DATA, "example_from_paper.tg"))) == 9
end

@testset "transformations" begin
    g = example_from_paper()
    trs = to_trs_graph(g)
    @test num_trs_nodes(trs) == 9
    @test num_edges(trs) == 5 + 7
    @test [length(trs_neighbors(trs, i)[1]) for i in 1:9] == [2, 2, 1, 1, 2, 0, 2, 2, 0]
    heads, tts = trs_neighbors(trs, 1)
    @test heads == [2, 7] && tts == [0, 5]

    @test normalize_graph(g, true) == g
    @test normalize_graph(example_from_paper_with_loops_and_multi_edges(), true) == g
    nl = normalize_graph(example_from_paper_with_loops_and_multi_edges(), false)
    @test edges(nl)[3].u == edges(nl)[3].v
    rng = MersenneTwister(1)
    r = OrderedEdgeList(100, [TemporalEdge(rand(rng, 1:100), rand(rng, 1:100), rand(rng, 0:99), rand(rng, 1:100)) for _ in 1:1000])
    @test issorted(edges(normalize_graph(r)))

    s = scale_timestamps(g, 10)
    @test num_nodes(s) == num_nodes(g) && num_edges(s) == num_edges(g)
    @test edges(s)[1].t == 10 * edges(g)[1].t
    @test edges(scale_timestamps(g, 0.25))[1].t == 0      # round(0.25) == 0
    @test edges(scale_timestamps(g, 0.5))[1].t == 1       # ties away from zero
    @test all(e -> e.tt == 1, edges(unit_transition_times(normalize_graph(r, false), 1)))

    u = make_undirected(g)
    @test num_edges(u) == 14
    @test all(e -> TemporalEdge(e.v, e.u, e.t, e.tt) in Set(edges(u)), edges(u))

    agg = to_aggregated_edge_list(g)
    @test [(e.u, e.v, e.weight) for e in agg] == [(1, 2, 2), (1, 4, 1), (2, 4, 1), (3, 2, 1), (4, 3, 2)]

    il = to_incident_lists(g)
    @test num_edges(il) == 7 && num_nodes(il) == 4
    @test to_ordered_edge_list(il) == g
    @test out_degree(il, 1) == 3 && issorted([e.t for e in out_edges(il, 1)])

    # reversed graph: durations are transposed
    rng = MersenneTwister(2)
    for (n, m) in ((30, 300), (8, 60))
        rg = random_general_graph(rng, n, m; min_tt=0)
        rev = reverse(rg)
        @test time_interval(rev) == time_interval(rg)
        ds = [minimum_durations(rg, u) for u in 1:n]
        rds = [minimum_durations(rev, u) for u in 1:n]
        @test all(ds[u][v] == rds[v][u] for u in 1:n, v in 1:n)
        @test num_edges(reverse(to_incident_lists(rg))) == num_edges(rev)
    end

    dlg = to_directed_line_graph(g)
    @test num_nodes(dlg) == 7
    for i in 1:7, j in 1:7
        e, f = edges(g)[i], edges(g)[j]
        @test (j in out_edges(dlg, i)) == (e.v == f.u && f.t >= e.t + e.tt)
    end
    dil = to_directed_line_graph(il)
    @test num_edges(dil) == num_edges(dlg)
end

@testset "snapshot graphs" begin
    rng = MersenneTwister(31)
    for _ in 1:80
        n = rand(rng, 1:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:20), rand(rng, 0:2)) for _ in 1:rand(rng, 0:15)])
        ti = rand(rng, Bool) ? time_interval(g) : (rand(rng, 0:4), rand(rng, 8:22))
        a, b = ti
        for Δ in (nothing, 1, 3, 7), directed in (true, false)
            S = snapshots(g, ti; resolution=Δ, directed=directed)
            D = Δ === nothing ? TemporalGraphs._time_grid(g, nothing)[3] : Δ
            @test S.times == collect(a:D:b)
            key(u, v) = directed ? (Int(u), Int(v)) : minmax(Int(u), Int(v))
            @test all(Set(key(e.u, e.v) for e in edges(g) if t <= e.t < t + D && a <= e.t <= b) ==
                      Set(key(Graphs.src(e), Graphs.dst(e)) for e in Graphs.edges(S.graphs[k])) for (k, t) in enumerate(S.times))
            U = directed ? Graphs.SimpleDiGraph(n) : Graphs.SimpleGraph(n)
            for h in S.graphs, e in Graphs.edges(h)
                Graphs.add_edge!(U, Graphs.src(e), Graphs.dst(e))
            end
            @test static_graph(g, ti; directed=directed) == U
        end
        S = snapshots(g)
        @test Set((e.u, e.v, e.t) for e in edges(OrderedEdgeList(S.graphs, S.times))) == Set((e.u, e.v, e.t) for e in edges(g))
        for Δ in (1, 2, 5)
            A = aggregate_time(g, Δ, ti)
            @test Set((e.u, e.v, e.t, e.tt) for e in edges(A)) == Set((e.u, e.v, a + Δ * fld(e.t - a, Δ), e.tt) for e in edges(g) if a <= e.t <= b)
            @test allunique(edges(A))
            @test num_edges(aggregate_time(g, Δ, ti; multiedges=true)) == count(e -> a <= e.t <= b, edges(g))
            @test all(e -> e.tt == Δ, edges(aggregate_time(g, Δ, ti; transition_time=Δ)))
        end
    end
    h = OrderedEdgeList([Graphs.path_graph(3), Graphs.SimpleGraph(3)], [10, 20])
    @test Set((e.u, e.v, e.t) for e in edges(h)) == Set([(1, 2, 10), (2, 1, 10), (2, 3, 10), (3, 2, 10)]) && time_interval(h) == (10, 21)
    f = OrderedEdgeList(3, [(1, 2, 0.5, 0.0), (2, 3, 0.5, 0.0), (1, 3, 1.25, 0.0)])
    @test snapshots(f).times == [0.5, 1.25] && Graphs.ne.(snapshots(f).graphs) == [2, 1]
    @test_throws ArgumentError OrderedEdgeList([Graphs.SimpleGraph(2), Graphs.SimpleGraph(3)])
    @test_throws ArgumentError OrderedEdgeList([Graphs.SimpleGraph(2)], [1, 2])
    @test_throws ArgumentError aggregate_time(f, 0)
end

@testset "generic node and time types" begin
    rng = MersenneTwister(29)
    for _ in 1:20
        g = random_general_graph(rng, rand(rng, 2:9), rand(rng, 0:50); min_tt=0)
        es = edges(g)
        gf = OrderedEdgeList(num_nodes(g), [TemporalEdge{Int32,Float64}(e.u, e.v, e.t, e.tt) for e in es])
        g64 = OrderedEdgeList(num_nodes(g), [TemporalEdge{Int64,Int64}(e.u, e.v, e.t, e.tt) for e in es])
        @test time_type(gf) == Float64 && node_type(g64) == Int64
        conv(d) = [x == typemax(x) ? Inf : x == typemin(x) ? -Inf : Float64(x) for x in d]
        for s in 1:num_nodes(g), dt in (EarliestArrival(), Fastest(), LatestDeparture(), MinimumTransitionTimes())
            d = temporal_distances(g, s, dt)
            for x in (gf, to_incident_lists(gf))
                @test temporal_distances(x, s, dt) == conv(d)
            end
            @test temporal_distances(g64, s, dt) == d
            dt isa LatestDeparture || @test temporal_distances(to_trs_graph(gf), s, dt) == conv(d)
        end
        @test minimum_hops(gf, 1) == minimum_hops(g, 1)
        @test temporal_closeness(gf, Fastest()) ≈ temporal_closeness(g, Fastest())
        @test temporal_walk_centrality(gf, 0.3, 0.7) ≈ temporal_walk_centrality(g, 0.3, 0.7)
        @test temporal_katz_centrality(gf, 0.3) ≈ temporal_katz_centrality(g, 0.3)
    end
    # real valued times
    g = OrderedEdgeList(3, [(1, 2, 0.5, 0.25), (2, 3, 0.75, 0.1)])
    @test time_type(g) == Float64
    @test earliest_arrival_times(g, 1) == [0.0, 0.75, 0.85]
    @test minimum_durations(g, 1) ≈ [0.0, 0.25, 0.35]
    @test minimum_hops(g, 3) == [typemax(Int), typemax(Int), 0]
    mktempdir() do dir
        path = joinpath(dir, "f.tg")
        write(path, "1 2 0.5 0.25\n2 3 0.75 0.1\n")
        @test load_ordered_edge_list(path; time_type=Float64) == g
    end
end

@testset "Tables.jl and Graphs.jl integration" begin
    tbl = (src=[10, 20, 10], dst=[20, 30, 30], time=[1, 2, 5], dur=[1, 1, 2])
    g = OrderedEdgeList(tbl; u=:src, v=:dst, t=:time, tt=:dur)
    @test num_nodes(g) == 3 && num_edges(g) == 3
    @test original_ids(g) == [10, 20, 30]
    @test edges(g) == [TemporalEdge(1, 2, 1, 1), TemporalEdge(2, 3, 2, 1), TemporalEdge(1, 3, 5, 2)]
    @test OrderedEdgeList(tbl; directed=false) |> num_edges == 6
    @test time_type(OrderedEdgeList((a=[1], b=[2], c=[0.5]))) == Float64
    @test_throws ArgumentError OrderedEdgeList(42)
    sg = static_graph(g)
    @test sg isa Graphs.SimpleDiGraph && Graphs.ne(sg) == 3 && Graphs.has_edge(sg, 1, 3)
    @test Graphs.ne(static_graph(make_undirected(g); directed=false)) == 3
    @test Graphs.nv(g) == 3 && Graphs.ne(g) == 3 && Graphs.edges(g) === edges(g)
end
