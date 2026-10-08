# Generators of temporal graph classes and temporal graph parameters.

# The static edges (ordered or unordered pairs) of g with their time stamps.
labels_of(g; directed=true) = Set((directed ? (e.u, e.v) : minmax(e.u, e.v), e.t) for e in edges(g))
# Whether every edge of g has its reverse with the same times.
symmetric(g) = Set((e.u, e.v, e.t, e.tt) for e in edges(g)) == Set((e.v, e.u, e.t, e.tt) for e in edges(g))
# Whether g, undirected, is temporally connected and every undirected edge is needed.
function minimal_connected(g)
    is_temporally_connected(g) || return false
    pairs = unique([(minmax(e.u, e.v), e.t) for e in edges(g)])
    return all(pairs) do (p, t)
        rest = [e for e in edges(g) if (minmax(e.u, e.v), e.t) != (p, t)]
        !is_temporally_connected(OrderedEdgeList(num_nodes(g), rest))
    end
end

@testset "generators of temporal graph classes" begin
    rng = MersenneTwister(5)
    for directed in (false, true), n in (0, 1, 2, 7)
        g = random_temporal_graph(n, 0.4, 6; directed=directed, rng=rng)
        @test num_nodes(g) == n && time_interval(g) == (1, 7)
        @test all(e -> e.u != e.v && 1 <= e.t <= 6 && e.tt == 1, edges(g))
        @test allunique(edges(g)) && (directed || symmetric(g))
        @test num_edges(random_temporal_graph(n, 0, 6; directed=directed)) == 0
        @test num_edges(random_temporal_graph(n, 1, 6; directed=directed)) == 6 * n * (n - 1)
    end
    # the expected number of edges
    m = sum(num_edges(random_temporal_graph(30, 0.1, 20; rng=rng)) for _ in 1:50) / 50
    @test abs(m - 0.1 * 30 * 29 * 20) < 0.05 * 0.1 * 30 * 29 * 20
    @test random_temporal_graph(10, 0.3, 4; rng=MersenneTwister(1)) == random_temporal_graph(10, 0.3, 4; rng=MersenneTwister(1))
    @test all(e -> e.tt == 0, edges(random_temporal_graph(5, 0.5, 3; transition_time=0, rng=rng)))
    for directed in (false, true), p in (0.3, 1.0)
        g = random_simple_temporal_graph(9, p; directed=directed, rng=rng)
        ts = sort(unique(e.t for e in edges(g)))
        @test ts == 1:length(ts) && is_simple(g; directed=directed) && is_proper(g; directed=directed)
        @test num_edges(g) == length(ts) * (directed ? 1 : 2) && (directed || symmetric(g))
        p == 1 && @test length(ts) == (directed ? 72 : 36) && is_temporally_connected(g)  # a temporal clique
    end
    for G in (Graphs.path_graph(6), Graphs.star_graph(5), Graphs.cycle_digraph(4)), k in (1, 3)
        g = random_temporal_labeling(G, 5; labels=k, rng=rng)
        @test num_nodes(g) == Graphs.nv(G) && time_interval(g) == (1, 6)
        L = labels_of(g; directed=Graphs.is_directed(G))
        @test length(L) == k * Graphs.ne(G) && all(((p, t),) -> 1 <= t <= 5, L)
        @test Graphs.is_directed(G) ? num_edges(g) == k * Graphs.ne(G) : symmetric(g)
        k == 1 && @test is_simple(g; directed=Graphs.is_directed(G))
    end
    @test_throws ArgumentError random_temporal_labeling(Graphs.path_graph(3), 2; labels=3)
    @test_throws ArgumentError random_temporal_graph(3, 1.5, 2)
    @test_throws ArgumentError random_temporal_graph(3, 0.5, 2; transition_time=-1)
    for n in 0:9
        g = round_robin_temporal_clique(n)
        L = collect(labels_of(g; directed=false))
        @test sort(first.(L)) == [(u, v) for u in 1:n for v in u+1:n]  # every pair exactly once
        @test isempty(L) || maximum(last.(L)) == (isodd(n) ? n : n - 1)
        @test symmetric(g) && is_proper(g; directed=false) && is_temporally_connected(g)
    end
    for d in 0:5
        for g in (temporal_hypercube(d), temporal_knodel_graph(d))
            @test num_nodes(g) == 2^d && length(labels_of(g; directed=false)) == (d == 0 ? 0 : d * 2^(d - 1))
            @test symmetric(g) && is_simple(g; directed=false) && is_proper(g; directed=false)
            @test !is_proper(g) || d == 0
            @test d > 4 ? is_temporally_connected(g) : minimal_connected(g)
        end
    end
end

@testset "temporal graph parameters" begin
    rng = MersenneTwister(8)
    for _ in 1:100
        n = rand(rng, 1:8)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:10), rand(rng, 0:2)) for _ in 1:rand(rng, 0:15)])
        es = edges(g)
        times = 0:10
        active(ts, t) = !isempty(ts) && minimum(ts) <= t <= maximum(ts)
        vim = maximum([count(v -> active([e.t for e in es if v in (e.u, e.v)], t), 1:n) for t in times])
        @test vertex_interval_membership_width(g) == vim
        for directed in (true, false)
            key(e) = directed ? (e.u, e.v) : minmax(e.u, e.v)
            pairs = unique(key.(es))
            eim = maximum([count(p -> active([e.t for e in es if key(e) == p], t), pairs) for t in times])
            @test edge_interval_membership_width(g; directed=directed) == eim
            @test is_simple(g; directed=directed) == all(p -> length(unique(e.t for e in es if key(e) == p)) <= 1, pairs)
            proper = all(es) do e
                all(f -> key(f) == key(e) || f.t != e.t || isempty(intersect((e.u, e.v), (f.u, f.v))), es)
            end
            @test is_proper(g; directed=directed) == proper
        end
        s = get_statistics(g)
        @test s.vertex_interval_membership_width == vim
        @test s.edge_interval_membership_width == edge_interval_membership_width(g)
        @test s.temporality == maximum([length(unique(e.t for e in es if (e.u, e.v) == p)) for p in unique((e.u, e.v) for e in es)]; init=0)
        @test s.max_edges_per_time_stamp == maximum([count(e -> e.t == t, es) for t in times])
    end
    # the edges (1, 2, 1) and (2, 1, 1) are adjacent when directed, the same edge otherwise
    g = make_undirected(OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1)]))
    @test !is_proper(g) && is_proper(g; directed=false) && is_simple(g)
    @test !is_proper(OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 1, 1)]))
    @test !is_simple(OrderedEdgeList(2, [(1, 2, 1, 1), (1, 2, 3, 1)]))
end
