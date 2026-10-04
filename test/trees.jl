# Temporal spanning trees and spanners.

@testset "spanning trees and spanners" begin
    # Huang, Fu and Liu, Figure 1 (root 0 is node 1, weights = transition times)
    g = OrderedEdgeList(6, [(1, 2, 1, 2), (1, 3, 1, 4), (1, 3, 3, 3), (1, 2, 4, 1), (2, 4, 4, 2), (2, 5, 5, 3), (3, 6, 6, 2), (3, 5, 7, 2)])
    @test sort([(e.u, e.v, e.t) for e in earliest_arrival_tree(g, 1)]) == [(1, 2, 1), (1, 3, 1), (2, 4, 4), (2, 5, 5), (3, 6, 6)]
    @test minimum_weight_spanning_tree(g, 1; level=1).weight == 12
    t = minimum_weight_spanning_tree(g, 1)
    @test t.weight == 11 && sort([(e.u, e.v, e.t) for e in t.edges]) == [(1, 2, 1), (1, 3, 3), (2, 4, 4), (3, 5, 7), (3, 6, 6)]
    # Figure 3: transition times 0, where a single pass misses node 3
    g = OrderedEdgeList(5, [(1, 2, 1, 0), (3, 1, 2, 0), (4, 2, 2, 0), (2, 5, 3, 0), (4, 3, 4, 0), (5, 4, 4, 0)])
    @test sort([(e.u, e.v, e.t) for e in earliest_arrival_tree(g, 1)]) == [(1, 2, 1), (2, 5, 3), (4, 3, 4), (5, 4, 4)]
    rng = MersenneTwister(13)
    for _ in 1:80
        n = rand(rng, 2:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:8), rand(rng, 0:2)) for _ in 1:rand(rng, 1:11)])
        es = edges(g)
        arr = oracle_tree_arrivals(g, 1, earliest_arrival_tree(g, 1))
        ea = earliest_arrival_times(g, 1)
        @test arr !== false && all(arr[v] == ea[v] for v in keys(arr))
        w = rand(rng, 0:5, length(es))
        reach = [v for v in 2:n if ea[v] < typemax(Int)]
        cands = [[i for (i, e) in enumerate(es) if e.v == v && e.u != v] for v in reach]
        opt = minimum((sum(w[collect(c)]; init=0) for c in Iterators.product(cands...)
                       if oracle_tree_arrivals(g, 1, es[collect(c)]) !== false); init=0)
        for level in (1, 2)
            res = minimum_weight_spanning_tree(g, 1; weights=w, level=level)
            @test oracle_tree_arrivals(g, 1, res.edges) !== false && res.weight >= opt && res.weight == sum(w[findall(in(res.edges), es)])
        end
        R = temporal_reachability(g)
        @test temporal_reachability(OrderedEdgeList(n, temporal_spanner(g), time_interval(g))) == R
        M = temporal_spanner(g; minimal=true)
        @test temporal_reachability(OrderedEdgeList(n, M, time_interval(g))) == R
        @test all(temporal_reachability(OrderedEdgeList(n, M[setdiff(eachindex(M), k)], time_interval(g))) != R for k in eachindex(M))
    end
    for n in (4, 9, 40, 100)
        es = TemporalEdge{Int32,Int64}[]
        for u in 1:n, v in u+1:n
            t = rand(rng, 1:3n)
            push!(es, TemporalEdge(u, v, t, 0), TemporalEdge(v, u, t, 0))
        end
        g = OrderedEdgeList(n, es)
        S = temporal_clique_spanner(g)
        @test all(temporal_reachability(OrderedEdgeList(n, S, time_interval(g))))
        n >= 40 && @test length(S) ÷ 2 <= n * log2(n)
    end
    @test_throws ArgumentError temporal_clique_spanner(OrderedEdgeList(3, [(1, 2, 1, 0), (2, 3, 1, 0)]))
    @test_throws ArgumentError minimum_weight_spanning_tree(g, 1; level=3)
end
