# Temporal flows, cuts, disjoint paths and separators.

@testset "flows, cuts, disjoint paths and separators" begin
    rng = MersenneTwister(21)
    for _ in 1:80
        n = rand(rng, 2:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:6), rand(rng, 0:2)) for _ in 1:rand(rng, 1:10)])
        es = edges(g)
        m = length(es)
        s, z = 1, n
        P = oracle_sz_walks(es, s, z, 10^9, n)
        cap = rand(rng, 1:3, m)
        f = temporal_max_flow(g, s, z; capacity=cap)
        @test f.value == oracle_min_cut(P, m, cap) && all(0 .<= f.flow .<= cap)
        c = temporal_min_cut(g, s, z; capacity=cap)
        @test c.value == f.value == sum(cap[c.edges]; init=0) && all(w -> any(in(c.edges), w), P)
        D = temporal_edge_disjoint_paths(g, s, z)
        @test length(D) == temporal_max_flow(g, s, z).value == oracle_packing(P, (a, b) -> !isdisjoint(a, b))
        @test all(p -> p[1].u == s && p[end].v == z && allunique([p[1].u; [e.v for e in p]]) &&
                       all(p[k].v == p[k+1].u && p[k+1].t >= p[k].t + p[k].tt for k in 1:length(p)-1), D)
        @test allunique(reduce(vcat, D; init=TemporalEdge{Int32,Int64}[])) || length(unique(es)) < m
        for β in (0, 1)
            W = oracle_sz_walks(es, s, z, β, n; paths=false)
            @test temporal_max_flow(g, s, z; β=β).value == temporal_min_cut(g, s, z; β=β).value == oracle_min_cut(W, m, ones(Int, m))
            @test length(temporal_edge_disjoint_paths(g, s, z; β=β)) == temporal_max_flow(g, s, z; β=β).value
        end
        @test temporal_max_flow(g, s, z; buffers=zeros(Int, n)).value == temporal_max_flow(g, s, z; β=0).value
        O = temporal_out_disjoint_paths(g, s, z)
        deps(w) = [(es[i].u, es[i].t) for i in w]
        @test length(O) == oracle_packing(P, (a, b) -> !isdisjoint(deps(a), deps(b)))
        @test all(isdisjoint([(e.u, e.t) for e in O[i]], [(e.u, e.t) for e in O[j]]) for i in eachindex(O) for j in i+1:length(O))
        S = temporal_vertex_separator(g, s, z)
        if any(e -> e.u == s && e.v == z, es)
            @test S === nothing
        else
            inner = [v for v in 1:n if v != s && v != z]
            hits(X) = all(w -> any(i -> es[i].v in X, w[1:end-1]), P)
            @test hits(S) && length(S) == minimum(length(X) for X in
                [[inner[i] for i in eachindex(inner) if (mask >> (i - 1)) & 1 == 1] for mask in 0:(2^length(inner)-1)] if hits(X))
        end
    end
    und(es) = reduce(vcat, [[(u, v, t, tt), (v, u, t, tt)] for (u, v, t, tt) in es])
    # Akrida et al., Figures 3, 4 and 6
    g = OrderedEdgeList(4, [(1, 2, 1, 1), (3, 4, 1, 1), (1, 3, 2, 1), (3, 2, 3, 1), (2, 4, 5, 1)])
    cap = [5, 5, 5, 5, 8]
    @test temporal_max_flow(g, 1, 4; capacity=cap).value == 8
    @test edges(g)[temporal_min_cut(g, 1, 4; capacity=cap).edges] == [TemporalEdge(2, 4, 5, 1)]
    g = OrderedEdgeList(3, [(1, 2, 1, 1), (1, 2, 7, 1), (2, 3, 8, 1), (1, 2, 9, 1)])
    @test temporal_max_flow(g, 1, 3; capacity=[10, 10, 2, 10]).value == 2
    for tt in (0, 1)
        raw = [(1, 2, 1, 2), (1, 2, 2, 2), (1, 3, 2, 2), (3, 2, 4, 4), (2, 4, 2, 2), (2, 4, 8, 2), (4, 3, 3, 2),
               (2, 5, 3, 3), (2, 5, 5, 3), (4, 5, 10, 2)]
        g = OrderedEdgeList(5, [(u, v, t, tt) for (u, v, t, _) in raw])
        c = Dict((u, v, t) => x for (u, v, t, x) in raw)
        @test temporal_max_flow(g, 1, 5; capacity=[c[(e.u, e.v, e.t)] for e in edges(g)]).value == 6
    end
    # Kempe, Kleinberg and Kumar: Menger's theorem fails for vertex separators
    for tt in (0, 1)
        g = OrderedEdgeList(5, und([(1, 2, 1, tt), (2, 3, 2, tt), (3, 5, 3, tt), (2, 4, 4, tt), (1, 3, 5, tt), (3, 4, 6, tt), (4, 5, 7, tt)]))
        @test length(temporal_vertex_separator(g, 1, 5)) == 2
        @test length(temporal_edge_disjoint_paths(g, 1, 5)) == 2 == length(temporal_out_disjoint_paths(g, 1, 5))
    end
    # Mertzios, Michail and Spirakis, Figure 4
    g = OrderedEdgeList(7, [(1, 2, 1, 1), (2, 3, 2, 1), (1, 3, 3, 1), (3, 4, 5, 1), (3, 6, 5, 1), (4, 5, 6, 1), (4, 5, 8, 1),
                            (5, 6, 9, 1), (5, 7, 7, 1), (6, 7, 10, 1)])
    @test length(temporal_out_disjoint_paths(g, 1, 7)) == 1
    @test length(temporal_edge_disjoint_paths(g, 1, 7)) == 2
    @test temporal_vertex_separator(g, 1, 7) == [3]
    # Zschoche et al., Figure 1: strict and non-strict separators differ
    raw = [(1, 2, 2), (1, 2, 4), (2, 3, 3), (1, 4, 1), (1, 4, 2), (4, 2, 1), (4, 2, 2), (2, 5, 1), (2, 6, 1), (2, 6, 3),
           (6, 3, 4), (4, 5, 2), (5, 6, 2), (4, 6, 1)]
    @test temporal_vertex_separator(OrderedEdgeList(6, und([(u, v, t, 1) for (u, v, t) in raw])), 1, 3) == [2]
    @test length(temporal_vertex_separator(OrderedEdgeList(6, und([(u, v, t, 0) for (u, v, t) in raw])), 1, 3)) == 2
    @test temporal_vertex_separator(OrderedEdgeList(2, [(1, 2, 1, 1)]), 1, 2) === nothing
    @test_throws ArgumentError temporal_max_flow(g, 1, 1)
    @test_throws ArgumentError temporal_max_flow(g, 1, 7; buffers=zeros(7), β=1)
    @test_throws ArgumentError temporal_max_flow(g, 1, 7; capacity=[1])
end
