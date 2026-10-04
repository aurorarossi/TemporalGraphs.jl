# Temporal graph isomorphisms.

@testset "temporal graph isomorphisms (Heeg, Sauer, Mutzel, Scholtes)" begin
    # Figure 2 of the paper: δ = 2 with unit transition times, i.e., β = 1; a..e = 1..5
    G(es) = OrderedEdgeList(5, [(u, v, t, 1) for (u, v, t) in es])
    G1 = G([(1, 2, 1), (2, 4, 2), (3, 4, 3), (4, 5, 4)])
    G2 = G([(1, 2, 2), (2, 4, 3), (3, 4, 4), (4, 5, 5)])
    G3 = G([(1, 2, 1), (2, 4, 3), (3, 4, 2), (4, 5, 4)])
    G4 = G([(1, 2, 1), (2, 4, 2), (3, 4, 3), (4, 5, 1)])
    G5 = G([(1, 2, 1), (2, 4, 2), (3, 4, 3), (4, 5, 4), (4, 5, 5)])
    iso(H, kind) = is_temporally_isomorphic(G1, H; kind=kind, β=1)
    @test [iso(H, :concatenated) for H in (G2, G3, G4, G5)] == [true, false, false, false]
    @test [iso(H, :event) for H in (G2, G3, G4, G5)] == [true, true, false, false]
    @test [iso(H, :aggregated) for H in (G2, G3, G4, G5)] == [true, true, true, false]
    @test temporal_isomorphism(G1, G3; β=1) == (nodes=[1, 2, 3, 4, 5], edges=[1, 3, 2, 4])
    @test temporal_wl_equivalent(G1, G3; β=1) && !temporal_wl_equivalent(G1, G4; β=1)
    A = augmented_event_graph(G1; β=1)
    @test Graphs.nv(A.graph) == 9 && Graphs.ne(A.graph) == 3 + 2 * 4 && A.labels == [0, 0, 0, 0, 0, 1, 1, 1, 1]
    # against the definitions
    rng = MersenneTwister(23)
    for _ in 1:150
        n = rand(rng, 1:5)
        es = [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:6), 1) for _ in 1:rand(rng, 0:6)]
        g1 = OrderedEdgeList(n, es)
        p = randperm(rng, n)
        mode = rand(rng, 1:4)
        shift = rand(rng, 0:3)
        es2 = [TemporalEdge(p[e.u], p[e.v], mode == 1 ? e.t + shift : mode == 2 ? e.t + rand(rng, 0:1) : mode == 3 ? rand(rng, 0:6) : e.t, 1)
               for e in es]
        g2 = OrderedEdgeList(n, es2)
        for (kind, β) in ((:event, 0), (:event, 1), (:event, 3), (:aggregated, 0), (:concatenated, 0))
            m = temporal_isomorphism(g1, g2; kind=kind, β=β)
            o = oracle_temporal_isomorphic(g1, g2, kind, β)
            @test (m !== nothing) == o
            o && @test temporal_wl_equivalent(g1, g2; kind=kind, β=β)
            if m !== nothing && kind == :event
                e1, e2 = edges(g1), edges(g2)
                @test isperm(m.nodes) && isperm(m.edges)
                @test all((e2[m.edges[i]].u, e2[m.edges[i]].v) == (m.nodes[e1[i].u], m.nodes[e1[i].v]) for i in eachindex(e1))
                @test all(_oracle_follows(e1[i], e1[j], β) == _oracle_follows(e2[m.edges[i]], e2[m.edges[j]], β)
                          for i in eachindex(e1), j in eachindex(e1))
            end
        end
        oracle_temporal_isomorphic(g1, g2, :concatenated, 0) && @test oracle_temporal_isomorphic(g1, g2, :event, 0)
        oracle_temporal_isomorphic(g1, g2, :event, 0) && @test oracle_temporal_isomorphic(g1, g2, :aggregated, 0)
    end
    # the WL kernel: isomorphic graphs have the same features; the event graph alone
    # (Oettershagen et al.) cannot tell which nodes the temporal edges join
    star = OrderedEdgeList(4, [(1, 2, 1, 1), (1, 3, 1, 1)])
    disjoint = OrderedEdgeList(4, [(1, 2, 1, 1), (3, 4, 1, 1)])
    K = temporal_wl_kernel([G1, G2, G3, G4]; β=1)
    @test K.kernel == K.kernel' && K.features[1] == K.features[2] == K.features[3] && K.features[1] != K.features[4]
    @test temporal_wl_kernel([star, disjoint]; augmented=false).features[1] == temporal_wl_kernel([star, disjoint]; augmented=false).features[2]
    @test temporal_wl_kernel([star, disjoint]).features[1] != temporal_wl_kernel([star, disjoint]).features[2]
    @test_throws ArgumentError temporal_isomorphism(G1, G2; kind=:other)
end
