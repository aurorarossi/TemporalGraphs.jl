# Temporal motifs.

@testset "temporal motifs (Paranjape, Benson, Leskovec)" begin
    # the example graph of SNAP and the counts computed by its implementation
    raw = [(8, 0, 5), (8, 0, 10), (1, 0, 15), (0, 1, 20), (8, 1, 25), (8, 1, 30), (7, 0, 35), (0, 1, 40), (7, 1, 45),
           (5, 8, 50), (5, 8, 52), (5, 8, 54), (5, 9, 56), (5, 8, 58), (5, 9, 60), (8, 9, 65), (9, 8, 70), (9, 5, 75),
           (8, 5, 80), (8, 5, 85), (9, 5, 90), (5, 8, 99), (5, 9, 102), (8, 9, 150), (9, 8, 201)]
    g = OrderedEdgeList(10, [(u + 1, v + 1, t, 0) for (u, v, t) in raw])
    @test temporal_motif_counts(g, 10) == [0 0 0 0 2 2; 0 0 1 0 0 0; 0 0 0 0 0 0; 4 0 3 0 2 1; 0 0 0 1 1 0; 4 0 9 0 0 3]
    @test temporal_motif_counts(g, 30) == [4 13 4 9 9 2; 0 10 9 8 3 6; 1 9 0 12 0 1; 4 9 7 9 14 10; 0 3 2 8 4 1; 4 9 11 11 10 13]
    # one instance of every motif: first edge a → b, row from the 2nd edge, column from the 3rd
    E = ((1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2))
    for (p, e2) in enumerate(E), (q, e3) in enumerate(E)
        M = temporal_motif_counts(OrderedEdgeList(3, [(1, 2, 1, 0), (e2..., 2, 0), (e3..., 3, 0)]), 2)
        @test M[7-p, q] == 1 && sum(M) == 1
        @test temporal_motif_count(OrderedEdgeList(3, [(1, 2, 1, 0), (e2..., 2, 0), (e3..., 3, 0)]),
                                   [(1, 2), e2, e3], 2) == 1
    end
    # brute force, with ties, loops and multi-edges
    rng = MersenneTwister(5)
    motifs = [[(1, 2), (2, 3), (3, 1)], [(1, 2), (2, 1), (1, 2), (2, 1)], [(1, 2), (1, 3), (1, 2)], [(1, 2), (2, 3), (3, 4)],
              [(1, 2)], [(1, 2), (1, 2)], [(1, 2), (3, 2), (2, 3), (1, 3)]]
    for _ in 1:100
        n = rand(rng, 2:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:10), 1) for _ in 1:rand(rng, 1:12)])
        for δ in (0, 2, 5, 100), strict in (true, false)
            @test temporal_motif_counts(g, δ; strict=strict) == oracle_motif_counts(g, δ; strict=strict)
        end
        for mot in motifs[rand(rng, 1:length(motifs), 2)], δ in (1, 100), strict in (true, false)
            @test temporal_motif_count(g, mot, δ; strict=strict) == oracle_motif_count(g, mot, δ; strict=strict)
        end
    end
    # the general algorithm agrees with the fast ones on a larger graph
    g = random_graph(MersenneTwister(3), 25, 2000; ti=(0, 2000))
    M = temporal_motif_counts(g, 150)
    for (p, e2) in enumerate(E), (q, e3) in enumerate(E)
        (p, q) in ((1, 1), (3, 4), (5, 2), (6, 6)) || continue
        @test temporal_motif_count(g, [(1, 2), e2, e3], 150) == M[7-p, q]
    end
    @test_throws ArgumentError temporal_motif_counts(g, -1)
    @test_throws ArgumentError temporal_motif_count(g, [(1, 2), (3, 4)], 10)
    @test_throws ArgumentError temporal_motif_count(g, [(1, 1)], 10)
    @test_throws ArgumentError temporal_motif_count(g, Tuple{Int,Int}[], 10)
end
