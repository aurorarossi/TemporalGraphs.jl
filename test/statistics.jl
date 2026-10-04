# Burstiness, clustering coefficient, topological overlap and cores.

@testset "burstiness" begin
    g = paper_file()
    b = edge_burstiness(g)
    @test b == Dict((1, 2) => -1.0, (3, 4) => -1.0)
    @test node_burstiness(g) ≈ [-1 / 3, 0, -1, 0]
    @test all(<(-0.8), values(edge_burstiness(non_bursty(1000))))
    @test all(<(-0.8), node_burstiness(non_bursty(1000)))
    rng = MersenneTwister(17)
    bg = bursty(rng, 1000)
    @test all(>(0.8), values(edge_burstiness(bg)))
    @test all(>(0.8), node_burstiness(bg))
    rg = random_graph(rng, 2, 1000)
    @test all(x -> -0.5 < x < 0.5, values(edge_burstiness(rg)))
end

@testset "clustering coefficient and topological overlap" begin
    il = to_incident_lists(paper_file())
    @test temporal_clustering_coefficient(il, 1) ≈ 0.045454545454545456
    @test temporal_clustering_coefficient(il) ≈ [0.045454545454545456, 0, 0, 0]
    @test temporal_clustering_coefficient(to_incident_lists(chain()), 3) == 0
    k5 = complete_graph()
    @test all(temporal_clustering_coefficient(k5, v, (1, 2)) == 1 for v in 1:5)
    @test temporal_clustering_coefficient(k5, (1, 2)) == ones(5)

    @test topological_overlap(il, 1) == 0.5
    @test topological_overlap(no_topological_overlap(200), 1) == 0
    @test topological_overlap(full_topological_overlap(10), 1) == 1
    @test topological_overlap(il) == [0.5, 0.0, 1.0, 0.0]
end

@testset "cores" begin
    g = lkcore_tg()
    # TGLib's static k-core only decreases the degrees of out-neighbors and returns
    # [3, 3, 3] for the first nodes; these are the correct undirected core numbers.
    @test kcores(to_aggregated_edge_list(g), 5) == [2, 2, 2, 1, 2]
    @test temporal_khcores(g, 1) == [2, 2, 2, 1, 2]
    @test temporal_khcores(g, 2) == [2, 2, 2, 0, 0]
    @test temporal_lkcores(g, 2, 2) == [1, 2, 3]
    @test temporal_lkcores(g, 1, 2) == [1, 2, 3]
    @test isempty(temporal_lkcores(g, 5, 1))
    rng = MersenneTwister(19)
    for _ in 1:50
        n = rand(rng, 1:12)
        pairs = [(rand(rng, 1:n), rand(rng, 1:n)) for _ in 1:rand(rng, 0:40)]
        @test kcores(pairs, n) == oracle_kcores(pairs, n)
    end
end
