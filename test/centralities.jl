# Closeness, edge betweenness, walk-based centralities and top-k closeness.

@testset "closeness" begin
    g = paper_file()
    @test temporal_closeness(g, Fastest()) ≈ [1.3928571428571428, 0.5, 0.5, 1.3333333333333333]
    @test temporal_closeness(chain(), 1, Fastest()) >= 2.5928
    ex = example_from_paper()
    il = to_incident_lists(ex)
    trs = to_trs_graph(ex)
    for v in 1:4
        c = temporal_closeness(ex, v, Fastest())
        @test temporal_closeness(il, v, Fastest()) == c
        @test temporal_closeness(trs, v, Fastest()) == c
    end
    rng = MersenneTwister(5)
    g = random_general_graph(rng, 30, 400)
    for dt in (EarliestArrival(), Fastest(), LatestDeparture(), MinimumTransitionTimes(), MinimumHops())
        @test temporal_closeness(g, dt) == [temporal_closeness(g, v, dt) for v in 1:30]
    end
    @test temporal_closeness(to_incident_lists(g), MinimumHops()) ≈ temporal_closeness(g, MinimumHops())

    # top-k
    tk = compute_topk_closeness(to_incident_lists(g), 2, Fastest())
    @test [p[1] for p in compute_topk_closeness(to_incident_lists(paper_file()), 2, Fastest())] == [1, 4]
    c = compute_topk_closeness(to_incident_lists(chain()), 2, Fastest())
    @test c[1][1] == 1 && c[1][2] >= 2.5928
    for seed in 1:5, dt in (Fastest(), MinimumTransitionTimes(), EarliestArrival())
        rng = MersenneTwister(seed)
        g = seed <= 3 ? random_graph(rng, 100, 1000) : random_general_graph(rng, 60, 900; tmax=100, min_tt=seed - 4)
        il = to_incident_lists(g)
        exact = temporal_closeness(g, dt)
        topk = compute_topk_closeness(il, 10, dt)
        top = sort(exact; rev=true)
        @test length(topk) >= 10
        @test all(topk[i][2] ≈ top[i] for i in eachindex(topk))
        @test all(exact[v] ≈ c for (v, c) in topk)
        @test compute_topk_closeness(g, 10, dt) == topk
    end

    # approximation
    rng = MersenneTwister(11)
    g = random_graph(rng, 10, 100)
    approx = temporal_closeness_approximation(g, 20000; rng=MersenneTwister(1))
    exact = temporal_closeness(g, Fastest())
    @test maximum(abs.(approx .- exact)) < 0.1 * maximum(exact)
    @test temporal_closeness_approximation(to_incident_lists(g), 10) isa Vector{Float64}
end

@testset "betweenness" begin
    @test temporal_edge_betweenness(paper_file()) == zeros(7)
    b = temporal_edge_betweenness(betweenness_test())
    @test all(>=(0), b)
    @test b ≈ oracle_betweenness(to_directed_line_graph(betweenness_test()))
    rng = MersenneTwister(9)
    for _ in 1:10
        g = random_general_graph(rng, 6, 30; min_tt=0)
        dlg = to_directed_line_graph(g)
        @test temporal_edge_betweenness(dlg) ≈ oracle_betweenness(dlg)
        # same values for the incident lists numbering of the edges
        il = to_incident_lists(g)
        bil = temporal_edge_betweenness(to_directed_line_graph(il))
        @test sort(bil) ≈ sort(temporal_edge_betweenness(dlg))
    end
end

@testset "walk based centralities" begin
    g = paper_file()
    @test temporal_katz_centrality(g, 1) == [0, 3, 5, 8]
    @test temporal_pagerank(g, 1, 1, 1) == [0, 0, 0, 0]
    @test temporal_walk_centrality(g, 1, 1) == [0, 6, 6, 0]
    rng = MersenneTwister(13)
    for _ in 1:30
        g = random_general_graph(rng, 5, rand(rng, 0:14); tmax=8)
        α, β = rand(rng), rand(rng)
        @test temporal_katz_centrality(g, β) ≈ oracle_katz(g, β)
        @test temporal_walk_centrality(g, α, β) ≈ oracle_walk_centrality(g, α, β) atol = 1e-12
    end
    pr = temporal_pagerank(random_general_graph(rng, 20, 200), 0.85, 0.5, 0.9)
    @test all(isfinite, pr) && all(>=(0), pr)
end

@testset "top-k result" begin
    k, total = 50, 1000
    rng = MersenneTwister(23)
    r = TemporalGraphs.TopkResult(k)
    values = Tuple{Int,Float64}[]
    perm = randperm(rng, 10000)
    for i in 1:total
        x = Float64(perm[i])
        push!(values, (i, x))
        insert!(r, i, x)
        @test (i < k && !TemporalGraphs._is_full(r)) || (i >= k && TemporalGraphs._is_full(r))
    end
    sort!(values; by=p -> (-p[2], p[1]))
    res = sort(TemporalGraphs.results(r); by=p -> (-p[2], p[1]))
    @test res[1:k] == values[1:k]
end

@testset "closeness integrated over time (Crescenzi, Magnien, Marino)" begin
    rng = MersenneTwister(47)
    for _ in 1:80
        n = rand(rng, 2:7)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:10), rand(rng, 1:3))
                                for _ in 1:rand(rng, 1:14)])
        a, b = time_interval(g)
        for ti in ((a, b), (a + 1, b + 2))
            ti[2] > ti[1] || continue
            o = oracle_harmonic_closeness(g, ti)
            @test temporal_harmonic_closeness(g, ti) ≈ o atol = 1e-12
            @test temporal_harmonic_closeness(g, 1, ti) ≈ o[1] atol = 1e-12
        end
    end
    # a single edge: leaving 1 at τ ∈ [0, 1], the latency to 2 is 1 - τ + 1
    g = OrderedEdgeList(2, [(1, 2, 1, 1)], (0, 2))
    @test temporal_harmonic_closeness(g) ≈ [log(2) / 2, 0.0]
    gf = OrderedEdgeList(2, [(1, 2, 1.0, 1.0)], (0.0, 2.0))
    @test temporal_harmonic_closeness(gf) ≈ [log(2) / 2, 0.0]
    # the sampling estimate is unbiased
    g = random_general_graph(MersenneTwister(53), 20, 600; tmax=200, min_tt=1)
    exact = temporal_harmonic_closeness(g)
    approx = temporal_harmonic_closeness_approximation(g, 20000; rng=MersenneTwister(1))
    @test maximum(abs.(approx .- exact)) < 0.02 * maximum(exact)
    @test_throws ArgumentError temporal_harmonic_closeness(OrderedEdgeList(2, [(1, 2, 1, 0)]))
end
