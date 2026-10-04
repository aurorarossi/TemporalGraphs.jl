# Temporal betweenness for optimal walks, its proxies and its approximations.

@testset "temporal betweenness (Brunelli, Crescenzi, Viennot)" begin
    # Figure 1 and Table 2 of the paper, β = 1
    g = OrderedEdgeList(8, [(1, 2, 1, 1), (2, 3, 2, 2), (1, 2, 4, 1), (3, 4, 4, 1), (2, 5, 4, 2), (4, 2, 5, 1),
                            (1, 5, 2, 5), (5, 7, 7, 1), (2, 5, 6, 3), (5, 6, 6, 3), (5, 6, 9, 1), (7, 8, 9, 1), (5, 8, 10, 1)])
    @test temporal_betweenness(g, MinimumHops(); β=1) ≈ [0, 9.5, 2, 4, 10, 0, 0.5, 0]
    @test temporal_betweenness(g, EarliestArrival(); β=1) ≈ [0, 9.5, 2.5, 4.5, 10, 0, 3, 0]
    @test temporal_betweenness(g, Fastest(); β=1) ≈ [0, 10.5, 2, 4, 10, 0, 0, 0]
    @test temporal_betweenness(g, ShortestForemost(); β=1) ≈ [0, 9, 2, 4, 10, 0, 3, 0]
    @test temporal_betweenness(g, ShortestFastest(); β=1) ≈ [0, 10, 2, 4, 10, 0, 0, 0]
    # the latest rows of Table 2 do not follow the definition (see the documentation)
    @test temporal_betweenness(g, LatestDeparture(); β=1) ≈ oracle_temporal_betweenness(g, LatestDeparture(), 1)
    @test temporal_betweenness(g, ShortestLatest(); β=1) ≈ oracle_temporal_betweenness(g, ShortestLatest(), 1)

    criteria = (MinimumHops(), EarliestArrival(), LatestDeparture(), Fastest(),
                ShortestForemost(), ShortestLatest(), ShortestFastest())
    rng = MersenneTwister(37)
    for _ in 1:60
        n = rand(rng, 2:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:8), rand(rng, 1:3))
                                for _ in 1:rand(rng, 1:12)])
        for c in criteria, β in (Inf, 0, 1, 3)
            @test temporal_betweenness(g, c; β=β) ≈ oracle_temporal_betweenness(g, c, isfinite(β) ? β : 10^9) atol = 1e-9
        end
        @test temporal_betweenness(g, PrefixForemost()) ≈ oracle_prefix_foremost_betweenness(g) atol = 1e-9
        @test temporal_betweenness(g, MinimumHops(); count_type=Rational{BigInt}) ≈ temporal_betweenness(g)
    end

    # counts that overflow Float64: 2^1100 foremost walks through a chain of diamonds
    k = 1100
    es = TemporalEdge{Int32,Int64}[]
    for i in 0:k-1  # diamond from node 3i + 1 to node 3i + 4 through 3i + 2 and 3i + 3
        push!(es, TemporalEdge(3i + 1, 3i + 2, 2i, 1), TemporalEdge(3i + 1, 3i + 3, 2i, 1),
                  TemporalEdge(3i + 2, 3i + 4, 2i + 1, 1), TemporalEdge(3i + 3, 3i + 4, 2i + 1, 1))
    end
    chain = OrderedEdgeList(3k + 1, es)
    exact = @test_logs (:warn, r"overflow") temporal_betweenness(chain, EarliestArrival())
    @test all(isfinite, exact)
    # all walks between the 3j nodes before and the 3(k - j) nodes after 3j + 1 pass through it
    @test all(exact[3j+1] ≈ 9j * (k - j) for j in (1, k ÷ 2, k - 1))
    @test_throws OverflowError temporal_betweenness(chain, EarliestArrival(); count_type=Float64)

    @test_throws ArgumentError temporal_betweenness(OrderedEdgeList(2, [(1, 2, 1, -1)]))
    @test_throws ArgumentError temporal_betweenness(g, PrefixForemost(); β=2)
    @test_throws ArgumentError temporal_distances(g, 1, ShortestForemost())
    @test_throws ArgumentError temporal_betweenness(g, MinimumTransitionTimes())
end

@testset "non-strict temporal betweenness and ATBC (Zhang et al.)" begin
    # transition times 0 give non-strict walks; cycles of such edges at one time stamp
    # give infinitely many optimal walks except for the shortest criteria
    criteria = (MinimumHops(), EarliestArrival(), LatestDeparture(), Fastest(),
                ShortestForemost(), ShortestLatest(), ShortestFastest())
    shortest = (MinimumHops, ShortestForemost, ShortestLatest, ShortestFastest)
    rng = MersenneTwister(71)
    errors = 0
    for _ in 1:60
        n = rand(rng, 2:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:6), rand(rng, 0:2))
                                for _ in 1:rand(rng, 1:10)])
        m = length(edges(g))
        for β in (Inf, 0, 2)
            b = isfinite(β) ? β : 10^9
            # a walk longer than m repeats an edge, so it goes around a cycle at one time stamp
            cycle = any(w -> length(w) > m, all_restless_walks(edges(g), b; maxlen=m + 1))
            for c in criteria
                if cycle && !any(S -> c isa S, shortest)
                    @test_throws ArgumentError temporal_betweenness(g, c; β=β)
                    errors += 1
                else
                    @test temporal_betweenness(g, c; β=β) ≈ oracle_temporal_betweenness(g, c, b; maxlen=m + 1) atol = 1e-9
                end
            end
        end
    end
    @test errors > 0
    cyc = OrderedEdgeList(3, [(1, 2, 1, 0), (2, 3, 1, 0), (3, 2, 1, 0)])
    @test_throws ArgumentError temporal_betweenness(cyc, EarliestArrival())
    @test temporal_betweenness(cyc, MinimumHops()) ≈ [0, 1, 0]  # 1 → 3 through 2
    @test_throws ArgumentError temporal_betweenness(cyc, PrefixForemost())

    # ATBC: normalized betweenness up to ε with probability 1 - δ, strict and non-strict
    strict = random_graph(MersenneTwister(43), 40, 1200)
    n = num_nodes(strict)
    nonstrict = OrderedEdgeList(n, [TemporalEdge(e.u, e.v, e.t, 0) for e in edges(strict)])
    for (g, cs) in ((strict, (MinimumHops(), EarliestArrival(), Fastest(), ShortestForemost())),
                    (nonstrict, (MinimumHops(), ShortestForemost())))
        for c in cs
            exact = temporal_betweenness(g, c) ./ (n * (n - 1))
            # foremost and fastest walks can revisit nodes: a warning says the bound may not hold
            r = @test_logs match_mode=:any temporal_betweenness_atbc(g, c; ε=0.02, rng=MersenneTwister(3))
            @test r.samples > 0 && r.error_bound <= 0.02
            @test maximum(abs.(r.estimates .- exact)) <= 0.02
        end
    end
    r = temporal_betweenness_atbc(strict; ε=0.02, max_samples=100, rng=MersenneTwister(4))
    @test r.samples <= 100 && r.error_bound > 0.02
    @test_throws ArgumentError temporal_betweenness_atbc(strict, LatestDeparture())
    @test_throws ArgumentError temporal_betweenness_atbc(strict; ε=0)
end

@testset "betweenness proxies (Becker, Crescenzi, Cruciani, Kodric)" begin
    rng = MersenneTwister(41)
    for _ in 1:40
        n = rand(rng, 2:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:8), rand(rng, 1:2))
                                for _ in 1:rand(rng, 1:12)])
        @test temporal_pass_through_degree(g) ≈ oracle_pass_through_degree(g)
        @test temporal_ego_betweenness(g, MinimumHops()) ≈ oracle_ego_betweenness(g, MinimumHops()) atol = 1e-9
        @test temporal_ego_betweenness(g, PrefixForemost()) ≈ oracle_ego_betweenness(g, PrefixForemost()) atol = 1e-9
    end
    g = OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1), (3, 1, 5, 1), (2, 1, 0, 1)])
    @test temporal_pass_through_degree(g) == [1.0, 1.0, 1.0]  # (2→1→2), (1→2→3), (2→3→1)
    @test temporal_pass_through_degree(g) ≈ oracle_pass_through_degree(g)

    # ONBRA: unbiased estimates of the normalized betweenness with a valid error bound
    g = random_graph(MersenneTwister(43), 30, 600)
    n = num_nodes(g)
    exact = temporal_betweenness(g) ./ (n * (n - 1))
    est, ε = temporal_betweenness_approximation(g, 20000; rng=MersenneTwister(1))
    @test maximum(abs.(est .- exact)) <= ε
    @test maximum(abs.(est .- exact)) < 0.01
    est2, _ = temporal_betweenness_approximation(g, 200, Fastest(); β=50, rng=MersenneTwister(2))
    @test length(est2) == n && all(>=(0), est2)
    @test_throws ArgumentError temporal_betweenness_approximation(g, 1)
end

@testset "MANTRA and temporal distance statistics (Cruciani)" begin
    TG = TemporalGraphs
    rng = MersenneTwister(59)
    criteria = (MinimumHops(), ShortestForemost(), PrefixForemost())
    for _ in 1:30
        n = rand(rng, 3:7)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:8), rand(rng, 1:3))
                                for _ in 1:rand(rng, 2:14)])
        for c in criteria
            # both estimators are unbiased: their means over the whole sample space are exact
            exact = temporal_betweenness(g, c) ./ (n * (n - 1))
            sm = TG._sampler(g, c)
            by_source = zeros(n)
            by_pair = zeros(n)
            for s in 1:n
                TG._sample!(sm, s, 0)
                by_source .+= sm.x ./ n
                for z in 1:n
                    z == s && continue
                    TG._sample!(sm, s, z)
                    by_pair .+= sm.x ./ (n * (n - 1))
                end
            end
            @test by_source ≈ exact atol = 1e-12
            @test by_pair ≈ exact atol = 1e-12
            # distance statistics
            d = oracle_hop_distances(g, c)
            pairs = sort!([d[s, t] for s in 1:n, t in 1:n if s != t && d[s, t] < typemax(Int)])
            st = temporal_distance_statistics(g, c; τ=0.5)
            if isempty(pairs)
                @test st.connectivity_rate == 0
            else
                @test st.diameter == pairs[end]
                @test st.effective_diameter == pairs[ceil(Int, length(pairs) / 2)]
                @test st.connectivity_rate ≈ length(pairs) / (n * (n - 1))
                @test st.average_distance ≈ sum(pairs) / length(pairs)
            end
        end
    end
    g = random_graph(MersenneTwister(43), 40, 1200)
    n = num_nodes(g)
    for c in criteria, est in (:pairs, :sources)
        exact = temporal_betweenness(g, c) ./ (n * (n - 1))
        r = temporal_betweenness_mantra(g, c; ε=0.01, estimator=est, rng=MersenneTwister(1))
        @test r.samples > 0 && r.error_bound <= 0.01
        @test maximum(abs.(r.estimates .- exact)) <= 0.01
    end
    st = temporal_distance_statistics(g; samples=200, rng=MersenneTwister(2))
    ex = temporal_distance_statistics(g)
    @test ex.diameter == temporal_diameter(g, MinimumHops())
    @test abs(st.connectivity_rate - ex.connectivity_rate) < 0.1
    @test st.diameter <= ex.diameter
    @test_throws ArgumentError temporal_betweenness_mantra(g, Fastest())
    @test_throws ArgumentError temporal_betweenness_mantra(g; estimator=:paths)
end
