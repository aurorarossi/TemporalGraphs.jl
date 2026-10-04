# Temporal distances, optimal paths, diameter and efficiency.

@testset "distances: expected values" begin
    g = paper_file()
    il = to_incident_lists(g)
    trs = to_trs_graph(g)
    for x in (g, il, trs)
        @test minimum_durations(x, 1) == [0, 1, 4, 7]
        @test earliest_arrival_times(x, 1) == [0, 3, 6, 8]
        @test minimum_transition_times(x, 1) == [0, 1, 3, 7]
    end
    @test [number_of_reachable_nodes(g, v) for v in 1:4] == [4, 2, 2, 3]
    @test earliest_arrival_times(chain(), 1) == [0, 2, 3, 4, 5, 6, 7, 8]
    @test latest_departure_times(g, 1) == [12, 5, 5, 1]
    @test temporal_distances(g, 1, LatestDeparture()) == latest_departure_times(g, 1)
    @test latest_departure_times(il, 1) == latest_departure_times(g, 1)
    @test_throws ArgumentError temporal_distances(to_trs_graph(g), 1, LatestDeparture())
    # restricted interval
    @test earliest_arrival_times(g, 1, (2, 9)) == [0, 3, 9, INF]
    @test minimum_durations(g, 1, (3, 20)) == [0, 2, 4, INF]
    @test minimum_hops(g, 1, (1, 5)) == [0, 1, INF, INF]
end

@testset "distances: all representations agree with brute force" begin
    rng = MersenneTwister(42)
    for trial in 1:120
        n = rand(rng, 2:8)
        m = rand(rng, 0:40)
        g = random_general_graph(rng, n, m; tmax=rand(rng, 1:20), ttmax=rand(rng, 0:4))
        il = to_incident_lists(g)
        trs = to_trs_graph(g)
        a, b = time_interval(g)
        intervals = [time_interval(g), (a + rand(rng, 0:5), b - rand(rng, 0:5))]
        for ti in intervals, s in 1:n
            o = oracle_distances(g, s, ti)
            for x in (g, il, trs)
                @test earliest_arrival_times(x, s, ti) == o.ea
                @test minimum_durations(x, s, ti) == o.fastest
                @test minimum_hops(x, s, ti) == o.hops
                @test minimum_transition_times(x, s, ti) == o.mtt
            end
            @test latest_departure_times(g, s, ti) == o.ld
            @test latest_departure_times(il, s, ti) == o.ld
            @test number_of_reachable_nodes(g, s, ti) == count(<(INF), o.ea)
        end
    end
end

@testset "distances: larger random graphs" begin
    rng = MersenneTwister(7)
    for (n, m) in ((10, 1000), (8, 40), (10, 100), (50, 2000))
        g = random_graph(rng, n, m)
        il = to_incident_lists(g)
        trs = to_trs_graph(g)
        for s in 1:n, f in (earliest_arrival_times, minimum_durations, minimum_hops, minimum_transition_times)
            d = f(g, s)
            @test f(il, s) == d
            @test f(trs, s) == d
        end
    end
end

@testset "paths" begin
    g = paper_file()
    il = to_incident_lists(g)
    @test minimum_duration_path(il, 1, 3) == [TemporalEdge(1, 2, 5, 2), TemporalEdge(2, 3, 7, 2)]
    @test earliest_arrival_path(il, 1, 3) == [TemporalEdge(1, 3, 1, 5)]
    @test minimum_transition_time_path(il, 1, 3) == [TemporalEdge(1, 2, 2, 1), TemporalEdge(2, 3, 7, 2)]
    @test minimum_hops_path(g, 1, 3) == [TemporalEdge(1, 3, 1, 5)]
    @test isempty(minimum_hops_path(il, 2, 1))
    @test isempty(earliest_arrival_path(il, 1, 1))

    # a returned path is a valid temporal path with the optimal value
    function valid(p, s, t, ti)
        isempty(p) && return false
        p[1].u == s && p[end].v == t || return false
        all(e -> ti[1] <= e.t && e.t + e.tt <= ti[2], p) || return false
        return all(p[i+1].u == p[i].v && p[i+1].t >= p[i].t + p[i].tt for i in 1:length(p)-1)
    end
    rng = MersenneTwister(3)
    for _ in 1:40
        n = rand(rng, 2:8)
        g = random_general_graph(rng, n, rand(rng, 1:40); min_tt=0)
        il = to_incident_lists(g)
        ti = time_interval(g)
        for s in 1:n
            o = oracle_distances(g, s, ti)
            for t in 1:n
                t == s && continue
                p1 = earliest_arrival_path(il, s, t)
                p2 = minimum_duration_path(il, s, t)
                p3 = minimum_transition_time_path(il, s, t)
                p4 = minimum_hops_path(il, s, t)
                if o.ea[t] == INF
                    @test isempty(p1) && isempty(p2) && isempty(p3) && isempty(p4)
                else
                    @test valid(p1, s, t, ti) && p1[end].t + p1[end].tt == o.ea[t]
                    @test valid(p2, s, t, ti) && p2[end].t + p2[end].tt - p2[1].t == o.fastest[t]
                    @test valid(p3, s, t, ti) && sum(e.tt for e in p3) == o.mtt[t]
                    @test valid(p4, s, t, ti) && length(p4) == o.hops[t]
                end
            end
        end
    end
end

@testset "diameter and efficiency" begin
    g = paper_file()
    @test temporal_efficiency(g, Fastest()) ≈ 0.310515873015873
    @test temporal_diameter(g, MinimumHops()) == 2
    @test temporal_eccentricity(g, 1, MinimumHops()) == 2
    @test temporal_eccentricity(g, 2, MinimumHops()) == 1
    c = chain()
    @test temporal_diameter(c, Fastest()) == temporal_diameter(to_incident_lists(c), Fastest()) == 7
    @test temporal_efficiency(c, Fastest()) == temporal_efficiency(to_incident_lists(c), Fastest())
    @test temporal_diameter(to_trs_graph(c), MinimumHops()) == 7
end
