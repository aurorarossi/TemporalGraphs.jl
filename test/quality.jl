# Type stability, allocations and code quality.

@testset "type stability and allocations" begin
    g = random_graph(MersenneTwister(31), 40, 800)
    gf = OrderedEdgeList(40, [TemporalEdge{Int32,Float64}(e.u, e.v, e.t, e.tt) for e in edges(g)])
    for x in (g, gf), y in (x, to_incident_lists(x), to_trs_graph(x))
        for dt in (EarliestArrival(), Fastest(), LatestDeparture(), MinimumTransitionTimes(), MinimumHops())
            (y isa TRSGraph && dt isa LatestDeparture) && continue
            @inferred temporal_distances(y, 1, dt)
            @inferred temporal_closeness(y, 1, dt)
            ws = distance_workspace(y, dt)
            dist = Vector{distance_eltype(y, dt)}(undef, num_nodes(y))
            for s in 1:num_nodes(y)  # let the buffers grow
                temporal_distances!(dist, ws, y, s, dt)
            end
            @test maximum(s -> @allocated(temporal_distances!(dist, ws, y, s, dt)), 1:num_nodes(y)) == 0
        end
    end
    @inferred temporal_walk_centrality(g, 0.5, 0.5)
    @inferred temporal_katz_centrality(g, 0.5)
    @inferred temporal_edge_betweenness(g)
    @inferred compute_topk_closeness(to_incident_lists(g), 3, Fastest())
    @inferred minimum_duration_path(to_incident_lists(g), 1, 2)
end

@testset "code quality (Aqua, JET)" begin
    Aqua.test_all(TemporalGraphs)
    JET.test_package(TemporalGraphs; target_modules=(TemporalGraphs,))
end
