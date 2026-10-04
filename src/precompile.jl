# Compile the main code paths for the default edge types, which reduces the latency
# of the first call (time to first result).
@setup_workload begin
    raw = [(1, 2, 1, 1), (1, 3, 2, 1), (2, 3, 3, 1), (3, 4, 5, 2), (4, 1, 7, 1), (2, 4, 8, 0)]
    @compile_workload begin
        for g in (OrderedEdgeList(4, raw), OrderedEdgeList(4, [(u, v, Float64(t), Float64(tt)) for (u, v, t, tt) in raw]))
            il = to_incident_lists(g)
            trs = to_trs_graph(g)
            for dt in (EarliestArrival(), Fastest(), LatestDeparture(), MinimumTransitionTimes(), MinimumHops())
                temporal_closeness(g, dt)
                temporal_distances(il, 1, dt)
                dt isa LatestDeparture || temporal_distances(trs, 1, dt)
            end
            compute_topk_closeness(il, 2, Fastest())
            minimum_duration_path(il, 1, 4)
            earliest_arrival_path(il, 1, 4)
            temporal_katz_centrality(g, 0.5)
            temporal_walk_centrality(g, 0.5, 0.5)
            temporal_edge_betweenness(g)
            get_statistics(g)
        end
    end
end
