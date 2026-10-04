@inline function _closeness_sum(dist::AbstractVector{C}) where {C}
    c = 0.0
    @inbounds for d in dist
        (d > zero(C) && d < typemax(C)) && (c += 1.0 / d)
    end
    return c
end

"""
    temporal_closeness(g, s, dt::DistanceType, ti = time_interval(g))

Temporal closeness of node `s`: the sum of `1/d(s, v)` over all nodes `v` with
`0 < d(s, v) < ∞`, where `d` is the temporal distance of type `dt` (for
[`EarliestArrival`](@ref) the arrival time, for [`LatestDeparture`](@ref) `ti[2]`
minus the latest departure time).
"""
function temporal_closeness(g::_AnyTemporalGraph, s::Integer, dt::DistanceType, ti=time_interval(g))
    dist = Vector{distance_eltype(g, dt)}(undef, num_nodes(g))
    _metric_distances!(dist, distance_workspace(g, dt), g, _check_node(g, s), dt, ti)
    return _closeness_sum(dist)
end

# Apply `f(dist)` to the distances from every node in parallel (OhMyThreads);
# result[s] = f(distances from s). Every task has its own buffers (`@local`).
function _map_sources(f::F, ::Type{R}, g, dt::DistanceType, ti) where {F,R}
    n = num_nodes(g)
    result = Vector{R}(undef, n)
    D = distance_eltype(g, dt)
    @tasks for s in 1:n
        @local begin
            dist = Vector{D}(undef, n)
            ws = distance_workspace(g, dt)
        end
        _metric_distances!(dist, ws, g, s, dt, ti)
        @inbounds result[s] = f(dist)
    end
    return result
end

"""
    temporal_closeness(g, dt::DistanceType, ti = time_interval(g))

Temporal closeness of all nodes. Runs in parallel when Julia is started with
several threads.
"""
temporal_closeness(g::_AnyTemporalGraph, dt::DistanceType, ti=time_interval(g)) =
    _map_sources(_closeness_sum, Float64, g, dt, _interval(time_type(g), ti))

"""
    temporal_closeness_approximation(g, h, ti = time_interval(g); rng = Random.default_rng())

Estimate the fastest path temporal closeness of all nodes from `h` samples: for
uniformly random nodes `w`, the durations towards `w` are computed on the reversed
temporal graph, and the sums are scaled by `n / h`, which gives an unbiased estimate.
"""
function temporal_closeness_approximation(g::Union{OrderedEdgeList,IncidentLists}, h::Integer,
                                          ti=time_interval(g); rng::AbstractRNG=Random.default_rng())
    n = num_nodes(g)
    result = zeros(Float64, n)
    (n == 0 || h <= 0) && return result
    rg = reverse(to_ordered_edge_list(g), ti)
    ws = distance_workspace(rg, Fastest())
    dist = Vector{time_type(g)}(undef, n)
    for _ in 1:h
        temporal_distances!(dist, ws, rg, rand(rng, 1:n), Fastest())
        @inbounds for v in 1:n
            d = dist[v]
            (d > 0 && d < typemax(d)) && (result[v] += 1.0 / d)
        end
    end
    result .*= n / h
    return result
end
