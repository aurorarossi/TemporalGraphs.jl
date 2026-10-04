const _AnyTemporalGraph = Union{OrderedEdgeList,IncidentLists,TRSGraph}

# criteria that only define optimal walks for the betweenness
const _BetweennessOnly = Union{ShortestForemost,ShortestFastest,ShortestLatest,PrefixForemost}
_betweenness_only(dt) = throw(ArgumentError("$(typeof(dt)) is only supported by temporal_betweenness"))
distance_workspace(::OrderedEdgeList, dt::_BetweennessOnly) = _betweenness_only(dt)
distance_workspace(::IncidentLists, dt::_BetweennessOnly) = _betweenness_only(dt)
distance_workspace(::TRSGraph, dt::_BetweennessOnly) = _betweenness_only(dt)

"""
    temporal_distances(g, s, dt::DistanceType, ti = time_interval(g); β = Inf)

Temporal distances of type `dt` from `s` to all nodes of `g` inside the time
interval `ti`, for an [`OrderedEdgeList`](@ref), [`IncidentLists`](@ref) or
[`TRSGraph`](@ref). Unreachable nodes get `typemax` of the element type (for
[`LatestDeparture`](@ref) `typemin`).

With a finite `β` (only for an `OrderedEdgeList`) the distances are those of
*β-restless walks* (Casteigts, Himmel, Molter and Zschoche, 2021), which wait at most
`β` at every intermediate node: every edge departs in `[a, a + β]`, where `a` is the
arrival time of the previous edge, and the departure from `s` is free. Such walks can
visit a node several times (finding restless *paths* is NP-hard, see
[`restless_path`](@ref)). They are computed in `O(m log m)` time with one scan of the
edges and a sliding window of the arrivals at every node.

For repeated computations use the in-place [`temporal_distances!`](@ref) with a
[`distance_workspace`](@ref), which avoids all allocations.

# References

- H. Wu, J. Cheng, Y. Ke, S. Huang, Y. Huang, and H. Wu. *Efficient algorithms for temporal path computation.* IEEE Transactions on Knowledge and Data Engineering 28(11), 2016. [DOI](https://doi.org/10.1109/TKDE.2016.2594065)
- A. Casteigts, A.-S. Himmel, H. Molter, and P. Zschoche. *Finding temporal paths under waiting time constraints.* Algorithmica 83(9), 2021. [DOI](https://doi.org/10.1007/s00453-021-00831-w), [arXiv](https://arxiv.org/abs/1909.06437)
"""
function temporal_distances(g::_AnyTemporalGraph, s::Integer, dt::DistanceType, ti=time_interval(g); β::Real=Inf)
    if isfinite(β)
        g isa OrderedEdgeList || throw(ArgumentError("waiting constraints need an OrderedEdgeList"))
        dt isa _RestlessCriterion || throw(ArgumentError("$(typeof(dt)) does not support waiting constraints"))
        β >= 0 || throw(ArgumentError("β must be non-negative"))
        return _restless_distances(g, _check_node(g, s), dt, β, ti)::Vector{distance_eltype(g, dt)}
    end
    dist = Vector{distance_eltype(g, dt)}(undef, num_nodes(g))
    return temporal_distances!(dist, distance_workspace(g, dt), g, s, dt, ti)
end

"""
    earliest_arrival_times(g, s, ti = time_interval(g))

Earliest arrival times at all nodes for paths starting at `s` inside `ti`.
Unreachable nodes get `typemax`, the source gets 0.
"""
earliest_arrival_times(g::_AnyTemporalGraph, s::Integer, ti=time_interval(g)) =
    temporal_distances(g, s, EarliestArrival(), ti)

"""
    minimum_durations(g, s, ti = time_interval(g))

Durations of fastest paths from `s` to all nodes inside `ti` (arrival time minus
departure time). Unreachable nodes get `typemax`, the source gets 0.
"""
minimum_durations(g::_AnyTemporalGraph, s::Integer, ti=time_interval(g)) =
    temporal_distances(g, s, Fastest(), ti)

"""
    latest_departure_times(g, s, ti = time_interval(g))

For every node `v`, the latest time at which one can leave `s` and still reach `v`
inside `ti`. Unreachable nodes get `typemin`, the source gets `ti[2]`.
Implemented for `OrderedEdgeList` and `IncidentLists`.
"""
latest_departure_times(g::_AnyTemporalGraph, s::Integer, ti=time_interval(g)) =
    temporal_distances(g, s, LatestDeparture(), ti)

"""
    minimum_hops(g, s, ti = time_interval(g))

Minimum number of edges of temporal paths from `s` to all nodes inside `ti`
(`typemax(Int)` for unreachable nodes).
"""
minimum_hops(g::_AnyTemporalGraph, s::Integer, ti=time_interval(g)) =
    temporal_distances(g, s, MinimumHops(), ti)

"""
    minimum_transition_times(g, s, ti = time_interval(g))

Minimum sum of transition times of temporal paths (shortest paths) from `s` to all
nodes inside `ti`. Unreachable nodes get `typemax`, the source gets 0.
"""
minimum_transition_times(g::_AnyTemporalGraph, s::Integer, ti=time_interval(g)) =
    temporal_distances(g, s, MinimumTransitionTimes(), ti)

# Distances used by closeness, eccentricity and efficiency, written into `dist`.
# For latest departure paths the distance of v is ti[2] minus the latest departure
# time towards v. Unreachable nodes get typemax.
function _metric_distances!(dist::AbstractVector{C}, ws, g, s::Int, dt::DistanceType, ti) where {C}
    temporal_distances!(dist, ws, g, s, dt, ti)
    if dt isa LatestDeparture
        b = convert(C, ti[2])
        @inbounds for i in eachindex(dist)
            dist[i] = dist[i] == typemin(C) ? typemax(C) : b - dist[i]
        end
    end
    return dist
end
