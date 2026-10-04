@inline function _eccentricity(dist::AbstractVector{C}) where {C}
    m = zero(C)
    @inbounds for d in dist
        (d < typemax(C) && d > m) && (m = d)
    end
    return m
end

"""
    temporal_eccentricity(g, s, dt::DistanceType, ti = time_interval(g))

The largest finite temporal distance of type `dt` from `s` (0 if no other node is
reachable), with distances as in [`temporal_closeness`](@ref).
"""
function temporal_eccentricity(g::_AnyTemporalGraph, s::Integer, dt::DistanceType, ti=time_interval(g))
    dist = Vector{distance_eltype(g, dt)}(undef, num_nodes(g))
    _metric_distances!(dist, distance_workspace(g, dt), g, _check_node(g, s), dt, ti)
    return _eccentricity(dist)
end

"""
    temporal_diameter(g, dt::DistanceType, ti = time_interval(g))

The largest temporal eccentricity over all nodes. Runs in parallel when Julia is
started with several threads.
"""
function temporal_diameter(g::_AnyTemporalGraph, dt::DistanceType, ti=time_interval(g))
    C = distance_eltype(g, dt)
    ecc = _map_sources(_eccentricity, C, g, dt, _interval(time_type(g), ti))
    return isempty(ecc) ? zero(C) : maximum(ecc)
end

"""
    temporal_efficiency(g, dt::DistanceType, ti = time_interval(g))

Global temporal efficiency: the sum of the temporal closeness of all nodes divided
by `n(n - 1)`.
"""
function temporal_efficiency(g::_AnyTemporalGraph, dt::DistanceType, ti=time_interval(g))
    n = num_nodes(g)
    return sum(temporal_closeness(g, dt, ti)) / (n * (n - 1))
end
