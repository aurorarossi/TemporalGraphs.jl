function _label_path(g::IncidentLists{V,T}, ws::IncidentListsWorkspace, l::Int) where {V,T}
    path = TemporalEdge{V,T}[]
    @inbounds while l != 0 && ws.edge[l] != 0
        push!(path, g.edges[ws.edge[l]])
        l = ws.parent[l]
    end
    return reverse!(path)
end

function _optimal_path(g::IncidentLists{V,T}, s::Integer, target::Integer, dt::DistanceType, ti) where {V,T}
    s = _check_node(g, s)
    target = _check_node(g, target)
    ws = distance_workspace(g, dt)
    dist = Vector{distance_eltype(g, dt)}(undef, num_nodes(g))
    l = _il_labels!(dist, ws, g, s, _interval(T, ti), dt, target)
    return _label_path(g, ws, l)
end

"""
    earliest_arrival_path(g, s, target, ti = time_interval(g))

An earliest arrival path from `s` to `target` inside `ti` as a vector of temporal
edges; empty if `target` is not reachable or equal to `s`.
"""
function earliest_arrival_path(g::IncidentLists{V,T}, s::Integer, target::Integer, ti=time_interval(g)) where {V,T}
    s = _check_node(g, s)
    target = _check_node(g, target)
    ws = EarliestArrivalWorkspace{T}(num_nodes(g))
    _il_ea!(Vector{T}(undef, num_nodes(g)), ws, g, s, _interval(T, ti), target)
    path = TemporalEdge{V,T}[]
    (target == s || !ws.done[target]) && return path
    v = target
    while v != s
        e = g.edges[ws.pred[v]]
        push!(path, e)
        v = Int(e.u)
    end
    return reverse!(path)
end

"""
    minimum_duration_path(g, s, target, ti = time_interval(g))

A fastest path from `s` to `target` inside `ti`; empty if there is none.
"""
minimum_duration_path(g::IncidentLists, s::Integer, target::Integer, ti=time_interval(g)) =
    _optimal_path(g, s, target, Fastest(), ti)

"""
    minimum_transition_time_path(g, s, target, ti = time_interval(g))

A path with minimum sum of transition times from `s` to `target` inside `ti`; empty
if there is none.
"""
minimum_transition_time_path(g::IncidentLists, s::Integer, target::Integer, ti=time_interval(g)) =
    _optimal_path(g, s, target, MinimumTransitionTimes(), ti)

"""
    minimum_hops_path(g, s, target, ti = time_interval(g))

A path with the minimum number of edges from `s` to `target` inside `ti`; empty if
there is none.
"""
minimum_hops_path(g::IncidentLists, s::Integer, target::Integer, ti=time_interval(g)) =
    _optimal_path(g, s, target, MinimumHops(), ti)

for f in (:earliest_arrival_path, :minimum_duration_path, :minimum_transition_time_path, :minimum_hops_path)
    @eval $f(g::OrderedEdgeList, s::Integer, target::Integer, ti=time_interval(g)) =
        $f(to_incident_lists(g), s, target, ti)
end
