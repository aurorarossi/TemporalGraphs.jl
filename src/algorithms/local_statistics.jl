function _clustering!(mark::Vector{Bool}, nbrs::Vector{Int}, g::IncidentLists{V,T}, u::Int, a::T, b::T) where {V,T}
    empty!(nbrs)
    @inbounds for e in out_edges(g, u)
        e.t > b && break
        if e.t >= a && !mark[e.v]
            mark[e.v] = true
            push!(nbrs, e.v)
        end
    end
    count = 0
    @inbounds for v in nbrs, e in out_edges(g, v)
        e.t > b && break
        e.t >= a && mark[e.v] && (count += 1)
    end
    for v in nbrs
        mark[v] = false
    end
    k = length(nbrs)
    k < 2 && return 0.0
    return (1.0 / (b - a)) * (count / (k * (k - 1)))
end

"""
    temporal_clustering_coefficient(g::IncidentLists, s, ti = time_interval(g))

Temporal clustering coefficient of `s`: the number of temporal edges inside `ti`
between the out-neighbors of `s` divided by `k(k - 1)` and by the length of `ti`,
where `k` is the number of out-neighbors of `s` inside `ti` (0 if `k < 2`).
"""
function temporal_clustering_coefficient(g::IncidentLists{V,T}, s::Integer, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    return _clustering!(zeros(Bool, num_nodes(g)), Int[], g, _check_node(g, s), a, b)
end

"""
    temporal_clustering_coefficient(g::IncidentLists, ti = time_interval(g))

Temporal clustering coefficient of all nodes.
"""
function temporal_clustering_coefficient(g::IncidentLists{V,T}, ti::Tuple=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    mark = zeros(Bool, num_nodes(g))
    nbrs = Int[]
    return [_clustering!(mark, nbrs, g, u, a, b) for u in 1:num_nodes(g)]
end

"""
    topological_overlap(g::IncidentLists, s, ti = time_interval(g))

Average topological overlap of the out-neighborhoods of `s` at consecutive time
stamps inside `ti`: the mean over consecutive pairs of neighborhoods `(A, B)` of
`|A ∩ B| / sqrt(|A| |B|)`. Returns 0 if `s` has edges at fewer than two time stamps.
"""
function topological_overlap(g::IncidentLists{V,T}, s::Integer, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    groups = Vector{Vector{Int}}()
    last_t = typemin(T)
    for e in out_edges(g, _check_node(g, s))
        (e.t < a || e.t > b) && continue
        if isempty(groups) || e.t != last_t
            push!(groups, Int[])
            last_t = e.t
        end
        push!(groups[end], e.v)
    end
    length(groups) < 2 && return 0.0
    foreach(x -> unique!(sort!(x)), groups)
    to = 0.0
    for i in 1:length(groups)-1
        A, B = groups[i], groups[i+1]
        to += length(intersect(A, B)) / sqrt(length(A) * length(B))
    end
    return to / (length(groups) - 1)
end

"""
    topological_overlap(g::IncidentLists, ti = time_interval(g))

Topological overlap of all nodes.
"""
topological_overlap(g::IncidentLists, ti::Tuple=time_interval(g)) =
    [topological_overlap(g, u, ti) for u in 1:num_nodes(g)]

for f in (:temporal_clustering_coefficient, :topological_overlap)
    @eval $f(g::OrderedEdgeList, args...) = $f(to_incident_lists(g), args...)
end
