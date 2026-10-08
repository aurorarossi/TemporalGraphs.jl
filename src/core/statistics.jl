"""
    TemporalGraphStatistics

Basic statistics of a temporal graph, see [`get_statistics`](@ref).
"""
struct TemporalGraphStatistics{T<:Real}
    number_of_nodes::Int
    number_of_edges::Int
    number_of_static_edges::Int
    number_of_time_stamps::Int
    number_of_transition_times::Int
    minimal_time_stamp::T
    maximal_time_stamp::T
    minimal_transition_time::T
    maximal_transition_time::T
    min_temporal_in_degree::Int
    max_temporal_in_degree::Int
    min_temporal_out_degree::Int
    max_temporal_out_degree::Int
    temporality::Int
    max_edges_per_time_stamp::Int
    vertex_interval_membership_width::Int
    edge_interval_membership_width::Int
end

"""
    get_statistics(g::OrderedEdgeList)

Compute the basic statistics of `g`: numbers of nodes, temporal edges, static edges
(ordered pairs `(u, v)` with an edge) and distinct time stamps (the *lifetime*) and
transition times, their ranges, the ranges of the temporal in- and out-degrees, and the
temporal graph parameters:

- the *temporality*, the largest number of distinct time stamps of a static edge
  (Mertzios, Michail and Spirakis, 2019);
- the largest number of temporal edges with the same time stamp;
- the [`vertex_interval_membership_width`](@ref) and the
  [`edge_interval_membership_width`](@ref) of the static edges.

# References

- G. B. Mertzios, O. Michail, and P. G. Spirakis. *Temporal network optimization subject to connectivity constraints.* Algorithmica 81(4), 2019. [DOI](https://doi.org/10.1007/s00453-018-0478-6), [arXiv](https://arxiv.org/abs/1502.04382)
"""
function get_statistics(g::OrderedEdgeList{V,T}) where {V,T}
    n = num_nodes(g)
    es = edges(g)
    indeg = zeros(Int, n)
    outdeg = zeros(Int, n)
    for e in es
        indeg[e.v] += 1
        outdeg[e.u] += 1
    end
    ts = unique!(sort!([e.t for e in es]))
    tts = unique!(sort!([e.tt for e in es]))
    static = unique!(sort!([(e.u, e.v) for e in es]))
    labels = unique!(sort!([(e.u, e.v, e.t) for e in es]))
    ext(v, f, default) = isempty(v) ? default : f(v)
    return TemporalGraphStatistics{T}(
        n, length(es), length(static), length(ts), length(tts),
        ext(ts, first, typemax(T)), ext(ts, last, zero(T)), ext(tts, first, typemax(T)), ext(tts, last, zero(T)),
        ext(indeg, minimum, typemax(Int)), ext(indeg, maximum, 0),
        ext(outdeg, minimum, typemax(Int)), ext(outdeg, maximum, 0),
        _max_run(labels, x -> (x[1], x[2])), _max_run(es, e -> e.t),
        vertex_interval_membership_width(g), edge_interval_membership_width(g))
end

# Length of the longest run of consecutive elements of the sorted vector v with the
# same key.
function _max_run(v, key)
    best = 0
    i = 1
    while i <= length(v)
        j = i
        while j < length(v) && key(v[j+1]) == key(v[i])
            j += 1
        end
        best = max(best, j - i + 1)
        i = j + 1
    end
    return best
end

function Base.show(io::IO, ::MIME"text/plain", s::TemporalGraphStatistics)
    println(io, "number of nodes: ", s.number_of_nodes)
    println(io, "number of edges: ", s.number_of_edges)
    println(io, "number of static edges: ", s.number_of_static_edges)
    println(io, "number of time stamps: ", s.number_of_time_stamps)
    println(io, "number of transition times: ", s.number_of_transition_times)
    println(io, "min. time stamp: ", s.minimal_time_stamp)
    println(io, "max. time stamp: ", s.maximal_time_stamp)
    println(io, "min. transition time: ", s.minimal_transition_time)
    println(io, "max. transition time: ", s.maximal_transition_time)
    println(io, "min. temporal in-degree: ", s.min_temporal_in_degree)
    println(io, "max. temporal in-degree: ", s.max_temporal_in_degree)
    println(io, "min. temporal out-degree: ", s.min_temporal_out_degree)
    println(io, "max. temporal out-degree: ", s.max_temporal_out_degree)
    println(io, "temporality: ", s.temporality)
    println(io, "max. edges per time stamp: ", s.max_edges_per_time_stamp)
    println(io, "vertex-interval-membership width: ", s.vertex_interval_membership_width)
    print(io, "edge-interval-membership width: ", s.edge_interval_membership_width)
end

# The static edge of e: the ordered pair (u, v), or the unordered one.
_static_pair(e::TemporalEdge, directed::Bool) = directed ? (Int(e.u), Int(e.v)) : minmax(Int(e.u), Int(e.v))

# Largest number of the closed intervals [lo, hi] that contain a common point.
function _max_overlap(intervals::Vector{Tuple{T,T}}) where {T}
    events = Tuple{T,Int}[]
    for (lo, hi) in intervals
        push!(events, (lo, -1), (hi, 1))  # at equal times, openings come first
    end
    sort!(events)
    best = cur = 0
    for (_, x) in events
        cur -= x
        best = max(best, cur)
    end
    return best
end

# Interval [first time stamp, last time stamp] of the edges of every key (nodes or
# static edges); `keys(e)` returns the keys of type K of the edge e.
function _activity_intervals(g::OrderedEdgeList{V,T}, ::Type{K}, keys) where {V,T,K}
    span = Dict{K,Tuple{T,T}}()
    for e in edges(g), k in keys(e)
        lo, hi = get(span, k, (e.t, e.t))
        span[k] = (min(lo, e.t), max(hi, e.t))
    end
    return collect(values(span))
end

"""
    vertex_interval_membership_width(g::OrderedEdgeList)

The *vertex-interval-membership width* of `g` (Bumpus and Meeks, 2023): every node is
active from the first to the last time stamp of its edges, and the width is the
largest number of nodes active at the same time. Problems such as temporal Eulerian
walks, NP-hard in general, are fixed-parameter tractable with respect to it.
`O(m + n log n)` time.

# References

- B. M. Bumpus and K. Meeks. *Edge exploration of temporal graphs.* Algorithmica 85, 2023. [DOI](https://doi.org/10.1007/s00453-022-01018-7), [arXiv](https://arxiv.org/abs/2103.05387)
"""
vertex_interval_membership_width(g::OrderedEdgeList) =
    _max_overlap(_activity_intervals(g, Int, e -> e.u == e.v ? (Int(e.u),) : (Int(e.u), Int(e.v))))

"""
    edge_interval_membership_width(g::OrderedEdgeList; directed = true)

The *edge-interval-membership width* of `g` (Bumpus and Meeks, 2023): every static
edge is active from its first to its last time stamp, and the width is the largest
number of static edges active at the same time. The static edges are the ordered pairs
`(u, v)` with an edge, or the unordered pairs if `directed = false` (for undirected
graphs stored with both directions). See [`vertex_interval_membership_width`](@ref).
"""
edge_interval_membership_width(g::OrderedEdgeList; directed::Bool=true) =
    _max_overlap(_activity_intervals(g, Tuple{Int,Int}, e -> (_static_pair(e, directed),)))

"""
    is_simple(g::OrderedEdgeList; directed = true)

Whether the labeling of `g` is *simple*: every static edge has at most one distinct
time stamp. The static edges are the ordered pairs `(u, v)`, or the unordered pairs if
`directed = false` (for undirected graphs stored with both directions, which then
count once). Together with [`is_proper`](@ref), directedness and strictness (positive
transition times) these are the settings of temporal graphs compared by Casteigts,
Corsini and Sarkar (2024).

# References

- A. Casteigts, T. Corsini, and W. Sarkar. *Simple, strict, proper, happy: a study of reachability in temporal graphs.* Theoretical Computer Science 991, 2024. [DOI](https://doi.org/10.1016/j.tcs.2024.114434), [arXiv](https://arxiv.org/abs/2208.01720)
"""
function is_simple(g::OrderedEdgeList; directed::Bool=true)
    labels = unique!(sort!([(_static_pair(e, directed), e.t) for e in edges(g)]))
    return _max_run(labels, first) <= 1
end

"""
    is_proper(g::OrderedEdgeList; directed = true)

Whether the labeling of `g` is *proper*: adjacent static edges (static edges sharing a
node) never have the same time stamp, i.e., every snapshot is a matching. The static
edges are as in [`is_simple`](@ref); with `directed = true` the edges `(u, v, t)` and
`(v, u, t)` are adjacent, so undirected graphs stored with both directions need
`directed = false`.
"""
function is_proper(g::OrderedEdgeList{V,T}; directed::Bool=true) where {V,T}
    incidences = Tuple{Int,T,Tuple{Int,Int}}[]  # (node, time stamp, static edge)
    for e in edges(g)
        pair = _static_pair(e, directed)
        push!(incidences, (pair[1], e.t, pair))
        pair[1] == pair[2] || push!(incidences, (pair[2], e.t, pair))
    end
    unique!(sort!(incidences))
    return _max_run(incidences, x -> (x[1], x[2])) <= 1
end
