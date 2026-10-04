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
end

"""
    get_statistics(g::OrderedEdgeList)

Compute the basic statistics of `g`: numbers of nodes, temporal edges, static edges,
distinct time stamps and transition times, their ranges, and the ranges of the
temporal in- and out-degrees.
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
    ext(v, f, default) = isempty(v) ? default : f(v)
    return TemporalGraphStatistics{T}(
        n, length(es), length(static), length(ts), length(tts),
        ext(ts, first, typemax(T)), ext(ts, last, zero(T)), ext(tts, first, typemax(T)), ext(tts, last, zero(T)),
        ext(indeg, minimum, typemax(Int)), ext(indeg, maximum, 0),
        ext(outdeg, minimum, typemax(Int)), ext(outdeg, maximum, 0))
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
    print(io, "max. temporal out-degree: ", s.max_temporal_out_degree)
end
