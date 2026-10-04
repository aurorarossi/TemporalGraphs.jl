# Conversions between TemporalGraphs.jl and the graph containers of MLDatasets.jl.
module TemporalGraphsMLDatasetsExt

using TemporalGraphs: TemporalGraphs, OrderedEdgeList, TemporalEdge, edges, num_nodes, time_interval, snapshots
using MLDatasets: MLDatasets

# a sequence of snapshots, e.g. the brain networks of MLDatasets.TemporalBrains()
function TemporalGraphs.OrderedEdgeList(tg::MLDatasets.TemporalSnapshotsGraph, times=1:tg.num_snapshots;
                                        transition_time::Real=1, node_type::Type{<:Integer}=Int32)
    length(times) == tg.num_snapshots || throw(ArgumentError("one time is needed per snapshot"))
    n = maximum(tg.num_nodes; init=0)
    T = promote_type(eltype(times), typeof(transition_time))
    es = TemporalEdge{node_type,T}[]
    for (k, s) in enumerate(tg.snapshots)
        us, vs = s.edge_index
        for i in eachindex(us)
            push!(es, TemporalEdge{node_type,T}(us[i], vs[i], times[k], transition_time))
        end
    end
    isempty(times) && return OrderedEdgeList(n, es)
    return OrderedEdgeList(n, sort!(es), (convert(T, first(times)), convert(T, last(times)) + convert(T, transition_time)))
end

function MLDatasets.TemporalSnapshotsGraph(g::OrderedEdgeList, ti=nothing; resolution=nothing)
    S = snapshots(g, ti; resolution=resolution)
    n = num_nodes(g)
    graphs = map(S.graphs) do h
        us = [e.src for e in TemporalGraphs.Graphs.edges(h)]
        vs = [e.dst for e in TemporalGraphs.Graphs.edges(h)]
        MLDatasets.Graph(; num_nodes=n, edge_index=(us, vs))
    end
    K = length(graphs)
    return MLDatasets.TemporalSnapshotsGraph(fill(n, K), [s.num_edges for s in graphs], K, graphs, (times=S.times,))
end

_time_feature(data, time) = haskey(data, time) ? data[time] :
                            haskey(data, String(time)) ? data[String(time)] :
                            throw(ArgumentError("no edge feature $time with the times"))

# a static graph with the times of its edges as edge features
function TemporalGraphs.OrderedEdgeList(g::MLDatasets.Graph; time::Symbol=:timestamp, transition_time::Real=1,
                                        node_type::Type{<:Integer}=Int32)
    g.edge_data === nothing && throw(ArgumentError("the graph has no edge data"))
    ts = vec(_time_feature(g.edge_data, time))
    us, vs = g.edge_index
    T = promote_type(eltype(ts), typeof(transition_time))
    return OrderedEdgeList(g.num_nodes, sort!([TemporalEdge{node_type,T}(us[i], vs[i], ts[i], transition_time) for i in eachindex(us)]))
end

# one relation of a heterogeneous graph with timestamps, e.g. the user-movie ratings of
# MLDatasets.MovieLens: the source nodes are numbered first, then the target nodes
function TemporalGraphs.OrderedEdgeList(h::MLDatasets.HeteroGraph, relation::Tuple{String,String,String};
                                        time::Symbol=:timestamp, transition_time::Real=1, node_type::Type{<:Integer}=Int32)
    haskey(h.edge_indices, relation) || throw(ArgumentError("the graph has no relation $relation"))
    us, vs = h.edge_indices[relation]
    ts = vec(_time_feature(h.edge_data[relation], time))
    src, _, dst = relation
    offset = src == dst ? 0 : h.num_nodes[src]
    n = src == dst ? h.num_nodes[src] : h.num_nodes[src] + h.num_nodes[dst]
    T = promote_type(eltype(ts), typeof(transition_time))
    es = [TemporalEdge{node_type,T}(us[i], vs[i] + offset, ts[i], transition_time) for i in eachindex(us)]
    return OrderedEdgeList(n, sort!(es))
end

end
