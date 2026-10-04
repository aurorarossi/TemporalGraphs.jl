# Conversions between TemporalGraphs.jl and GNNGraphs.jl (GraphNeuralNetworks.jl).
module TemporalGraphsGNNGraphsExt

using TemporalGraphs: TemporalGraphs, OrderedEdgeList, TemporalEdge, edges, num_nodes, time_interval, snapshots,
                      augmented_event_graph
using GNNGraphs: GNNGraphs, GNNGraph, TemporalSnapshotsGNNGraph, edge_index

# Continuous-time representation: all temporal edges, with their times as edge features.
function GNNGraphs.GNNGraph(g::OrderedEdgeList)
    es = edges(g)
    return GNNGraph([Int(e.u) for e in es], [Int(e.v) for e in es]; num_nodes=num_nodes(g),
                    edata=(t=[e.t for e in es], tt=[e.tt for e in es]))
end

function TemporalGraphs.OrderedEdgeList(h::GNNGraph; t::Symbol=:t, tt::Symbol=:tt, transition_time::Real=1,
                                        node_type::Type{<:Integer}=Int32)
    haskey(h.edata, t) || throw(ArgumentError("the GNNGraph has no edge feature $t with the times"))
    s, d = edge_index(h)
    ts = vec(h.edata[t])
    tts = haskey(h.edata, tt) ? vec(h.edata[tt]) : fill(transition_time, length(s))
    T = promote_type(eltype(ts), eltype(tts))
    es = [TemporalEdge{node_type,T}(s[i], d[i], ts[i], tts[i]) for i in eachindex(s)]
    return OrderedEdgeList(h.num_nodes, es)  # sorted by time only if needed: GNNGraph(g) round trips
end

# Discrete-time representation: one GNNGraph per snapshot, the times in tgdata.
function GNNGraphs.TemporalSnapshotsGNNGraph(g::OrderedEdgeList, ti=nothing; resolution=nothing)
    S = snapshots(g, ti; resolution=resolution)
    n = num_nodes(g)
    gs = map(S.graphs) do h
        src = [e.src for e in TemporalGraphs.Graphs.edges(h)]
        dst = [e.dst for e in TemporalGraphs.Graphs.edges(h)]
        GNNGraph(src, dst; num_nodes=n)
    end
    tg = TemporalSnapshotsGNNGraph(gs)
    tg.tgdata.times = S.times
    return tg
end

function TemporalGraphs.OrderedEdgeList(tg::TemporalSnapshotsGNNGraph,
                                        times=(:times in keys(tg.tgdata) ? tg.tgdata.times : 1:tg.num_snapshots);
                                        transition_time::Real=1, node_type::Type{<:Integer}=Int32)
    length(times) == tg.num_snapshots || throw(ArgumentError("one time is needed per snapshot"))
    n = maximum(tg.num_nodes; init=0)
    T = promote_type(eltype(times), typeof(transition_time))
    es = TemporalEdge{node_type,T}[]
    for (k, h) in enumerate(tg.snapshots)
        s, d = edge_index(h)
        for i in eachindex(s)
            push!(es, TemporalEdge{node_type,T}(s[i], d[i], times[k], transition_time))
        end
    end
    isempty(times) && return OrderedEdgeList(n, es)
    return OrderedEdgeList(n, sort!(es), (convert(T, first(times)), convert(T, last(times)) + convert(T, transition_time)))
end

function TemporalGraphs.augmented_event_gnngraph(g::OrderedEdgeList, ti=time_interval(g); β::Real=Inf)
    A = augmented_event_graph(g, ti; β=β)
    src = [e.src for e in TemporalGraphs.Graphs.edges(A.graph)]
    dst = [e.dst for e in TemporalGraphs.Graphs.edges(A.graph)]
    x = zeros(Float32, 2, length(A.labels))
    for (i, l) in enumerate(A.labels)
        x[l+1, i] = 1
    end
    return GNNGraph(src, dst; num_nodes=length(A.labels), ndata=(label=A.labels, x=x))
end

end
