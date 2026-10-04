"""
    DirectedLineGraph

Directed line graph of a temporal graph. Node `i` represents the temporal edge
`edges(dlg)[i]`, and there is an edge from `e = (u, v, t, tt)` to `f` iff `f` leaves
`v` and `f.t ≥ t + tt`, i.e., `(e, f)` is a temporal path of length two.
Construct it with [`to_directed_line_graph`](@ref).
"""
struct DirectedLineGraph{V<:Integer,T<:Real}
    num_tg_nodes::Int
    edges::Vector{TemporalEdge{V,T}}
    offsets::Vector{Int}
    heads::Vector{Int}
end

num_nodes(g::DirectedLineGraph) = length(g.edges)
num_edges(g::DirectedLineGraph) = length(g.heads)
edges(g::DirectedLineGraph) = g.edges

"""
    out_edges(dlg::DirectedLineGraph, i)

The out-neighbors of node `i` of the directed line graph.
"""
out_edges(g::DirectedLineGraph, i::Integer) = @inbounds view(g.heads, g.offsets[i]:g.offsets[i+1]-1)

# `perm[k]` is the DLG node of the k-th edge of `il`
function _dlg(il::IncidentLists{V,T}, es::Vector{TemporalEdge{V,T}}, perm::AbstractVector{Int}) where {V,T}
    m = length(es)
    first_succ = Vector{Int}(undef, length(il.edges))
    offsets = Vector{Int}(undef, m + 1)
    offsets[1] = 1
    inv = Vector{Int}(undef, m)
    for k in eachindex(perm)
        inv[perm[k]] = k
    end
    @inbounds for i in 1:m
        e = es[i]
        f = _first_edge_from(il, e.v, e.t + e.tt)
        first_succ[inv[i]] = f
        offsets[i+1] = offsets[i] + (il.offsets[e.v+1] - f)
    end
    heads = Vector{Int}(undef, offsets[end] - 1)
    @inbounds for i in 1:m
        e = es[i]
        p = offsets[i]
        for k in first_succ[inv[i]]:il.offsets[e.v+1]-1
            heads[p] = perm[k]
            p += 1
        end
    end
    return DirectedLineGraph{V,T}(num_nodes(il), es, offsets, heads)
end

"""
    to_directed_line_graph(g)

Directed line graph of `g`. For an `OrderedEdgeList`, node `i` of the result is the
`i`-th edge of `edges(g)`; for `IncidentLists` it is the `i`-th edge of `edges(g)` in
incident lists order.
"""
function to_directed_line_graph(g::IncidentLists)
    return _dlg(g, g.edges, 1:length(g.edges))
end

function to_directed_line_graph(g::OrderedEdgeList)
    il = to_incident_lists(g)
    # position in `il` of each stream edge: stable counting sort by tail
    n = num_nodes(g)
    pos = il.offsets[1:n]
    perm = Vector{Int}(undef, num_edges(g))
    for (i, e) in enumerate(edges(g))
        perm[pos[e.u]] = i
        pos[e.u] += 1
    end
    return _dlg(il, edges(g), perm)
end

function Base.show(io::IO, g::DirectedLineGraph)
    print(io, "DirectedLineGraph with ", length(g.edges), " nodes and ", length(g.heads), " edges")
end
