"""
    OrderedEdgeList(n, edges, ti = span(edges); original_ids = 1:n)

Temporal graph with nodes `1:n` stored as a chronologically ordered list of temporal
edges (the "edge stream" representation). This is the main representation: most
algorithms run a single scan over the contiguous edge vector.

`edges` is a vector of [`TemporalEdge`](@ref)s or of tuples `(u, v, t, tt)`. It is
copied, and sorted by `(t, tt, u, v)` if it is not sorted by time. `ti` is the time
interval spanned by the graph and defaults to `(minimum t, maximum t + tt)`.
`original_ids[i]` is the id that node `i` had in the input file (see
[`load_ordered_edge_list`](@ref)).
"""
struct OrderedEdgeList{V<:Integer,T<:Real}
    num_nodes::Int
    edges::Vector{TemporalEdge{V,T}}
    ti::TimeInterval{T}
    original_ids::Vector{Int}
    has_zero_tt::Bool   # some edge has transition time 0
    max_arrival::T      # largest t + tt of an edge

    function OrderedEdgeList{V,T}(n::Integer, es::AbstractVector, ti=nothing;
                                  original_ids::AbstractVector{<:Integer}=1:n) where {V,T}
        n >= 0 || throw(ArgumentError("number of nodes must be non-negative"))
        length(original_ids) == n || throw(ArgumentError("original_ids must have length $n"))
        v = _chronological(TemporalEdge{V,T}, es)
        _check_edges(n, v)
        span = _span(v)
        interval = ti === nothing ? span : _interval(T, ti)
        return new{V,T}(Int(n), v, interval, collect(Int, original_ids), any(e -> iszero(e.tt), v), span[2])
    end
end

OrderedEdgeList(n::Integer, es::AbstractVector{TemporalEdge{V,T}}, ti=nothing; kwargs...) where {V,T} =
    OrderedEdgeList{V,T}(n, es, ti; kwargs...)

function OrderedEdgeList(n::Integer, es::AbstractVector{<:Tuple}, ti=nothing; kwargs...)
    return OrderedEdgeList(n, [TemporalEdge(e...) for e in es], ti; kwargs...)
end

num_nodes(g::OrderedEdgeList) = g.num_nodes
num_edges(g::OrderedEdgeList) = length(g.edges)
"""
    edges(g)

The temporal edges of `g` (this extends `Graphs.edges`).
"""
edges(g::OrderedEdgeList) = g.edges
time_interval(g::OrderedEdgeList) = g.ti
time_type(::OrderedEdgeList{V,T}) where {V,T} = T
node_type(::OrderedEdgeList{V,T}) where {V,T} = V

"""
    original_id(g, u)

The id node `u` had in the input file.
"""
original_id(g::OrderedEdgeList, u::Integer) = g.original_ids[u]

"""
    original_ids(g)

Vector mapping each node to the id it had in the input file.
"""
original_ids(g::OrderedEdgeList) = g.original_ids

"""
    node_map(g)

Dictionary mapping the ids of the input file to the node ids of `g`.
"""
node_map(g::OrderedEdgeList) = Dict(id => i for (i, id) in enumerate(g.original_ids))

function Base.:(==)(a::OrderedEdgeList, b::OrderedEdgeList)
    return a.num_nodes == b.num_nodes && a.edges == b.edges && a.ti == b.ti
end

function Base.show(io::IO, g::OrderedEdgeList{V,T}) where {V,T}
    print(io, "OrderedEdgeList{", V, ", ", T, "} with ", g.num_nodes, " nodes, ", length(g.edges),
          " temporal edges, time interval ", g.ti)
end

