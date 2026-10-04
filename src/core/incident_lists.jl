"""
    IncidentLists(n, edges, ti = span(edges))

Temporal graph in incident lists representation: for every node the outgoing
temporal edges in chronological order. The lists are stored contiguously in one
vector (compressed sparse row layout); `out_edges(g, u)` returns a view.
"""
struct IncidentLists{V<:Integer,T<:Real}
    num_nodes::Int
    offsets::Vector{Int}
    edges::Vector{TemporalEdge{V,T}}
    ti::TimeInterval{T}
end

function IncidentLists(n::Integer, es::AbstractVector{TemporalEdge{V,T}}, ti=nothing) where {V,T}
    sorted = issorted(es; by=e -> e.t) ? es : sort(es)
    _check_edges(n, sorted)
    offsets = zeros(Int, n + 1)
    for e in sorted
        offsets[e.u+1] += 1
    end
    offsets[1] = 1
    for i in 2:n+1
        offsets[i] += offsets[i-1]
    end
    pos = offsets[1:n]
    out = similar(sorted)
    @inbounds for e in sorted  # stable: keeps the chronological order inside each list
        out[pos[e.u]] = e
        pos[e.u] += 1
    end
    return IncidentLists{V,T}(Int(n), offsets, out, ti === nothing ? _span(sorted) : _interval(T, ti))
end

num_nodes(g::IncidentLists) = g.num_nodes
num_edges(g::IncidentLists) = length(g.edges)
edges(g::IncidentLists) = g.edges
time_interval(g::IncidentLists) = g.ti
time_type(::IncidentLists{V,T}) where {V,T} = T
node_type(::IncidentLists{V,T}) where {V,T} = V

"""
    out_edges(g::IncidentLists, u)

The outgoing temporal edges of `u` in chronological order (a view).
"""
@inline out_edges(g::IncidentLists, u::Integer) = @inbounds view(g.edges, g.offsets[u]:g.offsets[u+1]-1)

"""
    out_degree(g::IncidentLists, u)

Number of outgoing temporal edges of `u`.
"""
@inline out_degree(g::IncidentLists, u::Integer) = @inbounds g.offsets[u+1] - g.offsets[u]

# first position in the list of `u` of an edge with t ≥ tmin (binary search)
@inline function _first_edge_from(g::IncidentLists, u, tmin)
    lo = @inbounds g.offsets[u]
    hi = @inbounds g.offsets[u+1] - 1
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if g.edges[mid].t < tmin
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return lo
end

function Base.show(io::IO, g::IncidentLists{V,T}) where {V,T}
    print(io, "IncidentLists{", V, ", ", T, "} with ", g.num_nodes, " nodes, ", length(g.edges),
          " temporal edges, time interval ", g.ti)
end
