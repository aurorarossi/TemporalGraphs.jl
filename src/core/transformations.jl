# Graphs.jl style accessors
nv(g::Union{OrderedEdgeList,IncidentLists,TRSGraph}) = num_nodes(g)
ne(g::Union{OrderedEdgeList,IncidentLists,TRSGraph}) = num_edges(g)

"""
    to_incident_lists(g::OrderedEdgeList)

Transform `g` into its [`IncidentLists`](@ref) representation.
"""
to_incident_lists(g::OrderedEdgeList) = IncidentLists(num_nodes(g), edges(g), time_interval(g))
to_incident_lists(g::IncidentLists) = g

"""
    to_ordered_edge_list(g::IncidentLists)

Transform `g` into its [`OrderedEdgeList`](@ref) representation.
"""
to_ordered_edge_list(g::IncidentLists) = OrderedEdgeList(num_nodes(g), sort(g.edges), time_interval(g))
to_ordered_edge_list(g::OrderedEdgeList) = g

"""
    to_aggregated_edge_list(g::OrderedEdgeList)

The static graph underlying `g` in which every edge `(u, v)` is weighted with the
number of temporal edges from `u` to `v`. The result is sorted by `(u, v)`.
"""
function to_aggregated_edge_list(g::OrderedEdgeList)
    pairs = sort!([(e.u, e.v) for e in edges(g)])
    result = StaticWeightedEdge[]
    i = 1
    while i <= length(pairs)
        j = i
        while j < length(pairs) && pairs[j+1] == pairs[i]
            j += 1
        end
        push!(result, StaticWeightedEdge(pairs[i][1], pairs[i][2], j - i + 1))
        i = j + 1
    end
    return result
end

_rebuild(g::OrderedEdgeList, es::Vector{<:TemporalEdge}) =
    OrderedEdgeList(num_nodes(g), es; original_ids=g.original_ids)

"""
    normalize_graph(g::OrderedEdgeList, remove_loops::Bool = true)

Remove multiple (identical) edges and, if `remove_loops`, self-loops. The resulting
edges are sorted by `(t, tt, u, v)` and the time interval is recomputed.
"""
function normalize_graph(g::OrderedEdgeList, remove_loops::Bool=true)
    es = remove_loops ? filter(e -> e.u != e.v, edges(g)) : copy(edges(g))
    return _rebuild(g, unique!(sort!(es)))
end

"""
    scale_timestamps(g::OrderedEdgeList, factor)

Multiply all time stamps by `factor` (for integer times rounding half away from
zero) and remove the resulting multiple edges.
"""
function scale_timestamps(g::OrderedEdgeList{V,T}, factor::Real) where {V,T}
    es = [TemporalEdge{V,T}(e.u, e.v, _scale(e.t, factor), e.tt) for e in edges(g)]
    return _rebuild(g, unique!(sort!(es)))
end

_scale(t::T, f) where {T<:Integer} = round(T, t * f, RoundNearestTiesAway)
_scale(t::T, f) where {T<:AbstractFloat} = T(t * f)

"""
    unit_transition_times(g::OrderedEdgeList, val = 1)

Replace all transition times by `val`.
"""
function unit_transition_times(g::OrderedEdgeList{V,T}, val::Real=one(T)) where {V,T}
    es = [TemporalEdge{V,T}(e.u, e.v, e.t, val) for e in edges(g)]
    return _rebuild(g, sort!(es))
end

"""
    make_undirected(g::OrderedEdgeList)

Insert for each edge `(u, v, t, tt)` the reverse edge `(v, u, t, tt)` and remove
multiple edges.
"""
function make_undirected(g::OrderedEdgeList{V,T}) where {V,T}
    es = Vector{TemporalEdge{V,T}}(undef, 2 * num_edges(g))
    for (i, e) in enumerate(edges(g))
        es[2i-1] = e
        es[2i] = TemporalEdge{V,T}(e.v, e.u, e.t, e.tt)
    end
    return _rebuild(g, unique!(sort!(es)))
end

function _reversed_edges(es::Vector{TemporalEdge{V,T}}, a::T, b::T) where {V,T}
    return sort!([TemporalEdge{V,T}(e.v, e.u, a + b - e.t - e.tt, e.tt) for e in es if e.t >= a && e.t + e.tt <= b])
end

"""
    reverse(g, ti = time_interval(g))

Temporal transpose (reverse temporal graph) of `g` restricted to `ti = (a, b)`, see
Oettershagen and Mutzel, "Efficient top-k temporal closeness calculation in temporal
networks", ICDM 2020. Each edge `(u, v, t, tt)` inside `ti` becomes
`(v, u, a + b - t - tt, tt)`, so that a path from `u` to `w` with duration `d` becomes a
path from `w` to `u` with duration `d` and the interval is preserved.

# References

- L. Oettershagen and P. Mutzel. *Efficient top-k temporal closeness calculation in temporal networks.* ICDM, 2020. [DOI](https://doi.org/10.1109/ICDM50108.2020.00049)
"""
function Base.reverse(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    return OrderedEdgeList(num_nodes(g), _reversed_edges(edges(g), a, b), (a, b); original_ids=g.original_ids)
end

function Base.reverse(g::IncidentLists{V,T}, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    return IncidentLists(num_nodes(g), _reversed_edges(edges(g), a, b), (a, b))
end
