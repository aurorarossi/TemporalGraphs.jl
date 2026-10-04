_endpoints(e) = (Int(e.u), Int(e.v))
_endpoints(e::Tuple) = (Int(e[1]), Int(e[2]))

function _max_node(es)
    m = 0
    for e in es
        m = max(m, _endpoints(e)...)
    end
    return m
end

"""
    kcores(edges, n = maximum node id)

Core numbers of the nodes `1:n` of the undirected simple graph given by `edges`
(anything with fields `u` and `v`, or tuples `(u, v, ...)`). Multiple edges, edge
directions and self-loops are ignored. Computed with `Graphs.core_number` (O(m)).
"""
function kcores(es::AbstractVector, n::Integer=_max_node(es))
    sg = Graphs.SimpleGraph(n)
    for e in es
        u, v = _endpoints(e)
        u != v && Graphs.add_edge!(sg, u, v)
    end
    return Graphs.core_number(sg)
end

"""
    static_graph(g::OrderedEdgeList, ti = nothing; directed = true)

The static graph underlying `g` as a `Graphs.SimpleDiGraph` (or `SimpleGraph` if
`directed = false`), with an edge `(u, v)` iff there is a temporal edge from `u` to
`v`. With a time window `ti = (a, b)`, only the edges with time stamps in `[a, b]`
are used (the *window graph*; see also [`window_components`](@ref) and
[`snapshots`](@ref)). Use [`to_aggregated_edge_list`](@ref) for the contact counts.
"""
function static_graph(g::OrderedEdgeList{V,T}, ti=nothing; directed::Bool=true) where {V,T}
    a, b = ti === nothing ? (typemin(T), typemax(T)) : _interval(T, ti)
    sg = directed ? Graphs.SimpleDiGraph(num_nodes(g)) : Graphs.SimpleGraph(num_nodes(g))
    for e in edges(g)
        a <= e.t <= b && Graphs.add_edge!(sg, Int(e.u), Int(e.v))
    end
    return sg
end

"""
    temporal_khcores(g::OrderedEdgeList, h)

Core numbers of the `(k, h)`-cores (Wu et al., 2015): the core numbers of the
static graph that connects two nodes iff they have at least `h` temporal contacts.

# References

- H. Wu, J. Cheng, Y. Lu, Y. Ke, Y. Huang, D. Yan, and H. Wu. *Core decomposition in large temporal graphs.* IEEE International Conference on Big Data, 2015. [DOI](https://doi.org/10.1109/BigData.2015.7363809)
"""
function temporal_khcores(g::OrderedEdgeList, h::Integer)
    agg = [e for e in to_aggregated_edge_list(g) if e.weight >= h]
    return kcores(agg, num_nodes(g))
end

"""
    temporal_lkcores(g::OrderedEdgeList, L, k)

Nodes of the largest `(L, k)`-lasting core (Hung and Tseng, 2021): a set of nodes that
forms a `k`-core during `L` consecutive time stamps. For every window of `L`
consecutive distinct time stamps the static graph of the node pairs that have a
contact at each of the `L` time stamps is built, and the nodes with core number at
least `k` of the window with the most such nodes are returned (sorted). The contact
counts are updated incrementally while the window slides.

# References

- W.-C. Hung and C.-Y. Tseng. *Maximum (L, K)-lasting cores in temporal social networks.* DASFAA, 2021. [DOI](https://doi.org/10.1007/978-3-030-73216-5_23)
"""
function temporal_lkcores(g::OrderedEdgeList{V,T}, L::Integer, k::Integer) where {V,T}
    L >= 1 || throw(ArgumentError("L must be positive"))
    contacts = unique!([(e.t, e.u, e.v) for e in edges(g)])  # edges are sorted by time
    times = unique!([c[1] for c in contacts])
    best = Int[]
    length(times) < L && return best
    count = Dict{Tuple{V,V},Int}()
    lasting = Tuple{V,V}[]
    lo = 1   # first contact inside the window
    hi = 1   # first contact after the window
    for i in 1:length(times)-L+1
        while hi <= length(contacts) && contacts[hi][1] <= times[i+L-1]
            key = (contacts[hi][2], contacts[hi][3])
            count[key] = get(count, key, 0) + 1
            hi += 1
        end
        while contacts[lo][1] < times[i]
            key = (contacts[lo][2], contacts[lo][3])
            c = count[key] - 1
            c == 0 ? delete!(count, key) : (count[key] = c)
            lo += 1
        end
        empty!(lasting)
        for (key, c) in count
            c == L && push!(lasting, key)
        end
        cores = kcores(lasting, num_nodes(g))
        nodes = findall(>=(k), cores)
        length(nodes) > length(best) && (best = nodes)
    end
    return best
end
