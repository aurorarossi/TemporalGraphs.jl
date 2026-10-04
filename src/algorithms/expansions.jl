# Static expansions: directed static graphs whose paths are the temporal walks of a
# temporal graph, the standard tool to reduce temporal questions (reachability,
# optimal paths, flows, cuts) to static ones.
#
# The vertex expansion (time-expanded graph) has one node (v, τ) per node and time,
# waiting arcs between consecutive times of the same node and one arc
# (u, t) → (v, t + tt) per temporal edge. The edge expansion (temporal line graph)
# has one node per temporal edge and an arc e → f whenever f can follow e.

"""
    vertex_expansion(g::OrderedEdgeList, ti = time_interval(g); compressed = true)

The *vertex expansion* (time-expanded graph) of `g` inside `ti`, as a
`Graphs.SimpleDiGraph`: a node for every pair `(v, τ)` of a node and a time, an arc
`(v, τ) → (v, τ')` between consecutive times of `v` (waiting) and an arc
`(u, t) → (v, t + tt)` for every temporal edge (traversal). Temporal walks of `g` are
exactly the directed walks of the expansion, so `v` is reachable from `u` iff some
copy of `v` is reachable from the first copy of `u`.

With `compressed = true` a node has copies only at the times at which an edge leaves
or reaches it (`O(n + m)` nodes and arcs, as used for flows); with
`compressed = false` every node has a copy at every such time of the whole graph
(`n L` nodes for `L` distinct times, the classical layered expansion).

Returns `(graph, nodes)` where `nodes[i] = (v, τ)` is the pair of expansion node `i`.

# References

- H. Wu, J. Cheng, Y. Ke, S. Huang, Y. Huang, and H. Wu. *Efficient algorithms for temporal path computation.* IEEE Transactions on Knowledge and Data Engineering 28(11), 2016. [DOI](https://doi.org/10.1109/TKDE.2016.2594065)
- E. C. Akrida, J. Czyzowicz, L. Gąsieniec, Ł. Kuszner, and P. G. Spirakis. *Temporal flows in temporal networks.* Journal of Computer and System Sciences 103, 2019. [DOI](https://doi.org/10.1016/j.jcss.2019.02.003), [arXiv](https://arxiv.org/abs/1606.01091)
"""
function vertex_expansion(g::OrderedEdgeList{V,T}, ti=time_interval(g); compressed::Bool=true) where {V,T}
    a, b = _interval(T, ti)
    es = edges(g)
    n = num_nodes(g)
    ids = _window_edges(g, a, b)
    times = [T[] for _ in 1:n]
    if compressed
        for i in ids
            push!(times[es[i].u], es[i].t)
            push!(times[es[i].v], es[i].t + es[i].tt)
        end
        foreach(x -> unique!(sort!(x)), times)
    else
        layers = unique!(sort!(vcat([es[i].t for i in ids], [es[i].t + es[i].tt for i in ids])))
        foreach(v -> append!(times[v], layers), 1:n)
    end
    off = ones(Int, n + 1)
    for v in 1:n
        off[v+1] = off[v] + length(times[v])
    end
    nodes = [(v, τ) for v in 1:n for τ in times[v]]
    h = Graphs.SimpleDiGraph(length(nodes))
    for v in 1:n, k in off[v]:off[v+1]-2
        Graphs.add_edge!(h, k, k + 1)
    end
    node(v, τ) = off[v] + searchsortedfirst(times[v], τ) - 1
    for i in ids
        e = es[i]
        Graphs.add_edge!(h, node(e.u, e.t), node(e.v, e.t + e.tt))
    end
    return (graph=h, nodes=nodes)
end

"""
    edge_expansion(g::OrderedEdgeList, ti = time_interval(g); β = Inf, hubs = false)

The *edge expansion* (temporal line graph) of `g` inside `ti`, as a
`Graphs.SimpleDiGraph`: a node for every temporal edge and an arc `e → f` whenever
`f` can follow `e`, i.e., `f` leaves the head of `e` at time `arr(e) = e.t + e.tt` or
later (and at most `arr(e) + β` with a waiting constraint). Directed walks of the
expansion are exactly the temporal walks (`β`-restless walks) of `g`.

The expansion can have `Θ(m²)` arcs. With `hubs = true` (only without waiting
constraint) one extra *hub* node per departure time `(v, τ)` is used instead: hubs of
a node are chained by increasing time, an edge points to the first hub of its head at
or after its arrival and a hub to the edges leaving at its time, which gives `O(m)`
arcs with the same reachability between edge nodes.

Returns `(graph, edges, hubs)`: node `i ≤ length(edges)` is the temporal edge
`edges(g)[edges[i]]`, the following nodes are the hubs `hubs[j] = (v, τ)`.
"""
function edge_expansion(g::OrderedEdgeList{V,T}, ti=time_interval(g); β::Real=Inf, hubs::Bool=false) where {V,T}
    hubs && isfinite(β) && throw(ArgumentError("hubs cannot be combined with a waiting constraint"))
    isfinite(β) && β < 0 && throw(ArgumentError("β must be non-negative"))
    a, b = _interval(T, ti)
    es = edges(g)
    n = num_nodes(g)
    ids = _window_edges(g, a, b)
    k = length(ids)
    out = [Int[] for _ in 1:n]  # positions in ids of the edges leaving every node, by departure
    for (x, i) in enumerate(ids)
        push!(out[es[i].u], x)
    end
    # first position in o of an edge leaving at or after τ (after τ if strict)
    function first_from(o, τ, strict)
        lo, hi = 1, length(o)
        while lo <= hi
            mid = (lo + hi) >>> 1
            t = es[ids[o[mid]]].t
            (strict ? t <= τ : t < τ) ? (lo = mid + 1) : (hi = mid - 1)
        end
        return lo
    end
    hublist = Tuple{Int,T}[]
    if !hubs
        h = Graphs.SimpleDiGraph(k)
        for (x, i) in enumerate(ids)
            e = es[i]
            o = out[e.v]
            arr = e.t + e.tt
            last = isfinite(β) ? first_from(o, arr + convert(T, β), true) - 1 : length(o)
            for y in first_from(o, arr, false):last
                Graphs.add_edge!(h, x, o[y])
            end
        end
        return (graph=h, edges=ids, hubs=hublist)
    end
    # hubs of every node: its distinct departure times, numbered after the edge nodes
    hubtimes = [unique!([es[ids[y]].t for y in out[v]]) for v in 1:n]
    hubstart = ones(Int, n + 1) .* (k + 1)
    for v in 1:n
        hubstart[v+1] = hubstart[v] + length(hubtimes[v])
        append!(hublist, [(v, τ) for τ in hubtimes[v]])
    end
    hub(v, τ) = hubstart[v] + searchsortedfirst(hubtimes[v], τ) - 1
    h = Graphs.SimpleDiGraph(k + length(hublist))
    for v in 1:n
        for j in hubstart[v]:hubstart[v+1]-2  # waiting
            Graphs.add_edge!(h, j, j + 1)
        end
        for y in out[v]
            Graphs.add_edge!(h, hub(v, es[ids[y]].t), y)
        end
    end
    for (x, i) in enumerate(ids)
        e = es[i]
        ts = hubtimes[e.v]
        j = searchsortedfirst(ts, e.t + e.tt)  # first departure from the head at or after the arrival
        j <= length(ts) && Graphs.add_edge!(h, x, hubstart[e.v] + j - 1)
    end
    return (graph=h, edges=ids, hubs=hublist)
end
