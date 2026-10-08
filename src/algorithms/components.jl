# Further notions of temporal components: source- and sink-based components defined
# by temporal reachability, and snapshot-based components defined by the static
# connectivity of the snapshots (window, persistent and T-interval connected
# components). See the "Temporal Connected Components" family of the Temporal Graph
# Wiki (temporalgraph.notion.site) for the definitions.

function _check_set(g, S)
    X = sort!(unique!([_check_node(g, s) for s in S]))
    isempty(X) && throw(ArgumentError("the set of designated nodes must not be empty"))
    return X
end

# Nodes reached from every node of S (from = true) or reaching every node of S, in the
# subgraph induced by `inside` (all nodes if nothing).
function _common_reach(g::OrderedEdgeList, S::Vector{Int}, ti, from::Bool, inside)
    n = num_nodes(g)
    if inside === nothing
        h, X = g, collect(1:n)
    else
        X = findall(inside)
        h = _induced(g, X, ti)
    end
    id = zeros(Int, n)
    id[X] = eachindex(X)
    any(s -> id[s] == 0, S) && return Int[]
    R = temporal_reachability(h, ti)
    keep = trues(length(X))
    for s in S
        keep .&= from ? R[id[s], :] : R[:, id[s]]
    end
    return X[keep]
end

function _directed_component(g, S, ti, closed::Bool, from::Bool)
    S = _check_set(g, S)
    X = _common_reach(g, S, ti, from, nothing)
    closed || return X
    # greatest fixed point of X ↦ nodes reached from all of S inside X
    while true
        inside = falses(num_nodes(g))
        inside[X] .= true
        Y = _common_reach(g, S, ti, from, inside)
        Y == X && return X
        X = Y
    end
end

"""
    source_component(g::OrderedEdgeList, S, ti = time_interval(g); closed = false)

The *source component* of the nodes `S` (a node or a collection of nodes): the
largest set `X` of nodes such that every node of `S` reaches every node of `X` by a
temporal path inside `ti`. With `closed = true` the paths must stay inside `X`
(which then contains `S`, or is empty if the nodes of `S` cannot reach each other
that way). For a single source this is the set of nodes reachable from it; with
`S` equal to `X` the condition is temporal connectivity.

The union of two such sets is again one, so the component is unique: the open one is
the intersection of the sets reachable from the nodes of `S`, the closed one the
greatest fixed point of restricting this intersection to the subgraph induced by the
current set (each step is a reachability computation).

# References

- S. Bhadra and A. Ferreira. *Complexity of connected components in evolving graphs and the computation of multicast trees in wireless networks.* ADHOC-NOW, 2003. [DOI](https://doi.org/10.1007/978-3-540-39611-6_23)
"""
source_component(g::OrderedEdgeList, S, ti=time_interval(g); closed::Bool=false) =
    _directed_component(g, S isa Integer ? [S] : S, ti, closed, true)

"""
    sink_component(g::OrderedEdgeList, S, ti = time_interval(g); closed = false)

The *sink component* of the nodes `S`: the largest set `X` of nodes such that every
node of `X` reaches every node of `S` by a temporal path inside `ti` (inside `X` with
`closed = true`). See [`source_component`](@ref).
"""
sink_component(g::OrderedEdgeList, S, ti=time_interval(g); closed::Bool=false) =
    _directed_component(g, S isa Integer ? [S] : S, ti, closed, false)

"""
    window_components(g::OrderedEdgeList, ti = time_interval(g))

The *window connected components* of `g` for the time window `ti = (a, b)`: the
connected components of the static undirected graph of the edges with time stamps in
`[a, b]`. Unlike temporal components they ignore the order of the edges. Returns the
components sorted by decreasing size (isolated nodes are components of size 1).
"""
function window_components(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    h = Graphs.SimpleGraph(num_nodes(g))
    for e in edges(g)
        a <= e.t <= b && e.u != e.v && Graphs.add_edge!(h, Int(e.u), Int(e.v))
    end
    return _sorted_components(Graphs.connected_components(h))
end

# Time steps of the snapshots: the grid a, a + Δ, ... of the time interval.
function _time_steps(g::OrderedEdgeList{V,T}, ti, resolution) where {V,T}
    a, b = _interval(T, ti)
    Δ = resolution === nothing ? _time_grid(g, nothing)[3] : convert(T, resolution)
    Δ === nothing && throw(ArgumentError("floating point times need a resolution"))
    Δ > 0 || throw(ArgumentError("the resolution must be positive"))
    return a, b, Δ, fld(b - a, Δ) + 1
end

# Edges (pairs u < v, without loops) of every snapshot, by step index.
function _snapshot_pairs(g::OrderedEdgeList, a, b, Δ, L)
    snaps = [Tuple{Int,Int}[] for _ in 1:L]
    for e in edges(g)
        (a <= e.t <= b && e.u != e.v) || continue
        r = e.t - a
        rem(r, Δ) == 0 || throw(ArgumentError("the time stamps are not on the grid of the resolution"))
        push!(snaps[fld(r, Δ)+1], minmax(Int(e.u), Int(e.v)))
    end
    foreach(x -> unique!(sort!(x)), snaps)
    return snaps
end

# Non-trivial connected components (size ≥ 2) of the graph of the given edges, as
# sorted vectors.
function _components_of(pairs::Vector{Tuple{Int,Int}})
    isempty(pairs) && return Vector{Int}[]
    nodes = unique!(sort!([x for p in pairs for x in p]))
    id(x) = searchsortedfirst(nodes, x)
    h = Graphs.SimpleGraph(length(nodes))
    for (u, v) in pairs
        Graphs.add_edge!(h, id(u), id(v))
    end
    return [nodes[c] for c in Graphs.connected_components(h)]  # sorted, of size ≥ 2
end

"""
    persistent_components(g::OrderedEdgeList, ti = time_interval(g); resolution = nothing,
                          min_length = 1, min_size = 2)

The *persistent connected components* of `g` (Vernet, Pigné and Sanlaville, 2023): pairs `(X, (t₁, t₂))` such that `X` is a connected component of
every snapshot (the static graph of the edges with one time stamp) at the times
`t₁, t₁ + Δ, ..., t₂`, with the interval as long as possible. `Δ` is the
`resolution` of the time steps (by default the gcd of the time stamps minus the start
of `ti`), so a component does not persist across an empty time step.

Returns the components with at least `min_length` steps and `min_size` nodes, as
named tuples `(nodes, interval, length)` sorted by decreasing length and size. One
sweep over the snapshots: `O(L + m α(n))` time for `L` time steps.

# References

- M. Vernet, Y. Pigné, and É. Sanlaville. *A study of connectivity on dynamic graphs: computing persistent connected components.* 4OR 21, 2023. [DOI](https://doi.org/10.1007/s10288-022-00507-3)
"""
function persistent_components(g::OrderedEdgeList{V,T}, ti=time_interval(g); resolution=nothing,
                               min_length::Integer=1, min_size::Integer=2) where {V,T}
    a, b, Δ, L = _time_steps(g, ti, resolution)
    snaps = _snapshot_pairs(g, a, b, Δ, L)
    out = NamedTuple{(:nodes, :interval, :length),Tuple{Vector{Int},Tuple{T,T},Int}}[]
    active = Dict{Vector{Int},Int}()  # component -> first step
    emit(c, first, last) = (last - first + 1 >= min_length && length(c) >= min_size) &&
                           push!(out, (nodes=c, interval=(a + Δ * (first - 1), a + Δ * (last - 1)), length=last - first + 1))
    for k in 1:L
        current = Set(_components_of(snaps[k]))
        for (c, first) in collect(active)
            if !(c in current)
                emit(c, first, k - 1)
                delete!(active, c)
            end
        end
        for c in current
            haskey(active, c) || (active[c] = k)
        end
    end
    for (c, first) in active
        emit(c, first, L)
    end
    return sort!(out; by=x -> (-x.length, -length(x.nodes), x.interval[1], x.nodes))
end

"""
    interval_connected_components(g::OrderedEdgeList, T, ti = time_interval(g); resolution = nothing)

The *`T`-interval connected components* of `g`: the maximal sets of nodes `X` such
that for every window of `T` consecutive time steps the graph of the edges present
in all the snapshots of the window, restricted to `X`, is connected (Kuhn, Lynch and
Oshman, 2010, for the whole graph; the component version is described in the Temporal
Graph Wiki). Time steps are as in [`persistent_components`](@ref). `g` is
`T`-interval connected iff the result is the single set of all nodes.

The maximal sets connected in a family of graphs partition the nodes, and are
computed by refining the partition `{V}` with the connected components of every
window graph until nothing changes. Returns the components sorted by decreasing size.

# References

- F. Kuhn, N. Lynch, and R. Oshman. *Distributed computation in dynamic networks.* STOC, 2010. [DOI](https://doi.org/10.1145/1806689.1806760)
- A. Casteigts, R. Klasing, Y. M. Neggaz, and J. G. Peters. *Efficiently testing T-interval connectivity in dynamic graphs.* CIAC, 2015. [DOI](https://doi.org/10.1007/978-3-319-18173-8_6)
"""
function interval_connected_components(g::OrderedEdgeList{V,TT}, T::Integer, ti=time_interval(g);
                                       resolution=nothing) where {V,TT}
    T >= 1 || throw(ArgumentError("the window length T must be at least 1"))
    a, b, Δ, L = _time_steps(g, ti, resolution)
    snaps = _snapshot_pairs(g, a, b, Δ, L)
    n = num_nodes(g)
    # pairs present in all the snapshots of every window, from the runs of presence
    runs = Dict{Tuple{Int,Int},Vector{Int}}()  # pair -> steps where present
    for k in 1:L, p in snaps[k]
        push!(get!(runs, p, Int[]), k)
    end
    nw = max(L - T + 1, 1)
    windows = [Tuple{Int,Int}[] for _ in 1:nw]
    for (p, ks) in runs
        start = 1
        for i in eachindex(ks)
            if i == length(ks) || ks[i+1] != ks[i] + 1  # run ks[start]..ks[i]
                lo, hi = ks[start], ks[i]
                for w in lo:hi-T+1
                    w <= nw && push!(windows[w], p)
                end
                start = i + 1
            end
        end
    end
    T > L && (windows = [Tuple{Int,Int}[]])  # no complete window: every node alone
    part = zeros(Int, n)  # part id of every node
    nparts = 1
    fill!(part, 1)
    changed = true
    while changed
        changed = false
        for W in windows
            # split every part by the components of the window graph inside it
            h = Graphs.SimpleGraph(n)
            for (u, v) in W
                part[u] == part[v] && Graphs.add_edge!(h, u, v)
            end
            comps = Graphs.connected_components(h)
            length(comps) == nparts && continue
            for (i, c) in enumerate(comps)
                part[c] .= i
            end
            nparts = length(comps)
            changed = true
        end
    end
    groups = [Int[] for _ in 1:nparts]
    for v in 1:n
        push!(groups[part[v]], v)
    end
    return _sorted_components(groups)
end

"""
    delta_temporal_connected_components(g::OrderedEdgeList, Δ, ti = time_interval(g);
                                        resolution = nothing, unilateral = false)

The *Δ-temporal connected components* of `g` (Temporal Graph Wiki): the maximal sets
of nodes in which every node reaches every other node by a temporal path inside every
time window `(τ, τ + Δ)` of `ti = (a, b)`, with `τ` on the grid `a, a + r, ...` of the
time steps (`r` is the `resolution`, as in [`persistent_components`](@ref)) and
`τ + Δ ≤ b`; if `Δ ≥ b - a` the only window is `ti` and these are the
[`temporal_connected_components`](@ref). With `unilateral = true`, of every two nodes
at least one must reach the other in every window. For integer times `Δ` is rounded
down.

A window is only needed if it starts at `a` or at the first grid point after a time
stamp: any other window contains all the edges of the last such window before it,
hence all its paths. The components are the maximal cliques of the graph joining the
nodes that reach each other in every window, enumerated with the Bron–Kerbosch
algorithm as for the temporal connected components. `O(K m n / 64)` time for the
reachability in the `K ≤ m + 1` windows, plus the clique enumeration (exponential in
the worst case). Returns the components sorted by decreasing size.
"""
function delta_temporal_connected_components(g::OrderedEdgeList{V,T}, Δ::Real, ti=time_interval(g);
                                             resolution=nothing, unilateral::Bool=false) where {V,T}
    Δ >= 0 || throw(ArgumentError("Δ must be non-negative"))
    a, b, r, _ = _time_steps(g, ti, resolution)
    D = T <: Integer ? floor(T, Δ) : convert(T, Δ)
    last = b - D < a ? a : a + r * fld(b - D - a, r)  # start of the last window
    starts = T[a]
    for e in edges(g)
        a <= e.t < last && push!(starts, a + r * (fld(e.t - a, r) + 1))  # first grid point after e.t
    end
    n = num_nodes(g)
    R = trues(n, n)  # pairs connected in all the windows so far
    for τ in unique!(sort!(starts))
        Rw = temporal_reachability(g, (τ, min(τ + D, b)))
        # unilateral: the direction may change from one window to the next
        unilateral ? (R .&= Rw .| permutedims(Rw)) : (R .&= Rw)
        count(R) == n && break  # only the diagonal is left
    end
    return _sorted_components(Graphs.maximal_cliques(_component_graph(R, unilateral)))
end

"""
    stream_components(g::OrderedEdgeList, ti = time_interval(g); gap = 0)

The connected components of `g` seen as a *stream graph* (Latapy, Viard and Magnien,
2018), whose elements are temporal nodes instead of nodes. Every edge
`(u, v, t, tt)` inside `ti` is present during `[t, t + tt]`, and so are its endpoints:
a node is present while it has an edge, and its presence intervals at most `gap`
apart are merged (it stays present in between). The temporal nodes are the pairs
`(v, (s, f))` of a node and a maximal presence interval, and an edge joins the
temporal nodes of its endpoints that contain it. The order of time and the direction
of the edges are ignored: a component is a set of temporal nodes linked by edges
through time, e.g. a group of people that keeps meeting without long interruptions.

Returns the components (vectors of temporal nodes, sorted) sorted by decreasing
number of temporal nodes; nodes without edges inside `ti` belong to none.
`O(m log m)` time.

# References

- M. Latapy, T. Viard, and C. Magnien. *Stream graphs and link streams for the modeling of interactions over time.* Social Network Analysis and Mining 8, 2018. [DOI](https://doi.org/10.1007/s13278-018-0537-7), [arXiv](https://arxiv.org/abs/1710.04073)
"""
function stream_components(g::OrderedEdgeList{V,T}, ti=time_interval(g); gap::Real=0) where {V,T}
    gap >= 0 || throw(ArgumentError("gap must be non-negative"))
    a, b = _interval(T, ti)
    es = edges(g)
    ids = _window_edges(g, a, b)
    at = [Int[] for _ in 1:num_nodes(g)]  # positions in ids of the edges of every node
    for (x, i) in enumerate(ids)
        push!(at[es[i].u], x)
        es[i].v == es[i].u || push!(at[es[i].v], x)
    end
    tnodes = Tuple{Int,Tuple{T,T}}[]
    tu = zeros(Int, length(ids))  # temporal node of the tail and of the head of every edge
    tv = zeros(Int, length(ids))
    for v in eachindex(at)
        xs = sort!(at[v]; by=x -> es[ids[x]].t)
        k = 1
        while k <= length(xs)
            s, f = es[ids[xs[k]]].t, es[ids[xs[k]]].t + es[ids[xs[k]]].tt
            j = k
            while j < length(xs) && es[ids[xs[j+1]]].t <= f + gap
                j += 1
                f = max(f, es[ids[xs[j]]].t + es[ids[xs[j]]].tt)
            end
            push!(tnodes, (v, (s, f)))
            for x in xs[k:j]
                es[ids[x]].u == v && (tu[x] = length(tnodes))
                es[ids[x]].v == v && (tv[x] = length(tnodes))
            end
            k = j + 1
        end
    end
    h = Graphs.SimpleGraph(length(tnodes))
    for x in eachindex(ids)
        tu[x] == tv[x] || Graphs.add_edge!(h, tu[x], tv[x])
    end
    comps = [sort!(tnodes[c]) for c in Graphs.connected_components(h)]
    return sort!(comps; by=c -> (-length(c), c))
end
