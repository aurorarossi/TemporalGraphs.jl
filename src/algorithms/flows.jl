# Temporal flows, cuts, disjoint paths and separators.
#
# Flows (Akrida, Czyzowicz, Gąsieniec, Kuszner and Spirakis, "Temporal flows in
# temporal networks", JCSS 2019) are computed as static maximum flows in the
# time-expanded network: one node (v, τ) for every time τ at which an edge leaves or
# reaches v, holdover arcs (v, τ) → (v, τ') between consecutive times (capacity the
# buffer of v, unbounded by default) and one arc (u, t) → (v, t + tt) per temporal
# edge with its capacity. An arrival and a departure at the same time share a node,
# so an edge can be followed by an edge leaving at its arrival time. With unbounded
# buffers the maximum flow equals the minimum capacity of a set of temporal edges
# whose removal disconnects s from z (a temporal cut), and with unit capacities this
# is Menger's theorem for edge-disjoint temporal paths (Berman, 1996). With a
# waiting constraint β the flow follows restless walks, and the network has one
# node per temporal edge instead.

# Edges inside ti.
_window_edges(g::OrderedEdgeList{V,T}, a, b) where {V,T} = [i for (i, e) in enumerate(edges(g)) if a <= e.t && e.t + e.tt <= b]

_flow_inf(::Type{C}) where {C<:Integer} = typemax(C) ÷ 4
_flow_inf(::Type{C}) where {C<:AbstractFloat} = C(Inf)

# Time-expanded network. Returns the network, the super source and sink, the arc of
# every edge of `ids` (0 if it has none) and, with `departures`, uses one node per
# departure (u, t) between (u, t) and the edges, with capacity 1.
function _time_expanded(g::OrderedEdgeList{V,T}, ids::Vector{Int}, s::Int, z::Int, cap::AbstractVector{C};
                        buffers=nothing, departures::Bool=false) where {V,T,C}
    es = edges(g)
    n = num_nodes(g)
    times = [T[] for _ in 1:n]
    for i in ids
        e = es[i]
        push!(times[e.u], e.t)
        push!(times[e.v], e.t + e.tt)
    end
    off = ones(Int, n + 1)
    for v in 1:n
        unique!(sort!(times[v]))
        off[v+1] = off[v] + length(times[v])
    end
    node(v, τ) = off[v] + searchsortedfirst(times[v], τ) - 1
    N = off[n+1] - 1
    S, Z = N + 1, N + 2
    D = Dict{Tuple{Int,T},Int}()  # departure nodes
    if departures
        for i in ids
            e = es[i]
            get!(D, (Int(e.u), e.t), N + 2 + length(D) + 1)
        end
    end
    net = _FlowNetwork{C}(N + 2 + length(D))
    inf = _flow_inf(C)
    for v in 1:n, k in off[v]:off[v+1]-2
        c = buffers === nothing || v == s || v == z ? inf : convert(C, buffers[v])
        _add_arc!(net, k, k + 1, c)
    end
    for ((u, t), d) in D
        _add_arc!(net, node(u, t), d, one(C))
    end
    arc = zeros(Int, length(es))
    for i in ids
        e = es[i]
        from = departures ? D[(Int(e.u), e.t)] : node(e.u, e.t)
        arc[i] = _add_arc!(net, from, node(e.v, e.t + e.tt), departures ? inf : cap[i])
    end
    isempty(times[s]) || _add_arc!(net, S, off[s], inf)
    isempty(times[z]) || _add_arc!(net, off[z+1] - 1, Z, inf)
    return _build!(net), S, Z, arc
end

# Network of restless walks: nodes e_in = 2i - 1 and e_out = 2i for every edge i,
# e_in → e_out with the capacity of the edge, e_out → f_in when f can follow e
# (arr(e) ≤ f.t ≤ arr(e) + β), the source to the edges leaving s and the edges
# entering z to the sink.
function _restless_network(g::OrderedEdgeList{V,T}, ids::Vector{Int}, s::Int, z::Int, cap::AbstractVector{C}, β) where {V,T,C}
    es = edges(g)
    n = num_nodes(g)
    m = length(es)
    out = [Int[] for _ in 1:n]  # edges of ids leaving every node, by departure
    for i in ids
        push!(out[es[i].u], i)
    end
    net = _FlowNetwork{C}(2m + 2)
    S, Z = 2m + 1, 2m + 2
    inf = _flow_inf(C)
    arc = zeros(Int, m)
    βT = convert(T, β)
    for i in ids
        e = es[i]
        arc[i] = _add_arc!(net, 2i - 1, 2i, cap[i])
        Int(e.u) == s && _add_arc!(net, S, 2i - 1, inf)
        Int(e.v) == z && _add_arc!(net, 2i, Z, inf)
        a = e.t + e.tt
        o = out[e.v]
        lo, hi = 1, length(o)  # first edge of o leaving at or after a
        while lo <= hi
            mid = (lo + hi) >>> 1
            es[o[mid]].t < a ? (lo = mid + 1) : (hi = mid - 1)
        end
        for x in lo:length(o)
            f = o[x]
            es[f].t > a + βT && break
            _add_arc!(net, 2i, 2f - 1, inf)
        end
    end
    return _build!(net), S, Z, arc
end

function _flow_network(g::OrderedEdgeList{V,T}, s, z, ti, capacity, buffers, β) where {V,T}
    s, z = _check_node(g, s), _check_node(g, z)
    s == z && throw(ArgumentError("source and sink must be different"))
    a, b = _interval(T, ti)
    ids = _window_edges(g, a, b)
    cap = capacity === nothing ? ones(Int, num_edges(g)) : capacity
    length(cap) == num_edges(g) || throw(ArgumentError("capacity must have one entry per edge"))
    all(>=(0), cap) || throw(ArgumentError("capacities must be non-negative"))
    if isfinite(β)
        buffers === nothing || throw(ArgumentError("buffers and waiting constraints cannot be combined"))
        β >= 0 || throw(ArgumentError("β must be non-negative"))
        return (_restless_network(g, ids, s, z, cap, β)..., ids)
    end
    buffers === nothing || length(buffers) == num_nodes(g) || throw(ArgumentError("buffers must have one entry per node"))
    return (_time_expanded(g, ids, s, z, cap; buffers=buffers)..., ids)
end

"""
    temporal_max_flow(g::OrderedEdgeList, s, z, ti = time_interval(g);
                      capacity = nothing, buffers = nothing, β = Inf)

Maximum temporal flow from `s` to `z` inside `ti` (Akrida, Czyzowicz, Gąsieniec,
Kuszner and Spirakis, 2019): the largest amount that can be sent from `s` to `z`
when temporal edge `i` carries at most `capacity[i]` (1 by default) at its time and
nodes can store any amount between edges (at most `buffers[v]` at node `v` if
given; `s` and `z` always have unbounded storage). With a waiting constraint `β`,
the flow follows `β`-restless walks instead (nothing can be stored for longer than
`β`).

Returns `(value, flow)` where `flow[i]` is the amount on edge `i` of `edges(g)`. With
the default unit capacities the value is the maximum number of edge-disjoint
temporal paths (see [`temporal_edge_disjoint_paths`](@ref)). Computed as a static
maximum flow (Dinic's algorithm) in the time-expanded network, of size `O(n + m)`
(`O(Σ_v deg⁻(v) deg⁺(v))` with `β`).

# References

- E. C. Akrida, J. Czyzowicz, L. Gąsieniec, Ł. Kuszner, and P. G. Spirakis. *Temporal flows in temporal networks.* Journal of Computer and System Sciences 103, 2019. [DOI](https://doi.org/10.1016/j.jcss.2019.02.003), [arXiv](https://arxiv.org/abs/1606.01091)
"""
function temporal_max_flow(g::OrderedEdgeList, s::Integer, z::Integer, ti=time_interval(g);
                           capacity::Union{Nothing,AbstractVector{<:Real}}=nothing,
                           buffers::Union{Nothing,AbstractVector{<:Real}}=nothing, β::Real=Inf)
    net, S, Z, arc, ids = _flow_network(g, s, z, ti, capacity, buffers, β)
    C = eltype(net.cap)
    value = _max_flow!(net, S, Z)
    flow = zeros(C, num_edges(g))
    for i in ids
        flow[i] = net.cap[_reverse_arc(arc[i])]
    end
    return (value=value, flow=flow)
end

"""
    temporal_min_cut(g::OrderedEdgeList, s, z, ti = time_interval(g); capacity = nothing, β = Inf)

A minimum temporal cut between `s` and `z`: a set of temporal edges of minimum total
capacity (number of edges by default) whose removal leaves no temporal path from `s`
to `z` inside `ti` (no `β`-restless walk with a waiting constraint). Its capacity
equals the [`temporal_max_flow`](@ref) (Akrida et al., 2019, Theorem 1; Berman, 1996).
Returns `(value, edges)` with the indices of the cut edges in `edges(g)`.

# References

- E. C. Akrida, J. Czyzowicz, L. Gąsieniec, Ł. Kuszner, and P. G. Spirakis. *Temporal flows in temporal networks.* Journal of Computer and System Sciences 103, 2019. [DOI](https://doi.org/10.1016/j.jcss.2019.02.003), [arXiv](https://arxiv.org/abs/1606.01091)
- K. A. Berman. *Vulnerability of scheduled networks and a generalization of Menger's theorem.* Networks 28(3), 1996. [DOI](https://doi.org/10.1002/%28SICI%291097-0037%28199610%2928%3A3%3C125%3A%3AAID-NET1%3E3.0.CO%3B2-P)
"""
function temporal_min_cut(g::OrderedEdgeList, s::Integer, z::Integer, ti=time_interval(g);
                          capacity::Union{Nothing,AbstractVector{<:Real}}=nothing, β::Real=Inf)
    net, S, Z, arc, ids = _flow_network(g, s, z, ti, capacity, nothing, β)
    value = _max_flow!(net, S, Z)
    side = _residual_reachable(net, S)
    cut = [i for i in ids if side[net.tail[arc[i]]] && !side[net.head[arc[i]]]]
    return (value=value, edges=cut)
end

# Decompose the flow of a unit-capacity network into paths of arcs from S to Z.
function _flow_paths(net::_FlowNetwork, S::Int, Z::Int, orig::Vector)
    flow = [orig[a] - net.cap[a] for a in eachindex(orig)]  # flow on the forward arcs
    paths = Vector{Vector{Int}}()
    while true
        path = Int[]
        pos = Dict{Int,Int}(S => 0)  # node -> length of the path when reached
        u = S
        while u != Z
            a = 0
            for x in net.off[u]:net.off[u+1]-1
                b = net.arcs[x]
                if isodd(b) && flow[b] > 0
                    a = b
                    break
                end
            end
            a == 0 && return paths
            push!(path, a)
            u = net.head[a]
            if haskey(pos, u)  # cancel the cycle
                k = pos[u]
                for b in path[k+1:end]
                    flow[b] -= 1
                    delete!(pos, net.head[b])
                end
                resize!(path, k)
            end
            pos[u] = length(path)
        end
        for b in path
            flow[b] -= 1
        end
        push!(paths, path)
    end
end

# A temporal path from a walk: cut out the part between two visits of a node.
function _shortcut(es, walk::Vector{Int})
    out = Int[]
    at = Dict{Int,Int}()  # node -> length of out when the node was reached
    at[Int(es[walk[1]].u)] = 0
    for i in walk
        push!(out, i)
        v = Int(es[i].v)
        if haskey(at, v)
            k = at[v]
            for j in out[k+1:end]
                delete!(at, Int(es[j].v))
            end
            resize!(out, k)
        end
        at[v] = length(out)
    end
    return out
end

"""
    temporal_edge_disjoint_paths(g::OrderedEdgeList, s, z, ti = time_interval(g); β = Inf)

A maximum set of edge-disjoint temporal paths from `s` to `z` inside `ti`, as vectors
of edges. Their number equals the minimum number of temporal edges whose removal
disconnects `z` from `s` ([`temporal_min_cut`](@ref)): Menger's theorem holds for
edge-disjoint temporal paths (Berman, 1996), unlike for vertex-disjoint ones (Kempe,
Kleinberg and Kumar, 2000). With a waiting constraint `β` the result are
edge-disjoint `β`-restless walks (which may visit a node several times).

# References

- K. A. Berman. *Vulnerability of scheduled networks and a generalization of Menger's theorem.* Networks 28(3), 1996. [DOI](https://doi.org/10.1002/%28SICI%291097-0037%28199610%2928%3A3%3C125%3A%3AAID-NET1%3E3.0.CO%3B2-P)
- D. Kempe, J. Kleinberg, and A. Kumar. *Connectivity and inference problems for temporal networks.* Journal of Computer and System Sciences 64(4), 2002. [DOI](https://doi.org/10.1006/jcss.2002.1829)
"""
function temporal_edge_disjoint_paths(g::OrderedEdgeList{V,T}, s::Integer, z::Integer, ti=time_interval(g);
                                      β::Real=Inf) where {V,T}
    net, S, Z, arc, ids = _flow_network(g, s, z, ti, nothing, nothing, β)
    orig = copy(net.cap)
    _max_flow!(net, S, Z)
    edge_of = Dict(arc[i] => i for i in ids)
    es = edges(g)
    out = Vector{Vector{TemporalEdge{V,T}}}()
    for p in _flow_paths(net, S, Z, orig)
        walk = [edge_of[a] for a in p if haskey(edge_of, a)]
        push!(out, es[isfinite(β) ? walk : _shortcut(es, walk)])
    end
    return out
end

"""
    temporal_out_disjoint_paths(g::OrderedEdgeList, s, z, ti = time_interval(g))

A maximum set of *out-disjoint* temporal paths from `s` to `z` inside `ti`: no two
of them leave the same node at the same time. Their number equals the minimum number
of node departure times `(v, t)` whose removal (deleting all the edges leaving `v` at
time `t`) disconnects `z` from `s` (Mertzios, Michail and Spirakis, 2019, Theorem 3).
Computed as a maximum flow with one unit-capacity node per departure time.

# References

- G. B. Mertzios, O. Michail, and P. G. Spirakis. *Temporal network optimization subject to connectivity constraints.* Algorithmica 81(4), 2019. [DOI](https://doi.org/10.1007/s00453-018-0478-6), [arXiv](https://arxiv.org/abs/1502.04382)
"""
function temporal_out_disjoint_paths(g::OrderedEdgeList{V,T}, s::Integer, z::Integer, ti=time_interval(g)) where {V,T}
    s, z = _check_node(g, s), _check_node(g, z)
    s == z && throw(ArgumentError("source and sink must be different"))
    a, b = _interval(T, ti)
    ids = _window_edges(g, a, b)
    net, S, Z, arc = _time_expanded(g, ids, s, z, ones(Int, num_edges(g)); departures=true)
    orig = copy(net.cap)
    _max_flow!(net, S, Z)
    edge_of = Dict(arc[i] => i for i in ids)
    es = edges(g)
    return [es[_shortcut(es, [edge_of[a] for a in p if haskey(edge_of, a)])] for p in _flow_paths(net, S, Z, orig)]
end

# Minimum hop temporal path from s to z avoiding the nodes in `removed` (edge
# indices), or nothing. A minimum hop walk is a path.
function _min_hop_path(es::Vector{TemporalEdge{V,T}}, ids::Vector{Int}, n::Int, s::Int, z::Int, removed::BitVector) where {V,T}
    hops = Dict{Int,Int}()
    pred = Dict{Int,Int}()
    best = fill((typemax(Int), 0), n)  # best arrived edge into every node
    pending = BinaryMinHeap{Tuple{T,Int}}()
    k = 1
    L = length(ids)
    goal = 0
    @inbounds while k <= L
        τ = es[ids[k]].t
        j = k
        while j < L && es[ids[j+1]].t == τ
            j += 1
        end
        for phase in 1:2  # edges with transition time 0 first, then the others
            while !isempty(pending) && first(pending)[1] <= τ
                _, i = pop!(pending)
                v = Int(es[i].v)
                (hops[i], i) < best[v] && (best[v] = (hops[i], i))
            end
            block = [ids[x] for x in k:j if iszero(es[ids[x]].tt) == (phase == 1)]
            changed = true
            while changed
                changed = false
                for i in block
                    e = es[i]
                    (removed[e.u] || removed[e.v]) && continue
                    h, p = Int(e.u) == s ? (1, 0) : best[e.u][1] == typemax(Int) ? (typemax(Int), 0) : (best[e.u][1] + 1, best[e.u][2])
                    if h < get(hops, i, typemax(Int))
                        hops[i] = h
                        pred[i] = p
                        changed = true
                        if phase == 1  # usable at once by the other edges of the block
                            v = Int(e.v)
                            (h, i) < best[v] && (best[v] = (h, i))
                        end
                    end
                end
                phase == 2 && break
            end
            for i in block
                haskey(hops, i) || continue
                Int(es[i].v) == z && (goal == 0 || hops[i] < hops[goal]) && (goal = i)
                phase == 2 && push!(pending, (es[i].t + es[i].tt, i))
            end
        end
        k = j + 1
    end
    goal == 0 && return nothing
    path = Int[]
    i = goal
    while i != 0
        pushfirst!(path, i)
        i = pred[i]
    end
    return path
end

"""
    temporal_vertex_separator(g::OrderedEdgeList, s, z, ti = time_interval(g))

A minimum temporal `(s, z)`-separator: a smallest set of nodes other than `s` and `z`
whose removal leaves no temporal path from `s` to `z` inside `ti`, or `nothing` if
an edge goes directly from `s` to `z`. Strict and non-strict paths correspond to
positive and zero transition times.

Finding a minimum separator is NP-hard (Kempe, Kleinberg and Kumar, 2000; Zschoche,
Fluschnik, Molter and Niedermeier, 2020), and the minimum can exceed the maximum
number of vertex-disjoint paths. This exact algorithm is a bounded search tree: a
separator must contain an inner node of every path, so it branches on the inner
nodes of a minimum hop path, for increasing sizes `k` starting from the number of
vertex-disjoint paths found greedily. It runs in `O(ℓ^k m)` time, where `ℓ` is the
largest number of hops of the paths, which is fixed-parameter tractable in `k` and
the number of time stamps for strict paths.

# References

- D. Kempe, J. Kleinberg, and A. Kumar. *Connectivity and inference problems for temporal networks.* Journal of Computer and System Sciences 64(4), 2002. [DOI](https://doi.org/10.1006/jcss.2002.1829)
- P. Zschoche, T. Fluschnik, H. Molter, and R. Niedermeier. *The complexity of finding small separators in temporal graphs.* Journal of Computer and System Sciences 107, 2020. [DOI](https://doi.org/10.1016/j.jcss.2019.07.006), [arXiv](https://arxiv.org/abs/1711.00963)
"""
function temporal_vertex_separator(g::OrderedEdgeList{V,T}, s::Integer, z::Integer, ti=time_interval(g)) where {V,T}
    s, z = _check_node(g, s), _check_node(g, z)
    s == z && throw(ArgumentError("source and sink must be different"))
    a, b = _interval(T, ti)
    es = edges(g)
    ids = [i for i in _window_edges(g, a, b) if es[i].u != es[i].v]
    any(i -> es[i].u == s && es[i].v == z, ids) && return nothing
    n = num_nodes(g)
    removed = falses(n)
    keep = falses(n)  # nodes excluded from the separator in the current branch
    # lower bound: number of paths found greedily, pairwise disjoint in the removable
    # inner nodes (a separator needs one node of each); nothing if a path has no
    # removable inner node
    function lower_bound()
        extra = falses(n)
        k = 0
        while (p = _min_hop_path(es, ids, n, s, z, removed .| extra)) !== nothing
            inner = [Int(es[i].v) for i in p[1:end-1] if !keep[es[i].v]]
            isempty(inner) && return typemax(Int)
            k += 1
            extra[inner] .= true
        end
        return k
    end
    # hitting set branching: branch i removes the i-th removable inner node of a
    # minimum hop path and keeps the previous ones, so every set is generated once
    function branch(k)
        p = _min_hop_path(es, ids, n, s, z, removed)
        p === nothing && return Int[]
        lower_bound() > k && return nothing
        kept = Int[]
        result = nothing
        for i in p[1:end-1]
            v = Int(es[i].v)
            keep[v] && continue
            removed[v] = true
            S = branch(k - 1)
            removed[v] = false
            if S !== nothing
                result = push!(S, v)
                break
            end
            keep[v] = true
            push!(kept, v)
        end
        keep[kept] .= false
        return result
    end
    for k in lower_bound():n-2
        S = branch(k)
        S === nothing || return sort!(S)
    end
    return nothing  # unreachable: removing all other nodes separates s and z
end
