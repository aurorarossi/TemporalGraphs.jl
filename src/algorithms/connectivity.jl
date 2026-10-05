# Temporal reachability between all pairs of nodes.
#
# One forward scan of the edges for a block of sources: every node keeps the bitset of
# the sources of the block that have reached it, 64 sources per machine word. An edge
# (u, v, t, tt) with tt > 0 brings the sources that reached u by time t to v at time
# t + tt: the bitset of u is copied when the edge departs and merged into v when it
# arrives (copies that would not add anything are skipped). Edges with transition
# time 0 and the same time stamp are scanned until nothing changes, as for the
# distances. Blocks of sources are processed in parallel.

# With `A`, the time at which every source of the block first reaches every node is
# recorded in A[source, node] (the bits are set in non-decreasing time order).
@inline function _record!(A, new::UInt64, w::Int, first_source::Int, v::Int, x)
    A === nothing && return nothing
    while new != 0
        r = trailing_zeros(new)
        A[first_source+64*(w-1)+r, v] = x
        new &= new - 1
    end
    return nothing
end

function _reachability_block!(S::Matrix{UInt64}, g::OrderedEdgeList{V,T}, first_source::Int, a::T, b::T,
                              A::Union{Nothing,Matrix{T}}=nothing) where {V,T}
    k, n = size(S)
    fill!(S, 0)
    for s in first_source:min(n, first_source + 64k - 1)
        r = s - first_source
        S[r>>6+1, s] |= UInt64(1) << (r & 63)
    end
    es = edges(g)
    pending = BinaryMinHeap{Tuple{T,Int}}()   # (arrival, snapshot id)
    snaps = Vector{Vector{UInt64}}()
    targets = Int[]
    free = Int[]
    i = 1
    m = length(es)
    @inbounds while i <= m
        t = es[i].t
        t > b && break
        j = i
        while j < m && es[j+1].t == t
            j += 1
        end
        if t >= a
            # arrivals up to time t
            while !isempty(pending) && first(pending)[1] <= t
                x, id = pop!(pending)
                v = targets[id]
                snap = snaps[id]
                for w in 1:k
                    _record!(A, snap[w] & ~S[w, v], w, first_source, v, x)
                    S[w, v] |= snap[w]
                end
                push!(free, id)
            end
            # transition time 0: scan the block until nothing changes (they come first)
            z = i
            while z <= j && iszero(es[z].tt)
                z += 1
            end
            changed = z > i
            while changed
                changed = false
                for x in i:z-1
                    e = es[x]
                    u, v = Int(e.u), Int(e.v)
                    for w in 1:k
                        new = S[w, u] & ~S[w, v]
                        if new != 0
                            _record!(A, new, w, first_source, v, t)
                            S[w, v] |= new
                            changed = true
                        end
                    end
                end
            end
            for x in z:j
                e = es[x]
                e.t + e.tt <= b || continue
                u, v = Int(e.u), Int(e.v)
                adds = false
                for w in 1:k
                    (S[w, u] & ~S[w, v]) != 0 && (adds = true; break)
                end
                adds || continue
                if isempty(free)
                    push!(snaps, Vector{UInt64}(undef, k))
                    push!(targets, 0)
                    id = length(snaps)
                else
                    id = pop!(free)
                end
                snap = snaps[id]
                for w in 1:k
                    snap[w] = S[w, u]
                end
                targets[id] = v
                push!(pending, (e.t + e.tt, id))
            end
        end
        i = j + 1
    end
    @inbounds while !isempty(pending)
        x, id = pop!(pending)
        v = targets[id]
        for w in 1:k
            _record!(A, snaps[id][w] & ~S[w, v], w, first_source, v, x)
            S[w, v] |= snaps[id][w]
        end
    end
    return S
end

"""
    temporal_reachability(g::OrderedEdgeList, ti = time_interval(g))

Reachability matrix of `g`: `R[s, v]` is `true` iff there is a temporal path from `s`
to `v` inside `ti` (every node reaches itself). Edges with transition time 0 can be
chained at the same time stamp, as for the distances.

The sources are processed 64 at a time with bitsets, in one scan of the edges per
block of sources: `O(m n / 64)` time instead of `O(m n)` for one search per source,
and blocks are processed in parallel.
"""
function temporal_reachability(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    n = num_nodes(g)
    R = falses(n, n)
    n == 0 && return R
    words = cld(n, 64)
    k = clamp(cld(words, 2 * Threads.nthreads()), 1, 8)  # words per block
    starts = collect(1:64k:n)
    blocks = tmap(s0 -> _reachability_block!(Matrix{UInt64}(undef, k, n), g, s0, a, b), starts)
    # bits of different columns of a BitMatrix share words: assemble sequentially
    for (s0, S) in zip(starts, blocks), v in 1:n, s in s0:min(n, s0 + 64k - 1)
        r = s - s0
        (S[r>>6+1, v] >> (r & 63)) & 1 == 1 && (R[s, v] = true)
    end
    return R
end

"""
    earliest_arrival_matrix(g::OrderedEdgeList, ti = time_interval(g))

Earliest arrival times between all pairs of nodes: `A[s, v]` is the earliest arrival
time at `v` of a temporal path from `s` inside `ti`, as returned by
[`earliest_arrival_times`](@ref)`(g, s, ti)`: `typemax` (printed as `∞`) for
unreachable nodes and 0 for `v = s`. Returned as a [`TemporalDistances`](@ref) matrix. Computed with the bitset scan of [`temporal_reachability`](@ref),
recording when every source first reaches every node: `O(m n / 64 + n²)` time and
`O(n²)` memory.
"""
function earliest_arrival_matrix(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    n = num_nodes(g)
    A = fill(typemax(T), n, n)
    n == 0 && return TemporalDistances(A)
    k = clamp(cld(cld(n, 64), 2 * Threads.nthreads()), 1, 8)
    # blocks write disjoint rows of A
    tmap(s0 -> (_reachability_block!(Matrix{UInt64}(undef, k, n), g, s0, a, b, A); nothing), collect(1:64k:n))
    for s in 1:n
        A[s, s] = zero(T)
    end
    return TemporalDistances(A)
end

"""
    temporal_flooding_time(g::OrderedEdgeList, s, ti = time_interval(g))

Time needed to *flood* `g` from `s` inside `ti = (a, b)`: if `s` knows a piece of
information at time `a` and every node forwards it on every edge it can use, the time
from `a` until all nodes know it, i.e., the largest earliest arrival time from `s`
minus `a`, or `nothing` if some node cannot be reached. This is the broadcast time of
`s` in the flooding model (every informed node transmits on all its edges).

# References

- A. E. F. Clementi, C. Macci, A. Monti, F. Pasquale, and R. Silvestri. *Flooding time of edge-Markovian evolving graphs.* SIAM Journal on Discrete Mathematics 24, 2010. [DOI](https://doi.org/10.1137/090756053)
"""
function temporal_flooding_time(g::OrderedEdgeList{V,T}, s::Integer, ti=time_interval(g)) where {V,T}
    s = _check_node(g, s)
    a, _ = _interval(T, ti)
    x = _flooding(earliest_arrival_times(g, s, ti), s, a)
    return x == typemax(T) ? nothing : x
end

function _flooding(ea::AbstractVector{T}, s::Int, a::T) where {T}
    worst = a
    for (v, x) in enumerate(ea)
        v == s && continue
        x == typemax(T) && return typemax(T)
        worst = max(worst, x)
    end
    return worst - a
end

"""
    temporal_flooding_times(g::OrderedEdgeList, ti = time_interval(g))

[`temporal_flooding_time`](@ref) of every node, from the [`earliest_arrival_matrix`](@ref),
as a [`TemporalDistances`](@ref) vector: `typemax` (printed as `∞`) for the nodes that
cannot flood `g`.
"""
function temporal_flooding_times(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    a, _ = _interval(T, ti)
    A = parent(earliest_arrival_matrix(g, ti))
    return TemporalDistances([_flooding(view(A, s, :), s, a) for s in 1:num_nodes(g)])
end

"""
    temporal_gossip_time(g::OrderedEdgeList, ti = time_interval(g))

Time needed for *gossiping* in `g` inside `ti = (a, b)`: if every node knows its own
piece of information at time `a` and all nodes forward everything they know on every
edge, the time from `a` until every node knows every piece, i.e., the largest
[`temporal_flooding_time`](@ref), or `nothing` if `g` is not temporally connected
inside `ti`. Equivalently, `a` plus the gossip time is the earliest `t` such that
`g` restricted to `(a, t)` is temporally connected.
"""
function temporal_gossip_time(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    num_nodes(g) <= 1 && return zero(T)
    x = maximum(temporal_flooding_times(g, ti))
    return x == typemax(T) ? nothing : x
end

"""
    is_temporally_connected(g::OrderedEdgeList, ti = time_interval(g))

Whether every node of `g` reaches every other node by a temporal path inside `ti`.
"""
is_temporally_connected(g::OrderedEdgeList, ti=time_interval(g)) = all(temporal_reachability(g, ti))

# Undirected graph joining the pairs that reach each other (both directions, or at
# least one if unilateral).
function _component_graph(R::BitMatrix, unilateral::Bool)
    n = size(R, 1)
    h = Graphs.SimpleGraph(n)
    for u in 1:n, v in u+1:n
        (unilateral ? (R[u, v] || R[v, u]) : (R[u, v] && R[v, u])) && Graphs.add_edge!(h, u, v)
    end
    return h
end

_sorted_components(cs) = sort!([sort!(collect(Int, c)) for c in cs]; by=c -> (-length(c), c))

"""
    temporal_connected_components(g::OrderedEdgeList, ti = time_interval(g); unilateral = false)

The (open) temporal connected components of `g`: the maximal sets of nodes in which
every node reaches every other node by a temporal path inside `ti` (Bhadra and
Ferreira, 2003). The paths may leave the set. With `unilateral = true`, the maximal
sets in which of every two nodes at least one reaches the other.

Unlike in static graphs the components can overlap: they are the maximal cliques of
the graph joining the nodes that reach each other (Costa, Lopes, Marino and Silva,
2023), enumerated with the Bron–Kerbosch algorithm (exponential in the worst case;
finding the largest component is NP-hard). Returns the components sorted by
decreasing size.

# References

- S. Bhadra and A. Ferreira. *Complexity of connected components in evolving graphs and the computation of multicast trees in wireless networks.* ADHOC-NOW, 2003. [DOI](https://doi.org/10.1007/978-3-540-39611-6_23)
- I. L. Costa, R. Lopes, A. Marino, and A. Silva. *On computing large temporal (unilateral) connected components.* Journal of Computer and System Sciences 144, 2024. [DOI](https://doi.org/10.1016/j.jcss.2024.103548), [arXiv](https://arxiv.org/abs/2302.12068)
"""
function temporal_connected_components(g::OrderedEdgeList, ti=time_interval(g); unilateral::Bool=false)
    h = _component_graph(temporal_reachability(g, ti), unilateral)
    return _sorted_components(Graphs.maximal_cliques(h))
end

# Temporal subgraph of g induced by the nodes X, with nodes relabeled 1:length(X).
function _induced(g::OrderedEdgeList{V,T}, X::Vector{Int}, ti) where {V,T}
    a, b = _interval(T, ti)
    id = zeros(Int, num_nodes(g))
    id[X] = eachindex(X)
    es = [TemporalEdge{V,T}(id[e.u], id[e.v], e.t, e.tt) for e in edges(g)
          if id[e.u] != 0 && id[e.v] != 0 && a <= e.t && e.t + e.tt <= b]
    return OrderedEdgeList(length(X), es, (a, b))
end

# Bitsets of nodes as vectors of words.
@inline _bit(x::Vector{UInt64}, v::Int) = (x[(v-1)>>6+1] >> ((v - 1) & 63)) & 1 == 1
@inline _setbit!(x::Vector{UInt64}, v::Int) = (x[(v-1)>>6+1] |= UInt64(1) << ((v - 1) & 63); x)
@inline _clearbit!(x::Vector{UInt64}, v::Int) = (x[(v-1)>>6+1] &= ~(UInt64(1) << ((v - 1) & 63)); x)
function _firstbit(x::Vector{UInt64})
    @inbounds for (w, y) in enumerate(x)
        y != 0 && return 64 * (w - 1) + trailing_zeros(y) + 1
    end
    return 0
end

# Greedy coloring of the nodes of P (Tomita and Seki): nodes in order of increasing
# color; a clique inside P has at most as many nodes as colors.
function _color_sort!(order::Vector{Int}, colors::Vector{Int}, nbr::Vector{Vector{UInt64}}, P::Vector{UInt64})
    empty!(order)
    empty!(colors)
    Q = copy(P)
    Qk = similar(Q)
    k = 0
    while any(!iszero, Q)
        k += 1
        copyto!(Qk, Q)
        while (v = _firstbit(Qk)) != 0
            push!(order, v)
            push!(colors, k)
            _clearbit!(Q, v)
            _clearbit!(Qk, v)
            @inbounds for w in eachindex(Qk)
                Qk[w] &= ~nbr[v][w]
            end
        end
    end
    return k
end

# Maximum clique of the graph h restricted to the nodes X, by branch and bound with
# coloring bounds (Tomita and Seki, 2003); nodes are renumbered by decreasing degree.
function _maximum_clique(h::Graphs.SimpleGraph, X::Vector{Int})
    isempty(X) && return Int[]
    k = length(X)
    pos = Dict(v => i for (i, v) in enumerate(X))
    deg = [count(w -> haskey(pos, w), Graphs.neighbors(h, v)) for v in X]
    perm = sortperm(deg; rev=true)  # node i of the search is X[perm[i]]
    rank = invperm(perm)
    W = cld(k, 64)
    nbr = [zeros(UInt64, W) for _ in 1:k]
    for (i, v) in enumerate(X), w in Graphs.neighbors(h, v)
        j = get(pos, w, 0)
        j == 0 || _setbit!(nbr[rank[i]], rank[j])
    end
    best = Int[]
    for i in 1:k  # greedy initial clique
        all(j -> _bit(nbr[i], j), best) && push!(best, i)
    end
    R = Int[]
    function expand(P::Vector{UInt64})
        order, colors = Int[], Int[]
        _color_sort!(order, colors, nbr, P)
        for x in length(order):-1:1
            length(R) + colors[x] <= length(best) && return
            v = order[x]
            push!(R, v)
            NP = P .& nbr[v]
            if all(iszero, NP)
                length(R) > length(best) && (best = copy(R))
            else
                expand(NP)
            end
            pop!(R)
            _clearbit!(P, v)
        end
    end
    P = zeros(UInt64, W)
    for i in 1:k
        _setbit!(P, i)
    end
    expand(P)
    return sort!([X[perm[i]] for i in best])
end

# Upper bound on the cliques of h inside X: the number of colors of a greedy coloring.
function _clique_bound(h::Graphs.SimpleGraph, X::Vector{Int})
    color = Dict{Int,Int}()
    for v in sort(X; by=v -> -Graphs.degree(h, v))
        used = Set(get(color, w, 0) for w in Graphs.neighbors(h, v))
        c = 1
        while c in used
            c += 1
        end
        color[v] = c
    end
    return isempty(color) ? 0 : maximum(values(color))
end

"""
    largest_temporal_connected_component(g::OrderedEdgeList, ti = time_interval(g);
                                         closed = false, unilateral = false)

A largest temporal connected component of `g` (see
[`temporal_connected_components`](@ref)). With `closed = true` the paths between the
nodes of the component must stay inside it, i.e., the temporal subgraph induced by
the component is temporally connected.

Both problems are NP-hard (Bhadra and Ferreira, 2003), and so is deciding whether a
set is a closed component (Costa et al., 2023). The open version is a maximum clique
of the graph joining the nodes that reach each other, found exactly by branch and
bound with coloring bounds (Tomita and Seki, 2003). The closed version is searched
exactly by branching: if `u` cannot reach `v` inside a candidate set, one of them is
removed; candidate sets are pruned with the same coloring bound. Both are
exponential in the worst case.

# References

- S. Bhadra and A. Ferreira. *Complexity of connected components in evolving graphs and the computation of multicast trees in wireless networks.* ADHOC-NOW, 2003. [DOI](https://doi.org/10.1007/978-3-540-39611-6_23)
- I. L. Costa, R. Lopes, A. Marino, and A. Silva. *On computing large temporal (unilateral) connected components.* Journal of Computer and System Sciences 144, 2024. [DOI](https://doi.org/10.1016/j.jcss.2024.103548), [arXiv](https://arxiv.org/abs/2302.12068)
- E. Tomita and T. Seki. *An efficient branch-and-bound algorithm for finding a maximum clique.* DMTCS, 2003. [DOI](https://doi.org/10.1007/3-540-45066-1_22)
"""
function largest_temporal_connected_component(g::OrderedEdgeList, ti=time_interval(g);
                                              closed::Bool=false, unilateral::Bool=false)
    n = num_nodes(g)
    n == 0 && return Int[]
    h = _component_graph(temporal_reachability(g, ti), unilateral)
    closed || return _maximum_clique(h, collect(1:n))
    best = [1]
    seen = Set{Vector{Int}}()
    function search(X::Vector{Int})
        length(X) <= length(best) && return
        X in seen && return
        push!(seen, X)
        _clique_bound(h, X) <= length(best) && return
        R = temporal_reachability(_induced(g, X, ti), ti)
        k = length(X)
        # the node with the most unreachable partners inside X, and one of its partners
        bad = zeros(Int, k)
        for i in 1:k, j in 1:k
            i == j && continue
            if unilateral ? !(R[i, j] || R[j, i]) : !R[i, j]
                bad[i] += 1
                bad[j] += 1
            end
        end
        if all(iszero, bad)
            best = X
            return
        end
        i = argmax(bad)
        j = findfirst(j -> j != i && (unilateral ? !(R[i, j] || R[j, i]) : !(R[i, j] && R[j, i])), 1:k)
        search(deleteat!(copy(X), i))
        search(deleteat!(copy(X), j))
    end
    # every closed component lies inside a connected component of h
    for c in sort!(Graphs.connected_components(h); by=length, rev=true)
        search(sort(c))
    end
    return best
end
