# Temporal spanning trees and spanners.
#
# Spanning trees follow Huang, Fu and Liu, "Minimum spanning trees in temporal graphs"
# (SIGMOD 2015): a temporal spanning tree rooted at r contains, for every node reachable
# from r, one incoming edge, such that the tree path from r to every node is a temporal
# path. MST_a minimizes the arrival times (O(m)), MST_w the total weight (NP-hard,
# approximated through a reduction to the Directed Steiner Tree problem).
#
# Spanners keep the reachability between all pairs of nodes with a subset of the
# edges. Minimum spanners are NP-hard to find (and to approximate, Axiotis and Fotakis,
# 2016); temporal cliques admit spanners of size O(n log n) (Casteigts, Peters and
# Schoeters, 2021), computed with the algorithm of Angrick et al. (ESA 2024).

# Arrival-time optimal spanning tree of the edges `ids` (Algorithm 2 of Huang et al.,
# correct for any transition times): the parent edge of every node and its arrival.
function _earliest_arrival_tree(es::Vector{TemporalEdge{V,T}}, ids::Vector{Int}, n::Int, r::Int, a::T) where {V,T}
    out = [Int[] for _ in 1:n]
    for i in Iterators.reverse(ids)  # the edges are sorted by time: non-increasing departures
        push!(out[es[i].u], i)
    end
    pos = ones(Int, n)
    arrival = fill(typemax(T), n)
    parent = zeros(Int, n)
    stack = [(r, a, 0)]
    while !isempty(stack)
        v, x, e = pop!(stack)
        x < arrival[v] || continue
        arrival[v] = x
        parent[v] = e
        o = out[v]
        while pos[v] <= length(o) && x <= es[o[pos[v]]].t
            f = o[pos[v]]
            push!(stack, (Int(es[f].v), es[f].t + es[f].tt, f))
            pos[v] += 1
        end
    end
    return parent, arrival
end

"""
    earliest_arrival_tree(g::OrderedEdgeList, r, ti = time_interval(g))

A temporal spanning tree of the nodes reachable from `r` inside `ti` in which every
node is reached at its earliest arrival time: *MST_a* of Huang, Fu and Liu (2015).
Returns the tree edges, one entering every reachable node other than `r`, sorted by
time. Computed in `O(m)` time (after sorting the edges of every node) with the
stack-based algorithm of the paper, which is correct also for transition times 0.

# References

- S. Huang, A. W.-C. Fu, and R. Liu. *Minimum spanning trees in temporal graphs.* SIGMOD, 2015. [DOI](https://doi.org/10.1145/2723372.2723717)
"""
function earliest_arrival_tree(g::OrderedEdgeList{V,T}, r::Integer, ti=time_interval(g)) where {V,T}
    r = _check_node(g, r)
    a, b = _interval(T, ti)
    es = edges(g)
    parent, _ = _earliest_arrival_tree(es, _window_edges(g, a, b), num_nodes(g), r, a)
    return es[sort!([p for p in parent if p != 0])]
end

# Shortest paths with non-negative weights from `src` (forward) in a graph given by
# adjacency lists of (neighbor, weight, label).
function _dijkstra(adj::Vector{Vector{Tuple{Int,Float64,Int}}}, src::Int)
    N = length(adj)
    dist = fill(Inf, N)
    pred = zeros(Int, N)   # previous node
    plab = zeros(Int, N)   # label of the arc from the previous node
    dist[src] = 0.0
    heap = BinaryMinHeap{Tuple{Float64,Int}}([(0.0, src)])
    while !isempty(heap)
        d, u = pop!(heap)
        d > dist[u] && continue
        for (v, w, l) in adj[u]
            if d + w < dist[v]
                dist[v] = d + w
                pred[v] = u
                plab[v] = l
                push!(heap, (dist[v], v))
            end
        end
    end
    return dist, pred, plab
end

"""
    minimum_weight_spanning_tree(g::OrderedEdgeList, r, ti = time_interval(g);
                                 weights = nothing, level = 2)

An approximate *MST_w* of Huang, Fu and Liu (2015): a temporal spanning tree of the
nodes reachable from `r` inside `ti` (one incoming edge per node, every tree path a
temporal path) of small total weight, where edge `i` weighs `weights[i]` (its
transition time by default). Returns `(edges, weight)`.

Finding a minimum weight tree is NP-hard and MAX-SNP-hard. The algorithm of the paper
transforms the graph into an instance of the Directed Steiner Tree problem (one copy
of every node per arrival time) and solves it with the level-`level` greedy algorithm
of Charikar et al. (1999), with the improvements of the paper (density bounds); the
approximation ratio is ``i^2 (i-1) k^{1/i}`` for level ``i ≥ 2`` and ``k`` terminals.
`level = 1` connects every node by a cheapest temporal path. Levels 1 and 2 are
supported; level 2 takes `O(n m log m)` time for the shortest paths plus the greedy
steps.

# References

- S. Huang, A. W.-C. Fu, and R. Liu. *Minimum spanning trees in temporal graphs.* SIGMOD, 2015. [DOI](https://doi.org/10.1145/2723372.2723717)
- M. Charikar, C. Chekuri, T. Cheung, Z. Dai, A. Goel, S. Guha, and M. Li. *Approximation algorithms for directed Steiner problems.* Journal of Algorithms 33(1), 1999. [DOI](https://doi.org/10.1006/jagm.1999.1042)
"""
function minimum_weight_spanning_tree(g::OrderedEdgeList{V,T}, r::Integer, ti=time_interval(g);
                                      weights::Union{Nothing,AbstractVector{<:Real}}=nothing, level::Integer=2) where {V,T}
    r = _check_node(g, r)
    level in (1, 2) || throw(ArgumentError("level must be 1 or 2"))
    a, b = _interval(T, ti)
    es = edges(g)
    n = num_nodes(g)
    w = weights === nothing ? [Float64(e.tt) for e in es] : Float64.(weights)
    length(w) == length(es) || throw(ArgumentError("weights must have one entry per edge"))
    all(>=(0), w) || throw(ArgumentError("weights must be non-negative"))
    ids = _window_edges(g, a, b)
    _, arrival = _earliest_arrival_tree(es, ids, n, r, a)
    reach = [v != r && arrival[v] != typemax(T) for v in 1:n]
    # transformed graph: copies of every node per arrival time, a terminal per node
    times = [T[] for _ in 1:n]
    times[r] = [a]
    for i in ids
        e = es[i]
        (reach[e.v] && e.v != r) && push!(times[e.v], e.t + e.tt)
    end
    off = ones(Int, n + 1)
    for v in 1:n
        unique!(sort!(times[v]))
        off[v+1] = off[v] + length(times[v])
    end
    C = off[n+1] - 1
    terminal = zeros(Int, n)
    N = C
    for v in 1:n
        reach[v] && (terminal[v] = (N += 1))
    end
    adj = [Tuple{Int,Float64,Int}[] for _ in 1:N]
    for v in 1:n, k in off[v]:off[v+1]-1
        k < off[v+1] - 1 && push!(adj[k], (k + 1, 0.0, 0))
        k == off[v+1] - 1 && terminal[v] != 0 && push!(adj[k], (terminal[v], 0.0, 0))
    end
    for i in ids
        e = es[i]
        (e.v == r || !reach[e.v] || (e.u != r && !reach[e.u])) && continue
        tu = times[e.u]
        x = searchsortedlast(tu, e.t)  # latest copy of u by time t
        x == 0 && continue
        y = searchsortedfirst(times[e.v], e.t + e.tt)
        push!(adj[off[e.u]+x-1], (off[e.v] + y - 1, w[i], i))
    end
    root = off[r]
    X = [terminal[v] for v in 1:n if terminal[v] != 0]
    isempty(X) && return (edges=TemporalEdge{V,T}[], weight=0.0)
    dr, pr, lr = _dijkstra(adj, root)
    selected = Set{Int}()
    addpath!(pred, plab, from, to) = (x = to; while x != from; plab[x] != 0 && push!(selected, plab[x]); x = pred[x]; end)
    if level == 1
        for x in X
            addpath!(pr, lr, root, x)
        end
    else
        # distances from every node to every terminal (on the reversed graph)
        radj = [Tuple{Int,Float64,Int}[] for _ in 1:N]
        for u in 1:N, (v, wt, l) in adj[u]
            push!(radj[v], (u, wt, l))
        end
        toX = Dict{Int,Tuple{Vector{Float64},Vector{Int},Vector{Int}}}()
        for x in X
            toX[x] = _dijkstra(radj, x)
        end
        left = Set(X)
        τ = fill(-Inf, N)   # lower bounds of the densities (they only grow)
        buf = Tuple{Float64,Int}[]
        while !isempty(left)
            best = (Inf, 0, 0)  # density, center v, number of terminals
            for v in sortperm(τ)
                τ[v] >= best[1] && break
                if !isfinite(dr[v])
                    τ[v] = Inf
                    continue
                end
                # B_1: the cheapest prefix of the terminals sorted by distance from v
                empty!(buf)
                for x in left
                    d = toX[x][1][v]
                    isfinite(d) && push!(buf, (d, x))
                end
                if isempty(buf)
                    τ[v] = Inf
                    continue
                end
                sort!(buf)
                total = dr[v]
                bd, bj = Inf, 0
                for (j, (d, _)) in enumerate(buf)
                    total += d
                    total / j < bd && ((bd, bj) = (total / j, j))
                end
                τ[v] = bd
                bd < best[1] && (best = (bd, v, bj))
            end
            best[2] == 0 && break
            _, v, j = best
            empty!(buf)
            for x in left
                d = toX[x][1][v]
                isfinite(d) && push!(buf, (d, x))
            end
            sort!(buf)
            addpath!(pr, lr, root, v)
            for (_, x) in buf[1:j]
                # follow the shortest path from v to x (successors in the reversed search)
                _, succ, slab = toX[x]
                y = v
                while y != x
                    slab[y] != 0 && push!(selected, slab[y])
                    y = succ[y]
                end
                delete!(left, x)
            end
        end
    end
    # a spanning tree inside the selected edges (its weight is at most theirs)
    parent, _ = _earliest_arrival_tree(es, sort!(collect(selected)), n, r, a)
    tree = sort!([p for p in parent if p != 0])
    return (edges=es[tree], weight=sum(w[tree]; init=0.0))
end

"""
    temporal_spanner(g::OrderedEdgeList, ti = time_interval(g); minimal = false)

A *temporal spanner* of `g`: a subset of the edges inside `ti` with the same
reachability between all pairs of nodes. It is the union of an earliest arrival tree
(see [`earliest_arrival_tree`](@ref)) from every node, with at most `n (n - 1)` edges,
computed in `O(n m)` time. With `minimal = true` edges are then removed one at a time
as long as the reachability does not change (`O(k m n / 64)` time for `k` edges of
the union), so that the result is inclusion-minimal.

Minimum spanners are NP-hard to compute and to approximate (Axiotis and Fotakis,
2016), and some temporally connected graphs need `Ω(n²)` edges; for temporal cliques
see [`temporal_clique_spanner`](@ref).

# References

- K. Axiotis and D. Fotakis. *On the size and the approximability of minimum temporally connected subgraphs.* ICALP, 2016. [DOI](https://doi.org/10.4230/LIPIcs.ICALP.2016.149), [arXiv](https://arxiv.org/abs/1602.06411)
"""
function temporal_spanner(g::OrderedEdgeList{V,T}, ti=time_interval(g); minimal::Bool=false) where {V,T}
    a, b = _interval(T, ti)
    es = edges(g)
    n = num_nodes(g)
    ids = _window_edges(g, a, b)
    keep = falses(length(es))
    parents = tmap(r -> _earliest_arrival_tree(es, ids, n, r, a)[1], 1:n)
    for p in parents, i in p
        i != 0 && (keep[i] = true)
    end
    S = findall(keep)
    if minimal
        R = temporal_reachability(g, (a, b))
        for i in reverse(copy(S))
            trial = filter(!=(i), S)
            temporal_reachability(OrderedEdgeList(n, es[trial], (a, b)), (a, b)) == R && (S = trial)
        end
    end
    return es[S]
end

"""
    temporal_clique_spanner(g::OrderedEdgeList)

A temporal spanner with `O(n log n)` edges of a *temporal clique*: an undirected
temporal graph with exactly one time label on every pair of nodes, given with both
directions of every edge (or one of them) and transition time 0 (non-strict paths).
Casteigts, Peters and Schoeters (2021) proved that such spanners always exist; this
is the simpler algorithm of Angrick et al. (ESA 2024): the clique is reduced to a bipartite clique, nodes are removed
when they can delegate their journeys to another node ("dismountability"), and the
remaining instance is split in halves recursively. Returns the edges of the spanner
(both directions of every chosen pair), in `O(n²)` time per level of the recursion.

# References

- A. Casteigts, J. G. Peters, and J. Schoeters. *Temporal cliques admit sparse spanners.* Journal of Computer and System Sciences 121, 2021. [DOI](https://doi.org/10.1016/j.jcss.2021.04.004), [arXiv](https://arxiv.org/abs/1810.00104)
- S. Angrick, B. Bals, T. Friedrich, H. Gawendowicz, N. Hastrich, N. Klodt, P. Lenzner, J. Schmidt, G. Skretas, and A. Wells. *How to reduce temporal cliques to find sparse spanners.* ESA, 2024. [DOI](https://doi.org/10.4230/LIPIcs.ESA.2024.11), [arXiv](https://arxiv.org/abs/2402.13624)
"""
function temporal_clique_spanner(g::OrderedEdgeList{V,T}) where {V,T}
    n = num_nodes(g)
    es = edges(g)
    all(e -> iszero(e.tt), es) || throw(ArgumentError("a temporal clique has transition times 0"))
    label = fill(typemin(T), n, n)
    for e in es
        e.u == e.v && continue
        x = label[e.u, e.v]
        (x == typemin(T) || x == e.t) || throw(ArgumentError("every pair of nodes needs a single time label"))
        label[e.u, e.v] = label[e.v, e.u] = e.t
    end
    all(label[u, v] != typemin(T) for u in 1:n, v in 1:n if u != v) ||
        throw(ArgumentError("the graph is not a temporal clique"))
    n <= 3 && return copy(es)
    # injective labels refining the order (safe for non-strict paths), and the
    # bipartite clique A = V, B = copies of V with the copy edges {u, u'} last
    pairs = sort!([(label[u, v], u, v) for u in 1:n for v in u+1:n])
    λ = zeros(Int, n, n)
    for (k, (_, u, v)) in enumerate(pairs)
        λ[u, v] = λ[v, u] = k
    end
    for u in 1:n
        λ[u, u] = length(pairs) + 1
    end
    # λ[a, b]: label of the edge between a ∈ A and the copy b' ∈ B
    S = Set{Tuple{Int,Int}}()  # edges (a, b) of the bi-spanner
    # Dismount nodes of A and B as long as possible; returns the remaining sets.
    function dismount(A::Vector{Int}, B::Vector{Int})
        A, B = copy(A), copy(B)
        changed = true
        while changed
            changed = false
            # a can delegate to a' if a reaches π⁻(a') before a' leaves it
            if length(A) > 1
                first_b = Dict(x => B[argmin([λ[x, y] for y in B])] for x in A)
                for (k, x) in enumerate(A), y in A
                    y == x && continue
                    p = first_b[y]
                    if λ[x, p] <= λ[y, p]
                        push!(S, (x, p), (y, p))
                        deleteat!(A, k)
                        changed = true
                        break
                    end
                end
                changed && continue
            end
            # b can delegate to b' if π⁺(b') reaches b after b' is reached
            if length(B) > 1
                last_a = Dict(y => A[argmax([λ[x, y] for x in A])] for y in B)
                for (k, y) in enumerate(B), x in B
                    x == y && continue
                    p = last_a[x]
                    if λ[p, x] <= λ[p, y]
                        push!(S, (p, x), (p, y))
                        deleteat!(B, k)
                        changed = true
                        break
                    end
                end
            end
        end
        return A, B
    end
    function bispanner(A::Vector{Int}, B::Vector{Int})
        if length(A) == 1
            for y in B
                push!(S, (A[1], y))
            end
            return
        end
        h = length(A) ÷ 2
        for half in (A[1:h], A[h+1:end])
            A2, B2 = dismount(half, B)
            bispanner(A2, B2)
        end
    end
    bispanner(collect(1:n), collect(1:n))
    chosen = Set{Tuple{Int,Int}}()
    for (x, y) in S
        x == y || push!(chosen, minmax(x, y))
    end
    return [e for e in es if minmax(Int(e.u), Int(e.v)) in chosen]
end
