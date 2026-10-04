# Temporal walks and paths under a waiting time constraint β ("Δ-restless"), following
# Casteigts, Himmel, Molter and Zschoche, "Finding temporal paths under waiting time
# constraints" (Algorithmica 2021): the edge following an edge e must depart in
# [arr(e), arr(e) + β], where arr(e) = t + tt; the departure from the source is free.
#
# Optimal restless walks are found in polynomial time by one scan of the edges by
# departure time: for every edge f the best cost of a restless walk from s ending
# with f is obtained from the incoming edges of f.u that arrived in [f.t - β, f.t]
# (a sliding window per node, kept as a monotone deque of costs). Edges with
# transition time 0 and the same time stamp can follow each other: their costs are
# relaxed until nothing changes. Restless paths (no repeated node) are NP-hard to
# find; they are searched exactly by backtracking.

_rl_cost_type(::MinimumHops, ::Type{T}) where {T} = Int
_rl_cost_type(::DistanceType, ::Type{T}) where {T} = T
@inline _rl_first(::MinimumHops, e) = 1
@inline _rl_first(::MinimumTransitionTimes, e) = e.tt
@inline _rl_first(::EarliestArrival, e) = zero(e.t)
@inline _rl_first(::Union{LatestDeparture,Fastest}, e) = zero(e.t) - e.t
@inline _rl_extend(::MinimumHops, c, e) = c + 1
@inline _rl_extend(::MinimumTransitionTimes, c, e) = c + e.tt
@inline _rl_extend(::Union{EarliestArrival,LatestDeparture,Fastest}, c, e) = c

const _RestlessCriterion = Union{EarliestArrival,LatestDeparture,Fastest,MinimumHops,MinimumTransitionTimes}

# Edges inside ti, their incoming lists sorted by arrival (transition time 0 last)
# and the best cost C of a restless walk from s ending with every edge.
struct _RestlessScan{V,T,C}
    es::Vector{TemporalEdge{V,T}}
    inoff::Vector{Int}
    ins::Vector{Int}
    C::Vector{C}
end

function _restless_costs(g::OrderedEdgeList{V,T}, s::Int, dt, β, ti; skip::Union{Nothing,AbstractVector{Bool}}=nothing,
                         start_after=nothing) where {V,T}
    a, b = _interval(T, ti)
    C = _rl_cost_type(dt, T)
    inf = typemax(C)
    es = [e for e in edges(g) if a <= e.t && e.t + e.tt <= b &&
                                 (skip === nothing || !(skip[e.u] || skip[e.v]))]
    n, m = num_nodes(g), length(es)
    βT = convert(T, β)
    inoff = zeros(Int, n + 1)
    for e in es
        inoff[e.v+1] += 1
    end
    inoff[1] = 1
    for u in 1:n
        inoff[u+1] += inoff[u]
    end
    ins = Vector{Int}(undef, m)
    pos = inoff[1:n]
    for i in sortperm([(e.t + e.tt, iszero(e.tt)) for e in es])
        v = es[i].v
        ins[pos[v]] = i
        pos[v] += 1
    end
    cost = fill(inf, m)
    inptr = inoff[1:n]
    head = inoff[1:n]
    tail = inoff[1:n] .- 1
    dq_arr = Vector{T}(undef, m)
    dq_cost = Vector{C}(undef, m)
    # minimum cost of the walks that an edge leaving u at time τ can extend
    function window(u::Int, τ, zero_block::Bool)
        k = inptr[u]
        @inbounds while k < inoff[u+1]
            e = ins[k]
            x = es[e].t + es[e].tt
            (x < τ || (x == τ && !(zero_block && iszero(es[e].tt)))) || break
            k += 1
            c = cost[e]
            c == inf && continue
            t = tail[u]
            while t >= head[u] && dq_cost[t] >= c
                t -= 1
            end
            t += 1
            dq_arr[t] = x
            dq_cost[t] = c
            tail[u] = t
        end
        inptr[u] = k
        h = head[u]
        @inbounds while h <= tail[u] && dq_arr[h] < τ - βT
            h += 1
        end
        head[u] = h
        return h <= tail[u] ? @inbounds(dq_cost[h]) : inf
    end
    # a walk can start with an edge leaving s (at or after `start_after` if given)
    starts(e) = Int(e.u) == s && (start_after === nothing || start_after <= e.t <= start_after + βT)
    function candidate(e, p)
        c = p == inf ? inf : _rl_extend(dt, p, e)
        starts(e) && (c = min(c, _rl_first(dt, e)))
        return c
    end
    bestin = fill(inf, n)
    i = 1
    @inbounds while i <= m
        f = es[i]
        if iszero(f.tt)
            j = i
            while j < m && es[j+1].t == f.t && iszero(es[j+1].tt)
                j += 1
            end
            for x in i:j
                cost[x] = candidate(es[x], window(Int(es[x].u), f.t, true))
            end
            changed = true
            while changed  # relax the edges of the block until nothing changes
                changed = false
                for x in i:j
                    v = Int(es[x].v)
                    cost[x] < bestin[v] && (bestin[v] = cost[x])
                end
                for x in i:j
                    p = bestin[Int(es[x].u)]
                    p == inf && continue
                    c = _rl_extend(dt, p, es[x])
                    if c < cost[x]
                        cost[x] = c
                        changed = true
                    end
                end
            end
            for x in i:j
                bestin[Int(es[x].v)] = inf
            end
            i = j + 1
        else
            cost[i] = candidate(f, window(Int(f.u), f.t, false))
            i += 1
        end
    end
    return _RestlessScan{V,T,C}(es, inoff, ins, cost)
end

function _restless_distances(g::OrderedEdgeList{V,T}, s::Int, dt::_RestlessCriterion, β, ti) where {V,T}
    a, b = _interval(T, ti)
    sc = _restless_costs(g, s, dt, β, (a, b))
    D = distance_eltype(g, dt)
    n = num_nodes(g)
    C = eltype(sc.C)
    best = fill(typemax(D), n)  # minimized target value of every node
    @inbounds for (i, e) in enumerate(sc.es)
        c = sc.C[i]
        c == typemax(C) && continue
        v = Int(e.v)
        x = dt isa EarliestArrival ? e.t + e.tt :
            dt isa Fastest ? e.t + e.tt + c : c
        x < best[v] && (best[v] = x)
    end
    if dt isa LatestDeparture  # best = minus the latest departure
        dist = [x == typemax(D) ? typemin(D) : -x for x in best]
        dist[s] = b
        return dist
    end
    best[s] = zero(D)
    return best
end

"""
    restless_path(g::OrderedEdgeList, s, z, β, ti = time_interval(g))

A *β-restless temporal path* from `s` to `z` (Casteigts, Himmel, Molter and Zschoche,
2021): a temporal path that visits every node at most once and waits at most `β` at
every intermediate node, i.e., every edge departs in `[a, a + β]` where `a` is the
arrival time of the previous edge (the departure from `s` is free). Returns the edges
of a path, or `nothing` if there is none.

Unlike restless *walks* (see [`temporal_distances`](@ref) with `β`), finding a
restless path is NP-hard for every `β ≥ 1`, so this is an exact exponential-time
search: a depth-first search over the edges, pruned by checking that `z` is still
reachable by a restless walk avoiding the visited nodes (a polynomial relaxation).
It is fast when restless walks rarely revisit nodes.

# References

- A. Casteigts, A.-S. Himmel, H. Molter, and P. Zschoche. *Finding temporal paths under waiting time constraints.* Algorithmica 83(9), 2021. [DOI](https://doi.org/10.1007/s00453-021-00831-w), [arXiv](https://arxiv.org/abs/1909.06437)
"""
function restless_path(g::OrderedEdgeList{V,T}, s::Integer, z::Integer, β::Real, ti=time_interval(g)) where {V,T}
    s, z = _check_node(g, s), _check_node(g, z)
    β >= 0 || throw(ArgumentError("β must be non-negative"))
    a, b = _interval(T, ti)
    s == z && return TemporalEdge{V,T}[]
    n = num_nodes(g)
    es = [e for e in edges(g) if a <= e.t && e.t + e.tt <= b && e.u != e.v]
    out = [Int[] for _ in 1:n]  # outgoing edges sorted by departure
    for (i, e) in enumerate(es)
        push!(out[e.u], i)
    end
    visited = falses(n)
    visited[s] = true
    path = Int[]
    sub = OrderedEdgeList(n, es, (a, b))
    # can z still be reached from u, arriving there at time `arr`, avoiding visited nodes?
    function feasible(u, arr)
        skip = copy(visited)
        skip[u] = false
        sc = _restless_costs(sub, u, EarliestArrival(), β, (a, b); skip=skip, start_after=arr)
        return any(i -> sc.es[i].v == z && sc.C[i] != typemax(eltype(sc.C)), eachindex(sc.es))
    end
    function search(u, arr)
        for i in out[u]
            e = es[i]
            if arr !== nothing
                e.t < arr && continue
                e.t > arr + β && break
            end
            v = Int(e.v)
            visited[v] && continue
            push!(path, i)
            v == z && return true
            visited[v] = true
            if feasible(v, e.t + e.tt) && search(v, e.t + e.tt)
                return true
            end
            visited[v] = false
            pop!(path)
        end
        return false
    end
    feasible(s, nothing) || return nothing
    return search(s, nothing) ? es[path] : nothing
end
