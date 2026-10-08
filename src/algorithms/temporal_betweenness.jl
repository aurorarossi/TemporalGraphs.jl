# Temporal betweenness of nodes for optimal temporal walks, in the framework of
# Brunelli, Crescenzi and Viennot, "Making temporal betweenness computation faster
# and restless" (KDD 2024), implemented from the definitions of the paper.
#
# Walks are strict (transition times are positive) and, with a waiting constraint β,
# restless: the edge following an edge e must depart in [arr(e), arr(e) + β].
# An optimality criterion is given by a cost structure (costs of walks, combined
# edge by edge) and a target cost (which walks to a node are optimal).
#
# For every source s the algorithm makes three passes over the edges:
# 1. forward, by departure time: for every edge f the minimum cost C[f] of the
#    walks from s ending with f, their number Σ[f], and the minimum cost P[f] of the
#    walks that f can extend (its best predecessors);
# 2. target: for every node v the optimal target cost and the number of optimal walks;
# 3. backward, by decreasing departure time: the s-betweenness of every edge,
#    b[e] = Σ[e] · Σ_{f successor of e} b[f] / Σ[f] + [e ends an optimal walk] Σ[e] / σ*[v],
#    where f is a successor of e iff f extends e and C[e] = P[f] (e is a best
#    predecessor of f) and C[f] = C[e] ⊕ γ(f).
# The betweenness of v sums b over the edges entering v, minus 1 for every source
# that reaches v (Brunelli et al., Fact 3). Walks are counted with multiplicity.
#
# Without waiting constraint the best predecessor cost P is non-increasing along the
# outgoing edges of a node, so the successors of an edge are the end of a block of
# equal P and their contributions are read from sums inside the block in O(1). With a
# waiting constraint the predecessors of an edge form a sliding window of arrival
# times (a monotone deque gives their minimum and count), and the successors are
# found in buckets of the outgoing edges grouped by P. All sums only add non-negative
# values (no differences of prefix sums), so they are numerically stable.

const _WalkCriterion = Union{MinimumHops,EarliestArrival,LatestDeparture,Fastest,
                             ShortestForemost,ShortestLatest,ShortestFastest}

# cost structures (Brunelli et al., Tables 8-10)
_bw_cost_type(::Union{MinimumHops,ShortestForemost,EarliestArrival}, ::Type{T}) where {T} = Int
_bw_cost_type(::Union{LatestDeparture,Fastest}, ::Type{T}) where {T} = T
_bw_cost_type(::Union{ShortestLatest,ShortestFastest}, ::Type{T}) where {T} = Tuple{T,Int}
_bw_target_type(::MinimumHops, ::Type{T}) where {T} = Int
_bw_target_type(::Union{EarliestArrival,LatestDeparture,Fastest}, ::Type{T}) where {T} = T
_bw_target_type(::Union{ShortestForemost,ShortestLatest,ShortestFastest}, ::Type{T}) where {T} = Tuple{T,Int}

_bw_inf(::Type{T}) where {T<:Real} = typemax(T)
_bw_inf(::Type{Tuple{T,Int}}) where {T} = (typemax(T), typemax(Int))

# cost γ(e) of the walk ⟨e⟩ (latest criteria: minus the departure time)
@inline _bw_first(::Union{MinimumHops,ShortestForemost}, e) = 1
@inline _bw_first(::EarliestArrival, e) = 0
@inline _bw_first(::Union{LatestDeparture,Fastest}, e) = zero(e.t) - e.t
@inline _bw_first(::Union{ShortestLatest,ShortestFastest}, e) = (zero(e.t) - e.t, 1)
# cost c ⊕ γ(e) of a walk with cost c extended by e
@inline _bw_extend(::Union{MinimumHops,ShortestForemost}, c, e) = c + 1
@inline _bw_extend(::EarliestArrival, c, e) = 0
@inline _bw_extend(::Union{LatestDeparture,Fastest}, c, e) = c
@inline _bw_extend(::Union{ShortestLatest,ShortestFastest}, c, e) = (c[1], c[2] + 1)
# target cost of a walk with cost c ending with e
@inline _bw_target(::MinimumHops, c, e) = c
@inline _bw_target(::EarliestArrival, c, e) = e.t + e.tt
@inline _bw_target(::LatestDeparture, c, e) = c
@inline _bw_target(::Fastest, c, e) = e.t + e.tt + c
@inline _bw_target(::ShortestForemost, c, e) = (e.t + e.tt, c)
@inline _bw_target(::ShortestLatest, c, e) = c
@inline _bw_target(::ShortestFastest, c, e) = (e.t + e.tt + c[1], c[2])

# Edge lists shared by all sources: the edges sorted by departure time (their index
# is the edge id), the outgoing edges of each node sorted by departure time, the
# incoming edges sorted by arrival time (edges with transition time 0 last among equal
# arrivals), and for every edge the range of outgoing edges of its head that can
# follow it. `restless` selects the sliding window and bucket algorithm, used with a
# waiting constraint (`expire`) and for graphs with transition times 0 (non-strict
# walks).
struct _WalkGraph{V,T}
    n::Int
    es::Vector{TemporalEdge{V,T}}
    outoff::Vector{Int}
    outs::Vector{Int}        # outs[g] = id of the g-th outgoing edge (global position g)
    outpos::Vector{Int}      # global out position of each edge
    inoff::Vector{Int}
    ins::Vector{Int}         # incoming edges of each node, sorted by arrival
    first_succ::Vector{Int}  # first out position of e.v departing at or after arr(e)
    last_succ::Vector{Int}   # last out position of e.v departing at or before arr(e) + β
    restless::Bool
    expire::Bool
    has_zero::Bool
    β::T
end

function _walk_graph(g::OrderedEdgeList{V,T}, ti, β) where {V,T}
    a, b = _interval(T, ti)
    es = [e for e in edges(g) if a <= e.t && e.t + e.tt <= b]
    any(e -> e.tt < 0, es) && throw(ArgumentError("transition times must be non-negative"))
    expire = β !== nothing && isfinite(β)
    βT = expire ? convert(T, β) : zero(T)
    expire && βT < 0 && throw(ArgumentError("the waiting constraint β must be non-negative"))
    has_zero = any(e -> iszero(e.tt), es)
    restless = expire || has_zero
    n = num_nodes(g)
    m = length(es)
    outoff = zeros(Int, n + 1)
    inoff = zeros(Int, n + 1)
    for e in es
        outoff[e.u+1] += 1
        inoff[e.v+1] += 1
    end
    outoff[1] = inoff[1] = 1
    for i in 2:n+1
        outoff[i] += outoff[i-1]
        inoff[i] += inoff[i-1]
    end
    outs = Vector{Int}(undef, m)
    outpos = Vector{Int}(undef, m)
    pos = outoff[1:n]
    for (i, e) in enumerate(es)  # stable: sorted by departure inside each list
        outs[pos[e.u]] = i
        outpos[i] = pos[e.u]
        pos[e.u] += 1
    end
    byarr = sortperm([(e.t + e.tt, iszero(e.tt)) for e in es])
    ins = Vector{Int}(undef, m)
    pos = inoff[1:n]
    for i in byarr
        v = es[i].v
        ins[pos[v]] = i
        pos[v] += 1
    end
    first_succ = Vector{Int}(undef, m)
    last_succ = Vector{Int}(undef, m)
    for (i, e) in enumerate(es)
        lo, hi = outoff[e.v], outoff[e.v+1] - 1
        arr = e.t + e.tt
        first_succ[i] = _first_departure(es, outs, lo, hi, arr, false)
        last_succ[i] = expire ? _first_departure(es, outs, lo, hi, arr + βT, true) - 1 : hi
    end
    return _WalkGraph{V,T}(n, es, outoff, outs, outpos, inoff, ins, first_succ, last_succ, restless, expire,
                           has_zero, βT)
end

# first position g in lo:hi with es[outs[g]].t ≥ x (or > x if strict), hi + 1 if none
function _first_departure(es, outs, lo, hi, x, strict::Bool)
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        t = es[outs[mid]].t
        if strict ? t <= x : t < x
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return lo
end

# Sums of ranges of a vector of non-negative values with a segment tree: the sums
# only add values, so they do not suffer from cancellation even when the numbers of
# walks span hundreds of orders of magnitude (differences of prefix sums would).
struct _SumTree{R}
    size::Int
    tree::Vector{R}
end
_SumTree{R}(n::Integer) where {R} = _SumTree{R}(nextpow(2, max(n, 1)), fill(zero(R), 2 * nextpow(2, max(n, 1))))

@inline function _set!(t::_SumTree, i::Int, x)
    j = i + t.size - 1
    @inbounds t.tree[j] = x
    j >>= 1
    @inbounds while j >= 1
        t.tree[j] = t.tree[2j] + t.tree[2j+1]
        j >>= 1
    end
    return t
end

# sum of the values at i:k
@inline function _sum(t::_SumTree{R}, i::Int, k::Int) where {R}
    s = zero(R)
    l = i + t.size - 1
    r = k + t.size
    @inbounds while l < r
        if isodd(l)
            s += t.tree[l]
            l += 1
        end
        if isodd(r)
            r -= 1
            s += t.tree[r]
        end
        l >>= 1
        r >>= 1
    end
    return s
end

# Buffers of one task for the single-source computations.
struct _WalkWork{T,C,F,R}
    C::Vector{C}           # minimum cost of the walks from s ending with each edge
    Σ::Vector{R}           # number of such walks
    P::Vector{C}           # minimum cost of the walks each edge can extend
    valid::Vector{Bool}    # C[f] = P[f] ⊕ γ(f), i.e., f has best predecessors
    mincost::Vector{C}     # unrestricted: best cost and count of the walks into each node
    mincount::Vector{R}
    inptr::Vector{Int}     # next incoming edge of each node to be considered
    head::Vector{Int}      # restless: monotone deque of the incoming edges in the window
    tail::Vector{Int}
    dq_arr::Vector{T}
    dq_cost::Vector{C}
    dq_count::_SumTree{R}  # numbers of walks of the deque entries
    cstar::Vector{F}       # optimal target cost of each node
    σstar::Vector{R}       # number of optimal walks to each node
    S::Vector{R}           # sums of b[f] / Σ[f] to the end of each block of equal P
    bkey::Vector{Tuple{C,Int}}  # restless: outgoing edges grouped by (P, position)
    bidx::Vector{Int}
    bend::Vector{Int}
    SB::_SumTree{R}        # restless: b[f] / Σ[f] by position in bkey
    bc::Vector{Float64}    # betweenness accumulated over the sources of the task
    order::Vector{Int}     # edges in the order of the forward pass
    # blocks of edges with transition time 0 at the same time stamp
    pe::Vector{C}          # best predecessor cost outside the block
    ce::Vector{R}
    ib::Vector{C}          # best cost of the block edges entering each node
    ic::Vector{R}
    bdone::Vector{Bool}
    breach::Vector{Bool}
    bfirst::Vector{Int}    # block edges leaving each node
    blast::Vector{Int}
    bmark::Vector{Int}     # first edge of the block in which bfirst/blast were set
    indeg::Vector{Int}
    queue::Vector{Int}
    bheap::MinHeap{Tuple{C,Int}}
end

function _WalkWork(wg::_WalkGraph{V,T}, dt, ::Type{R}) where {V,T,R}
    C = _bw_cost_type(dt, T)
    F = _bw_target_type(dt, T)
    n, m = wg.n, length(wg.es)
    return _WalkWork{T,C,F,R}(Vector{C}(undef, m), Vector{R}(undef, m), Vector{C}(undef, m), zeros(Bool, m),
                              Vector{C}(undef, n), Vector{R}(undef, n), zeros(Int, n),
                              zeros(Int, n), zeros(Int, n),
                              Vector{T}(undef, m), Vector{C}(undef, m), _SumTree{R}(wg.restless ? m : 0),
                              Vector{F}(undef, n), Vector{R}(undef, n), Vector{R}(undef, m),
                              Vector{Tuple{C,Int}}(undef, wg.restless ? m : 0), zeros(Int, m),
                              zeros(Int, n), _SumTree{R}(wg.restless ? m : 0), zeros(Float64, n), Int[],
                              Vector{C}(undef, m), Vector{R}(undef, m), Vector{C}(undef, n), Vector{R}(undef, n),
                              zeros(Bool, m), zeros(Bool, m), zeros(Int, n), zeros(Int, n), zeros(Int, n),
                              zeros(Int, n), Int[], MinHeap{Tuple{C,Int}}())
end

# Best predecessors of the edges leaving u at time τ: the incoming edges of u that
# arrived by τ (and not before τ - β) are moved into the deque (restless) or into the
# running minimum. With `zero_block` the edges with transition time 0 arriving at τ
# are left out (they belong to the block of edges at time τ being resolved). Returns
# the minimum cost and the number of walks with that cost.
@inline function _window!(w::_WalkWork{T,C,F,R}, wg::_WalkGraph, u::Int, τ, zero_block::Bool) where {T,C,F,R}
    es, ins, inoff = wg.es, wg.ins, wg.inoff
    inf = _bw_inf(C)
    k = w.inptr[u]
    @inbounds if wg.restless
        # the deque costs are non-decreasing from head to tail
        while k < inoff[u+1]
            e = ins[k]
            a = es[e].t + es[e].tt
            (a < τ || (a == τ && !(zero_block && iszero(es[e].tt)))) || break
            k += 1
            c = w.C[e]
            c == inf && continue
            t = w.tail[u]
            while t >= w.head[u] && w.dq_cost[t] > c
                t -= 1
            end
            t += 1
            w.dq_arr[t] = a
            w.dq_cost[t] = c
            _set!(w.dq_count, t, w.Σ[e])
            w.tail[u] = t
        end
        w.inptr[u] = k
        h = w.head[u]
        if wg.expire  # drop the edges that arrived before τ - β
            while h <= w.tail[u] && w.dq_arr[h] < τ - wg.β
                h += 1
            end
            w.head[u] = h
        end
        h > w.tail[u] && return inf, zero(R)
        p = w.dq_cost[h]
        lo, hi = h, w.tail[u]  # last entry with the minimum cost
        while lo < hi
            mid = (lo + hi + 1) >>> 1
            if w.dq_cost[mid] == p
                lo = mid
            else
                hi = mid - 1
            end
        end
        return p, _sum(w.dq_count, h, lo)
    else
        while k < inoff[u+1] && es[ins[k]].t + es[ins[k]].tt <= τ
            e = ins[k]
            k += 1
            c = w.C[e]
            if c < w.mincost[u]
                w.mincost[u] = c
                w.mincount[u] = w.Σ[e]
            elseif c == w.mincost[u] && c != inf
                w.mincount[u] += w.Σ[e]
            end
        end
        w.inptr[u] = k
        return w.mincost[u], w.mincount[u]
    end
end

# C, Σ, P and valid of edge i from its best predecessors (p, cnt); returns true iff
# the count overflowed.
@inline function _set_edge!(w::_WalkWork{T,C,F,R}, s::Int, dt, i::Int, f, p, cnt) where {T,C,F,R}
    inf = _bw_inf(C)
    cf = inf
    sf = zero(R)
    if p != inf
        cf = _bw_extend(dt, p, f)
        sf = cnt
    end
    if Int(f.u) == s
        c0 = _bw_first(dt, f)
        if c0 < cf
            cf = c0
            sf = one(R)
        elseif c0 == cf
            sf += one(R)
        end
    end
    @inbounds begin
        w.C[i] = cf
        w.Σ[i] = sf
        w.P[i] = p
        w.valid[i] = p != inf && _bw_extend(dt, p, f) == cf
    end
    return !isfinite(sf)
end

# criteria whose costs grow with every edge: optimal walks never use a cycle of edges
# at the same time stamp
const _HopCriterion = Union{MinimumHops,ShortestForemost,ShortestLatest,ShortestFastest}

# best predecessor (cost, count) of block edge x: outside the block or in the block
@inline function _block_pred(w::_WalkWork{T,C,F,R}, x::Int, u::Int) where {T,C,F,R}
    @inbounds pe, ib = w.pe[x], w.ib[u]
    p = min(pe, ib)
    p == _bw_inf(C) && return p, zero(R)
    return p, (pe == p ? w.ce[x] : zero(R)) + (ib == p ? w.ic[u] : zero(R))
end

@inline function _block_update!(w::_WalkWork{T,C,F,R}, v::Int, c, σ) where {T,C,F,R}
    @inbounds if c < w.ib[v]
        w.ib[v] = c
        w.ic[v] = σ
    elseif c == w.ib[v] && c != _bw_inf(C)
        w.ic[v] += σ
    end
    return nothing
end

# Edges i:j have transition time 0 and the same time stamp τ, so they can follow each
# other. Their best predecessors before τ come from the windows; inside the block they
# are resolved by increasing cost (criteria whose costs grow with each edge) or in a
# topological order of the block (other criteria, for which a reachable cycle would
# give infinitely many optimal walks).
function _zero_block!(w::_WalkWork{T,C,F,R}, wg::_WalkGraph, s::Int, dt, i::Int, j::Int) where {T,C,F,R}
    es = wg.es
    inf = _bw_inf(C)
    τ = es[i].t
    overflow = false
    k = i
    @inbounds while k <= j  # the block is sorted by tail
        u = Int(es[k].u)
        l = k
        while l < j && es[l+1].u == u
            l += 1
        end
        p, cnt = _window!(w, wg, u, τ, true)
        for x in k:l
            w.pe[x] = p
            w.ce[x] = cnt
        end
        w.bfirst[u], w.blast[u], w.bmark[u] = k, l, i
        k = l + 1
    end
    @inbounds for x in i:j  # no block edge has been resolved yet
        for v in (Int(es[x].u), Int(es[x].v))
            w.ib[v] = inf
            w.ic[v] = zero(R)
        end
        w.bdone[x] = false
    end
    block_out(v) = @inbounds w.bmark[v] == i ? (w.bfirst[v]:w.blast[v]) : (1:0)
    function tentative(x)
        f = @inbounds es[x]
        p, _ = _block_pred(w, x, Int(f.u))
        c = p == inf ? inf : _bw_extend(dt, p, f)
        Int(f.u) == s && (c = min(c, _bw_first(dt, f)))
        return c
    end
    function finalize!(x)
        f = @inbounds es[x]
        u = Int(f.u)
        p, cnt = _block_pred(w, x, u)
        overflow |= _set_edge!(w, s, dt, x, f, p, cnt)
        @inbounds w.bdone[x] = true
        push!(w.order, x)
        @inbounds w.C[x] != inf && _block_update!(w, Int(f.v), w.C[x], w.Σ[x])
        return nothing
    end
    if dt isa _HopCriterion
        heap = empty!(w.bheap)
        @inbounds for x in i:j
            c = tentative(x)
            w.C[x] = c
            c != inf && push!(heap, (c, x))
        end
        @inbounds while !isempty(heap)
            c, x = pop!(heap)
            (w.bdone[x] || c != w.C[x]) && continue
            finalize!(x)
            for y in block_out(Int(es[x].v))
                w.bdone[y] && continue
                cy = tentative(y)
                if cy < w.C[y]
                    w.C[y] = cy
                    push!(heap, (cy, y))
                end
            end
        end
    else
        # block edges reachable from the source or from the edges before τ
        queue = empty!(w.queue)
        @inbounds for x in i:j
            w.breach[x] = w.pe[x] != inf || Int(es[x].u) == s
            w.breach[x] && push!(queue, x)
        end
        h = 1
        @inbounds while h <= length(queue)
            x = queue[h]
            h += 1
            for y in block_out(Int(es[x].v))
                w.breach[y] || (w.breach[y] = true; push!(queue, y))
            end
        end
        # Kahn's algorithm on the nodes, counting the reachable block edges entering them
        @inbounds for x in i:j
            w.indeg[es[x].v] = 0
            w.indeg[es[x].u] = 0
        end
        @inbounds for x in i:j
            w.breach[x] && (w.indeg[es[x].v] += 1)
        end
        empty!(queue)
        @inbounds for x in i:j
            u = Int(es[x].u)
            (w.breach[x] && w.indeg[u] == 0 && first(block_out(u)) == x) && push!(queue, u)
        end
        h = 1
        @inbounds while h <= length(queue)
            u = queue[h]
            h += 1
            for x in block_out(u)
                w.breach[x] || continue
                finalize!(x)
                v = Int(es[x].v)
                w.indeg[v] -= 1
                (w.indeg[v] == 0 && !isempty(block_out(v))) && push!(queue, v)
            end
        end
        @inbounds for x in i:j
            (w.breach[x] && !w.bdone[x]) && throw(ArgumentError(
                "edges with transition time 0 form a cycle at time $τ: there are infinitely many optimal " *
                "$(nameof(typeof(dt))) walks; use positive transition times (strict walks)"))
        end
    end
    @inbounds for x in i:j  # unreachable block edges
        w.bdone[x] || finalize!(x)
    end
    return overflow
end

# Forward pass from s: C, Σ, P and valid for every edge, and the processing order
# of the edges. Returns true iff a count is not finite (overflow of a floating point
# count type).
function _walk_forward!(w::_WalkWork{T,C,F,R}, wg::_WalkGraph{V,T}, s::Int, dt) where {T,C,F,R,V}
    es, inoff = wg.es, wg.inoff
    inf = _bw_inf(C)
    overflow = false
    @inbounds for u in 1:wg.n
        w.inptr[u] = inoff[u]
        w.mincost[u] = inf
        w.mincount[u] = zero(R)
        w.head[u] = inoff[u]
        w.tail[u] = inoff[u] - 1
        w.bmark[u] = 0
    end
    empty!(w.order)
    m = length(es)
    i = 1
    @inbounds while i <= m
        f = es[i]
        if wg.has_zero && iszero(f.tt)
            j = i
            while j < m && es[j+1].t == f.t && iszero(es[j+1].tt)
                j += 1
            end
            overflow |= _zero_block!(w, wg, s, dt, i, j)
            i = j + 1
            continue
        end
        p, cnt = _window!(w, wg, Int(f.u), f.t, false)
        overflow |= _set_edge!(w, s, dt, i, f, p, cnt)
        push!(w.order, i)
        i += 1
    end
    return overflow
end

# Target pass: optimal target cost and number of optimal walks for every node ≠ s.
function _walk_targets!(w::_WalkWork{T,C,F,R}, wg::_WalkGraph, s::Int, dt, target::Int=0) where {T,C,F,R}
    fill!(w.cstar, _bw_inf(F))
    fill!(w.σstar, zero(R))
    inf = _bw_inf(C)
    @inbounds for i in eachindex(wg.es)
        e = wg.es[i]
        v = Int(e.v)
        (v == s || w.C[i] == inf || (target != 0 && v != target)) && continue
        tc = _bw_target(dt, w.C[i], e)
        if tc < w.cstar[v]
            w.cstar[v] = tc
            w.σstar[v] = w.Σ[i]
        elseif tc == w.cstar[v]
            w.σstar[v] += w.Σ[i]
        end
    end
    return w
end

# Group the valid outgoing edges of every node by (P, position) for the restless
# backward pass.
function _walk_buckets!(w::_WalkWork, wg::_WalkGraph)
    @inbounds for v in 1:wg.n
        lo = wg.outoff[v]
        k = lo
        for gpos in lo:wg.outoff[v+1]-1
            f = wg.outs[gpos]
            if w.valid[f]
                w.bkey[k] = (w.P[f], gpos)
                k += 1
            end
            w.bidx[gpos] = 0
        end
        w.bend[v] = k - 1
        k > lo + 1 && sort!(view(w.bkey, lo:k-1))
        for j in lo:k-1
            w.bidx[w.bkey[j][2]] = j
        end
    end
    return w
end

# Backward pass: s-betweenness of every edge, accumulated on the heads (edges into
# `target` only end walks if target ≠ 0). Calls `f(v, value)` for every edge.
function _walk_backward!(onedge::G, w::_WalkWork{T,C,F,R}, wg::_WalkGraph, s::Int, dt) where {G,T,C,F,R}
    es, outs, outoff = wg.es, wg.outs, wg.outoff
    inf = _bw_inf(C)
    wg.restless && _walk_buckets!(w, wg)
    @inbounds for k in length(w.order):-1:1
        i = w.order[k]
        e = es[i]
        u, v = Int(e.u), Int(e.v)
        ce = w.C[i]
        bi = zero(R)
        if ce != inf
            g1 = wg.first_succ[i]
            if wg.restless
                g3 = wg.last_succ[i]
                if g1 <= g3
                    seg = view(w.bkey, outoff[v]:w.bend[v])
                    lo = outoff[v] - 1 + searchsortedfirst(seg, (ce, g1))
                    hi = outoff[v] - 1 + searchsortedlast(seg, (ce, g3))
                    lo <= hi && (bi = w.Σ[i] * _sum(w.SB, lo, hi))
                end
            else
                # P is non-increasing and at most ce from g1 on, so the successors are
                # the block of equal P starting at g1 if that P equals ce
                if g1 < outoff[v+1] && w.P[outs[g1]] == ce
                    bi = w.Σ[i] * w.S[g1]
                end
            end
            if v != s && w.σstar[v] > zero(R) && _bw_target(dt, ce, e) == w.cstar[v]
                bi += w.Σ[i] / w.σstar[v]
            end
        end
        v != s && onedge(v, bi)
        # e as a successor candidate of the edges entering u
        gpos = wg.outpos[i]
        val = (w.valid[i] && w.Σ[i] > zero(R)) ? bi / w.Σ[i] : zero(R)
        if wg.restless
            k = w.bidx[gpos]
            k != 0 && _set!(w.SB, k, val)
        else
            same = gpos < outoff[u+1] - 1 && w.P[outs[gpos+1]] == w.P[i]
            w.S[gpos] = val + (same ? w.S[gpos+1] : zero(R))
        end
    end
    return w
end

# Betweenness of all nodes and whether a walk count overflowed.
function _walk_betweenness(wg::_WalkGraph, dt, ::Type{R}; parallel::Bool=true) where {R}
    n = wg.n
    n == 0 && return Float64[], false
    function accumulate(sources)
        w = _WalkWork(wg, dt, R)
        bc = w.bc
        overflow = false
        for s in sources
            overflow |= _walk_forward!(w, wg, s, dt)
            _walk_targets!(w, wg, s, dt)
            _walk_backward!(w, wg, s, dt) do v, x
                @inbounds bc[v] += Float64(x)
            end
            @inbounds for v in 1:n
                (v != s && w.σstar[v] > zero(R)) && (bc[v] -= 1)
            end
        end
        return bc, overflow
    end
    parallel || return accumulate(1:n)
    return tmapreduce(accumulate, (x, y) -> (x[1] .+ y[1], x[2] | y[2]), chunks(1:n; n=min(n, 4 * Threads.nthreads())))
end

# Runs `f(R)` with the count type R and, if the walk counts overflow Float64, again
# with BigFloat.
function _with_count_type(f, ::Type{R}, explicit::Bool) where {R}
    result, overflow = _unwrap_task_errors(() -> f(R))
    overflow || return result
    (explicit || R !== Float64) && throw(OverflowError(
        "the numbers of optimal walks overflow $R; use count_type = BigFloat or Rational{BigInt}"))
    @warn "The numbers of optimal walks overflow Float64, recomputing with BigFloat (slower). " *
          "Pass count_type = BigFloat to skip the first attempt."
    return first(_unwrap_task_errors(() -> f(BigFloat)))
end

"""
    temporal_betweenness(g::OrderedEdgeList, criterion = MinimumHops(), ti = time_interval(g);
                         β = Inf, count_type = nothing)

Temporal betweenness of all nodes: for every node `u`, the sum over all pairs of
nodes `s ≠ t` different from `u`, with `t` reachable from `s`, of the fraction of the
optimal temporal walks from `s` to `t` that pass through `u` (walks passing several
times through `u` are counted with multiplicity).

The optimality `criterion` is one of

| criterion | optimal walks |
|---|---|
| [`MinimumHops`](@ref) | fewest edges (*shortest*) |
| [`EarliestArrival`](@ref) | earliest arrival (*foremost*) |
| [`LatestDeparture`](@ref) | latest departure (*latest*) |
| [`Fastest`](@ref) | minimum duration |
| [`ShortestForemost`](@ref), [`ShortestLatest`](@ref), [`ShortestFastest`](@ref) | the previous three, ties broken by fewest edges |
| [`PrefixForemost`](@ref) | prefix foremost paths (Buß et al., 2020) |

`β` is the maximum waiting time at a node (restless walks); `β = Inf` means no
constraint. Only edges inside the time window `ti` are used.

Positive transition times give *strict* walks. Edges with transition time 0 give
*non-strict* walks (Zhang et al., WWW 2024): such an edge can be followed at once by
another edge with the same time stamp. If the edges with transition time 0 and the
same time stamp form a cycle that a walk can reach, there are infinitely many optimal
walks for every criterion except the shortest ones (`MinimumHops`, `ShortestForemost`,
`ShortestLatest`, `ShortestFastest`), and an `ArgumentError` is thrown. Prefix
foremost betweenness requires positive transition times.

By default (`count_type = nothing`) the numbers of optimal walks are counted in `Float64`; if they
overflow (typical for foremost and fastest walks on large graphs, whose optimal
walks are very numerous) the computation is repeated with `BigFloat`. Pass
`count_type = BigFloat` to use it directly, or `count_type = Rational{BigInt}` to
count the walks and compute the contribution of every source exactly. The result is
a `Vector{Float64}` in every case: the contributions are rounded to `Float64` and
summed.

With `MinimumHops()` and no waiting constraint this is the *shortest temporal
betweenness* of Buß et al. (KDD 2020), since shortest strict walks are paths. The
algorithm follows the approach of Brunelli, Crescenzi and Viennot (KDD 2024) and runs
in `O(n M)` time without waiting constraint and `O(n M log M)` with it (`M` edges);
sources are processed in parallel.

# References

- F. Brunelli, P. Crescenzi, and L. Viennot. *Making temporal betweenness computation faster and restless.* KDD, 2024. [DOI](https://doi.org/10.1145/3637528.3671825), [arXiv](https://arxiv.org/abs/2501.12708)
- S. Buß, H. Molter, R. Niedermeier, and M. Rymar. *Algorithmic aspects of temporal betweenness.* KDD, 2020. [DOI](https://doi.org/10.1145/3394486.3403259), [arXiv](https://arxiv.org/abs/2006.08668)
- T. Zhang, Y. Gao, J. Zhao, L. Chen, L. Jin, Z. Yang, B. Cao, and J. Fan. *Efficient exact and approximate betweenness centrality computation for temporal graphs.* The Web Conference (WWW), 2024. [DOI](https://doi.org/10.1145/3589334.3645438)
"""
function temporal_betweenness(g::OrderedEdgeList, criterion::DistanceType=MinimumHops(), ti=time_interval(g);
                              β::Real=Inf, count_type::Union{Nothing,Type{<:Real}}=nothing)
    R = count_type === nothing ? Float64 : count_type
    if criterion isa PrefixForemost
        isfinite(β) && throw(ArgumentError("prefix foremost betweenness has no waiting constraint"))
        return _with_count_type(R -> _prefix_foremost_betweenness(g, ti, R), R, count_type !== nothing)
    end
    criterion isa _WalkCriterion ||
        throw(ArgumentError("temporal betweenness is not defined for $(typeof(criterion))"))
    wg = _walk_graph(g, ti, β)
    return _with_count_type(R -> _walk_betweenness(wg, criterion, R), R, count_type !== nothing)
end
