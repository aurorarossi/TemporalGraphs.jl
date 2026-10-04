# δ-temporal motifs, following Paranjape, Benson and Leskovec, "Motifs in temporal
# networks" (WSDM 2017). A k-node, l-edge δ-temporal motif is a sequence of l directed
# edges on k nodes whose static graph is connected; an instance is a sequence of l
# temporal edges of the graph, ordered in time and spanning at most δ time units,
# that maps onto the motif by a bijection of the nodes.
#
# Equal time stamps: with `strict = true` the time stamps of an instance must be
# strictly increasing (the definition of the paper); with `strict = false` edges
# with the same time stamp are ordered as in the edge list of the graph. All counts
# are computed on "blocks" of edges: maximal runs of equal time stamps (strict) or
# single edges (non-strict). A block enters and leaves the sliding windows at once,
# so that no instance contains two edges of the same block.
#
# Only the departure times `t` are used (transition times are ignored) and self
# loops are skipped.

# Last index of the block starting at i.
@inline function _block_last(ts::AbstractVector, i::Int, L::Int, strict::Bool)
    strict || return i
    j = i
    @inbounds while j < L && ts[j+1] == ts[i]
        j += 1
    end
    return j
end

# Generic sweep of Algorithm 2 of the paper: for every block, the "pre" window holds
# the earlier edges within δ, the "post" window the later edges within δ.
function _motif_sweep!(c, ts::AbstractVector, L::Int, δ, strict::Bool)
    start = 1
    stop = 1
    i = 1
    @inbounds while i <= L
        j = _block_last(ts, i, L, strict)
        τ = ts[i]
        while start < i && ts[start] + δ < τ
            s2 = _block_last(ts, start, L, strict)
            _pop_block!(c, true, start, s2)
            start = s2 + 1
        end
        while stop <= L && ts[stop] <= τ + δ
            s2 = _block_last(ts, stop, L, strict)
            _push_block!(c, false, stop, s2)
            stop = s2 + 1
        end
        _pop_block!(c, false, i, j)  # the current block is the earliest of the post window
        _process_block!(c, i, j)
        _push_block!(c, true, i, j)
        i = j + 1
    end
    return c
end

# ----------------------------------------------------------------------------------
# Grid of the 36 motifs with 3 edges on 2 or 3 nodes (Figure 3 of the paper, and the
# output of the SNAP implementation): with first edge a → b and third node c, the
# row is given by the second edge and the column by the third edge, in the order
# a → b, b → a, a → c, c → a, b → c, c → b (reversed for the rows).

const _MOTIF_EDGES = ((1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2))

# Grid cell of the motif formed by three edges on abstract nodes.
function _motif_cell(e1::Tuple{Int,Int}, e2::Tuple{Int,Int}, e3::Tuple{Int,Int})
    a, b = e1
    lab(x) = x == a ? 1 : x == b ? 2 : 3
    pos(e) = findfirst(==((lab(e[1]), lab(e[2]))), _MOTIF_EDGES)
    return (7 - pos(e2), pos(e3))
end

# ----------------------------------------------------------------------------------
# Stars and 2-node motifs: one sweep over the edges incident to each center node u.
# An incident edge has a neighbor k (local index) and a direction (1 out of u, 2 into
# u). The three classes of 3-node stars are pre (1st and 2nd edge with the same
# neighbor), post (2nd and 3rd) and mid (1st and 3rd). The counters keep the pairs of
# edges with the same neighbor in the windows, in total and per neighbor: the
# per-neighbor pairs give the 2-node motifs, which are subtracted from the stars.

mutable struct _StarCounter
    K::Vector{Int}          # local neighbor of each incident edge
    D::Vector{Int8}         # direction of each incident edge
    two::Vector{Bool}       # count 2-node motifs for this neighbor (u < neighbor)
    pre_nodes::Vector{Int}  # [d, k]: edges of the window with direction d and neighbor k
    post_nodes::Vector{Int}
    pre_sum::Vector{Int}    # [d1, d2]: pairs of the window with the same neighbor
    post_sum::Vector{Int}
    mid_sum::Vector{Int}    # pairs (before, after the current block) within δ
    pre_pair::Vector{Int}   # [d1, d2, k]: the same, for neighbor k
    post_pair::Vector{Int}
    mid_pair::Vector{Int}
    cnt_pre::Vector{Int}    # [d1, d2, d3]: stars of class pre with these directions
    cnt_post::Vector{Int}
    cnt_mid::Vector{Int}
    cnt_two::Vector{Int}    # 2-node motifs
end

_StarCounter() = _StarCounter(Int[], Int8[], Bool[], Int[], Int[], zeros(Int, 4), zeros(Int, 4), zeros(Int, 4),
                              Int[], Int[], Int[], zeros(Int, 8), zeros(Int, 8), zeros(Int, 8), zeros(Int, 8))

@inline _i2(d1, d2) = 2 * (d1 - 1) + d2
@inline _i3(d1, d2, d3) = 4 * (d1 - 1) + 2 * (d2 - 1) + d3
@inline _ik(d, k) = 2 * (k - 1) + d
@inline _ip(d1, d2, k) = 4 * (k - 1) + 2 * (d1 - 1) + d2

function _reset!(c::_StarCounter, L::Int, deg::Int)
    resize!(c.K, L)
    resize!(c.D, L)
    resize!(c.two, deg)
    for v in (c.pre_nodes, c.post_nodes)
        resize!(v, 2deg)
        fill!(v, 0)
    end
    for v in (c.pre_pair, c.post_pair, c.mid_pair)
        resize!(v, 4deg)
        fill!(v, 0)
    end
    fill!(c.pre_sum, 0)
    fill!(c.post_sum, 0)
    fill!(c.mid_sum, 0)
    return c
end

function _push_block!(c::_StarCounter, pre::Bool, lo::Int, hi::Int)
    nodes, sum, pair = pre ? (c.pre_nodes, c.pre_sum, c.pre_pair) : (c.post_nodes, c.post_sum, c.post_pair)
    @inbounds for e in lo:hi
        k, d = c.K[e], c.D[e]
        for d1 in 1:2
            x = nodes[_ik(d1, k)]
            sum[_i2(d1, d)] += x
            pair[_ip(d1, d, k)] += x
        end
    end
    @inbounds for e in lo:hi
        nodes[_ik(c.D[e], c.K[e])] += 1
    end
    return c
end

function _pop_block!(c::_StarCounter, pre::Bool, lo::Int, hi::Int)
    nodes, sum, pair = pre ? (c.pre_nodes, c.pre_sum, c.pre_pair) : (c.post_nodes, c.post_sum, c.post_pair)
    @inbounds for e in lo:hi
        nodes[_ik(c.D[e], c.K[e])] -= 1
    end
    @inbounds for e in lo:hi
        k, d = c.K[e], c.D[e]
        for d2 in 1:2
            x = nodes[_ik(d2, k)]
            sum[_i2(d, d2)] -= x
            pair[_ip(d, d2, k)] -= x
        end
    end
    return c
end

function _process_block!(c::_StarCounter, lo::Int, hi::Int)
    @inbounds for e in lo:hi  # pairs ending in the block no longer straddle it
        k, d = c.K[e], c.D[e]
        for d1 in 1:2
            x = c.pre_nodes[_ik(d1, k)]
            c.mid_sum[_i2(d1, d)] -= x
            c.mid_pair[_ip(d1, d, k)] -= x
        end
    end
    @inbounds for e in lo:hi
        k, d = c.K[e], c.D[e]
        for d1 in 1:2, d2 in 1:2
            p = c.pre_pair[_ip(d1, d2, k)]
            c.cnt_pre[_i3(d1, d2, d)] += c.pre_sum[_i2(d1, d2)] - p
            c.cnt_post[_i3(d, d1, d2)] += c.post_sum[_i2(d1, d2)] - c.post_pair[_ip(d1, d2, k)]
            c.cnt_mid[_i3(d1, d, d2)] += c.mid_sum[_i2(d1, d2)] - c.mid_pair[_ip(d1, d2, k)]
            c.two[k] && (c.cnt_two[_i3(d1, d2, d)] += p)
        end
    end
    @inbounds for e in lo:hi  # pairs starting in the block
        k, d = c.K[e], c.D[e]
        for d2 in 1:2
            x = c.post_nodes[_ik(d2, k)]
            c.mid_sum[_i2(d, d2)] += x
            c.mid_pair[_ip(d, d2, k)] += x
        end
    end
    return c
end

# ----------------------------------------------------------------------------------
# Triangles: one sweep per pair of nodes (a, b) over the edges between a and b and
# the edges between a or b and the third nodes w of the triangles assigned to the
# pair. A side edge has a side s (1 if it is incident to a, 2 to b), a direction
# (1 from a or b to w, 2 towards a or b) and a third node k (local index). The edges
# between a and b have direction 1 (a → b) or 2. The counters keep the pairs of side
# edges in the windows with the same third node and different sides.

mutable struct _TriangleCounter
    S::Vector{Int8}         # side, 0 for an edge between a and b
    D::Vector{Int8}
    K::Vector{Int}
    pre_nodes::Vector{Int}  # [s, d, k]
    post_nodes::Vector{Int}
    pre_sum::Vector{Int}    # [s, d1, d2]: pairs whose first edge has side s
    post_sum::Vector{Int}
    mid_sum::Vector{Int}
    cnt_pre::Vector{Int}    # [d_ab, s, d1, d2]: the edge between a and b is the 3rd
    cnt_post::Vector{Int}   # ... the 1st
    cnt_mid::Vector{Int}    # ... the 2nd
end

_TriangleCounter() = _TriangleCounter(Int8[], Int8[], Int[], Int[], Int[], zeros(Int, 8), zeros(Int, 8),
                                      zeros(Int, 8), zeros(Int, 16), zeros(Int, 16), zeros(Int, 16))

@inline _ikt(s, d, k) = 4 * (k - 1) + 2 * (s - 1) + d
@inline _i4(d, s, d1, d2) = 8 * (d - 1) + 4 * (s - 1) + 2 * (d1 - 1) + d2

function _reset!(c::_TriangleCounter, L::Int, deg::Int)
    resize!(c.S, L)
    resize!(c.D, L)
    resize!(c.K, L)
    for v in (c.pre_nodes, c.post_nodes)
        resize!(v, 4deg)
        fill!(v, 0)
    end
    fill!(c.pre_sum, 0)
    fill!(c.post_sum, 0)
    fill!(c.mid_sum, 0)
    return c
end

function _push_block!(c::_TriangleCounter, pre::Bool, lo::Int, hi::Int)
    nodes, sum = pre ? (c.pre_nodes, c.pre_sum) : (c.post_nodes, c.post_sum)
    @inbounds for e in lo:hi
        s = c.S[e]
        s == 0 && continue
        k, d, o = c.K[e], c.D[e], 3 - s
        for d1 in 1:2
            sum[_i3(o, d1, d)] += nodes[_ikt(o, d1, k)]
        end
    end
    @inbounds for e in lo:hi
        c.S[e] == 0 || (nodes[_ikt(c.S[e], c.D[e], c.K[e])] += 1)
    end
    return c
end

function _pop_block!(c::_TriangleCounter, pre::Bool, lo::Int, hi::Int)
    nodes, sum = pre ? (c.pre_nodes, c.pre_sum) : (c.post_nodes, c.post_sum)
    @inbounds for e in lo:hi
        c.S[e] == 0 || (nodes[_ikt(c.S[e], c.D[e], c.K[e])] -= 1)
    end
    @inbounds for e in lo:hi
        s = c.S[e]
        s == 0 && continue
        k, d, o = c.K[e], c.D[e], 3 - s
        for d2 in 1:2
            sum[_i3(s, d, d2)] -= nodes[_ikt(o, d2, k)]
        end
    end
    return c
end

function _process_block!(c::_TriangleCounter, lo::Int, hi::Int)
    @inbounds for e in lo:hi
        s = c.S[e]
        s == 0 && continue
        k, d, o = c.K[e], c.D[e], 3 - s
        for d1 in 1:2
            c.mid_sum[_i3(o, d1, d)] -= c.pre_nodes[_ikt(o, d1, k)]
        end
    end
    @inbounds for e in lo:hi
        c.S[e] == 0 || continue
        d = c.D[e]
        for s in 1:2, d1 in 1:2, d2 in 1:2
            x = _i3(s, d1, d2)
            y = _i4(d, s, d1, d2)
            c.cnt_pre[y] += c.pre_sum[x]
            c.cnt_post[y] += c.post_sum[x]
            c.cnt_mid[y] += c.mid_sum[x]
        end
    end
    @inbounds for e in lo:hi
        s = c.S[e]
        s == 0 && continue
        k, d, o = c.K[e], c.D[e], 3 - s
        for d2 in 1:2
            c.mid_sum[_i3(s, d, d2)] += c.post_nodes[_ikt(o, d2, k)]
        end
    end
    return c
end

# ----------------------------------------------------------------------------------
# Static structure.

# Incident edges of every node in time order, without self loops: (edge index,
# neighbor, direction 1 out / 2 in).
function _incidence(g::OrderedEdgeList)
    n = num_nodes(g)
    es = edges(g)
    off = zeros(Int, n + 1)
    for e in es
        e.u == e.v && continue
        off[e.u+1] += 1
        off[e.v+1] += 1
    end
    off[1] = 1
    for u in 1:n
        off[u+1] += off[u]
    end
    pos = off[1:n]
    idx = Vector{Int}(undef, off[n+1] - 1)
    nbr = Vector{Int}(undef, length(idx))
    dir = Vector{Int8}(undef, length(idx))
    for (i, e) in enumerate(es)
        u, v = Int(e.u), Int(e.v)
        u == v && continue
        idx[pos[u]], nbr[pos[u]], dir[pos[u]] = i, v, 1
        pos[u] += 1
        idx[pos[v]], nbr[pos[v]], dir[pos[v]] = i, u, 2
        pos[v] += 1
    end
    return off, idx, nbr, dir
end

# Undirected static graph: the pairs a < b joined by some edge, the indices of their
# edges in time order, and sorted adjacency lists with the pair ids.
function _static_pairs(g::OrderedEdgeList)
    n = num_nodes(g)
    es = edges(g)
    keys = Tuple{Int,Int}[]
    ids = Int[]
    for (i, e) in enumerate(es)
        e.u == e.v && continue
        push!(keys, minmax(Int(e.u), Int(e.v)))
        push!(ids, i)
    end
    p = sortperm(keys; alg=Base.Sort.DEFAULT_STABLE)
    pairs = Tuple{Int,Int}[]
    pair_off = Int[1]
    edge_pair = zeros(Int, length(es))  # 0 for self loops
    for (r, j) in enumerate(p)
        if r == 1 || keys[j] != keys[p[r-1]]
            r > 1 && push!(pair_off, r)
            push!(pairs, keys[j])
        end
        edge_pair[ids[j]] = length(pairs)
    end
    push!(pair_off, length(p) + 1)
    isempty(pairs) && (pair_off = [1])
    adj_off = zeros(Int, n + 1)
    for (a, b) in pairs
        adj_off[a+1] += 1
        adj_off[b+1] += 1
    end
    adj_off[1] = 1
    for u in 1:n
        adj_off[u+1] += adj_off[u]
    end
    pos = adj_off[1:n]
    adj = Vector{Int}(undef, 2 * length(pairs))
    adj_pair = Vector{Int}(undef, 2 * length(pairs))
    for (q, (a, b)) in enumerate(pairs)  # pairs are sorted, so are the lists
        adj[pos[a]], adj_pair[pos[a]] = b, q
        pos[a] += 1
        adj[pos[b]], adj_pair[pos[b]] = a, q
        pos[b] += 1
    end
    for u in 1:n
        r = adj_off[u]:adj_off[u+1]-1
        o = sortperm(view(adj, r))
        adj[r] = adj[r][o]
        adj_pair[r] = adj_pair[r][o]
    end
    return pairs, pair_off, edge_pair, adj_off, adj, adj_pair
end

# Triangles of the static graph, each assigned to its pair (a, b) with the most
# temporal edges (Algorithm 5 of the paper), grouped by pair: group G has pair
# gpair[G] and the third nodes x in goff[G]:goff[G+1]-1, with the pairs qa[x] = (a, w)
# and qb[x] = (b, w).
function _triangle_groups(n, pairs, pair_off, adj_off, adj, adj_pair)
    # orientation from lower to higher (degree, id): every node has O(√M) higher neighbors
    higher(u, v) = (adj_off[u+1] - adj_off[u], u) < (adj_off[v+1] - adj_off[v], v)
    hoff = ones(Int, n + 1)
    for u in 1:n
        hoff[u+1] = hoff[u] + count(i -> higher(u, adj[i]), adj_off[u]:adj_off[u+1]-1)
    end
    hadj = Vector{Int}(undef, hoff[n+1] - 1)
    hpair = similar(hadj)
    for u in 1:n
        x = hoff[u]
        for i in adj_off[u]:adj_off[u+1]-1
            if higher(u, adj[i])
                hadj[x], hpair[x] = adj[i], adj_pair[i]
                x += 1
            end
        end
    end
    σ(q) = pair_off[q+1] - pair_off[q]
    # (q, w, pair of w with min(q), pair of w with max(q))
    assign(q, w, qu, qv, u, v) = u < v ? (q, w, qu, qv) : (q, w, qv, qu)
    parts = tmap(chunks(1:n; n=min(max(n, 1), 8 * Threads.nthreads()))) do us
        mark = zeros(Int, n)  # pair id of (u, w) for the higher neighbors w of u
        out = NTuple{4,Int}[]
        for u in us
            for i in hoff[u]:hoff[u+1]-1
                mark[hadj[i]] = hpair[i]
            end
            for i in hoff[u]:hoff[u+1]-1
                v, quv = hadj[i], hpair[i]
                for j in hoff[v]:hoff[v+1]-1
                    w = hadj[j]
                    quw = mark[w]
                    quw == 0 && continue
                    qvw = hpair[j]
                    # the pair with the most edges, ties broken by pair id
                    best = max((σ(quv), quv), (σ(qvw), qvw), (σ(quw), quw))[2]
                    push!(out, best == quv ? assign(quv, w, quw, qvw, u, v) :
                               best == qvw ? assign(qvw, u, quv, quw, v, w) :
                                             assign(quw, v, quv, qvw, u, w))
                end
            end
            for i in hoff[u]:hoff[u+1]-1
                mark[hadj[i]] = 0
            end
        end
        out
    end
    # counting sort by pair
    P = length(pairs)
    cnt = zeros(Int, P + 1)
    for part in parts, t in part
        cnt[t[1]+1] += 1
    end
    gpair = Int[]
    goff = Int[1]
    pos = zeros(Int, P)
    for q in 1:P
        if cnt[q+1] > 0
            push!(gpair, q)
            pos[q] = goff[end]
            push!(goff, goff[end] + cnt[q+1])
        end
    end
    τ = goff[end] - 1
    W = Vector{Int}(undef, τ)
    qa = Vector{Int}(undef, τ)
    qb = Vector{Int}(undef, τ)
    for part in parts, (q, w, x, y) in part
        W[pos[q]], qa[pos[q]], qb[pos[q]] = w, x, y
        pos[q] += 1
    end
    return gpair, goff, W, qa, qb
end

# The edges of every group in time order, computed with one pass over the edges: an
# edge of pair p is appended to the groups listed for p. Entries are the edge index
# and a code: 0 for an edge between a and b, 2k + s - 1 for a side edge to the k-th
# third node of the group, on side s.
function _triangle_streams(edge_pair, pair_off, gpair, goff, qa, qb)
    P = length(pair_off) - 1
    σ(q) = pair_off[q+1] - pair_off[q]
    loff = zeros(Int, P + 1)
    for G in eachindex(gpair)
        loff[gpair[G]+1] += 1
        for x in goff[G]:goff[G+1]-1
            loff[qa[x]+1] += 1
            loff[qb[x]+1] += 1
        end
    end
    loff[1] = 1
    for q in 1:P
        loff[q+1] += loff[q]
    end
    lpos = loff[1:P]
    lgroup = Vector{Int32}(undef, loff[P+1] - 1)
    lcode = Vector{Int32}(undef, length(lgroup))
    soff = ones(Int, length(gpair) + 1)
    for G in eachindex(gpair)
        q = gpair[G]
        lgroup[lpos[q]], lcode[lpos[q]] = G, 0
        lpos[q] += 1
        size = σ(q)
        for (k, x) in enumerate(goff[G]:goff[G+1]-1)
            for (s, p) in ((1, qa[x]), (2, qb[x]))
                lgroup[lpos[p]], lcode[lpos[p]] = G, 2k + s - 1
                lpos[p] += 1
                size += σ(p)
            end
        end
        soff[G+1] = soff[G] + size
    end
    sidx = Vector{Int32}(undef, soff[end] - 1)
    scode = Vector{Int32}(undef, length(sidx))
    spos = soff[1:end-1]
    @inbounds for (i, p) in enumerate(edge_pair)
        p == 0 && continue
        for r in loff[p]:loff[p+1]-1
            G = lgroup[r]
            sidx[spos[G]], scode[spos[G]] = i, lcode[r]
            spos[G] += 1
        end
    end
    return soff, sidx, scode
end

function _triangle_motifs!(c::_TriangleCounter, ts::Vector, es, a::Int, b::Int, sidx, scode,
                           r::UnitRange{Int}, nw::Int, δ, strict::Bool)
    L = length(r)
    _reset!(c, L, nw)
    resize!(ts, L)
    @inbounds for (e, x) in enumerate(r)
        edge = es[sidx[x]]
        code = Int(scode[x])
        if code == 0
            c.S[e], c.D[e], c.K[e] = 0, edge.u == a ? 1 : 2, 0
        else
            k, s1 = divrem(code, 2)
            c.S[e], c.D[e], c.K[e] = s1 + 1, edge.u == (s1 == 0 ? a : b) ? 1 : 2, k
        end
        ts[e] = edge.t
    end
    _motif_sweep!(c, ts, L, δ, strict)
    return c
end

# ----------------------------------------------------------------------------------

function _star_motifs!(c::_StarCounter, ts::Vector, loc::Vector{Int}, off, idx, nbr, dir, es, u::Int, δ, strict::Bool)
    r = off[u]:off[u+1]-1
    L = length(r)
    L < 3 && return c
    deg = 0
    for x in r
        loc[nbr[x]] == 0 && (deg += 1; loc[nbr[x]] = deg)
    end
    _reset!(c, L, deg)
    resize!(ts, L)
    @inbounds for (e, x) in enumerate(r)
        k = loc[nbr[x]]
        c.K[e] = k
        c.D[e] = dir[x]
        c.two[k] = u < nbr[x]
        ts[e] = es[idx[x]].t
    end
    _motif_sweep!(c, ts, L, δ, strict)
    for x in r
        loc[nbr[x]] = 0
    end
    return c
end

_motif_edge(d, center, x) = d == 1 ? (center, x) : (x, center)

"""
    temporal_motif_counts(g::OrderedEdgeList, δ; strict = true)

Numbers of instances of all δ-temporal motifs with 3 edges on 2 or 3 nodes
(Paranjape, Benson and Leskovec, WSDM 2017), as a 6×6 matrix `M`: `M[i, j]` counts
the motif ``M_{i,j}`` of Figure 3 of the paper, in the same layout as the SNAP
implementation. If the first edge of the motif is ``a → b`` and ``c`` is the third
node, column `j` is given by the third edge and row `i` by the second edge:

| index | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| column `j` (3rd edge) | a → b | b → a | a → c | c → a | b → c | c → b |
| row `i` (2nd edge) | c → b | b → c | c → a | a → c | b → a | a → b |

so that `M[5:6, 1:2]` are the 2-node motifs, `M[1:2, 3:4]` and `M[3:4, 5:6]` the
triangles and the others the stars.

An instance is a sequence of three edges `e₁, e₂, e₃` of `g` with `t₃ - t₁ ≤ δ` that
maps onto the motif. With `strict = true` the time stamps must increase strictly
(the definition of the paper); with `strict = false` edges with equal time stamps are
ordered as in `edges(g)`. Only departure times are used and self loops are ignored.

Stars and 2-node motifs are counted in `O(m)` time with one sweep over the edges of
each node, and triangles in `O(m √τ)` time (`τ` triangles in the static graph) with
one sweep per pair of nodes, as in the paper; nodes and pairs are processed in
parallel.

# References

- A. Paranjape, A. R. Benson, and J. Leskovec. *Motifs in temporal networks.* WSDM, 2017. [DOI](https://doi.org/10.1145/3018661.3018731), [arXiv](https://arxiv.org/abs/1612.09259)
"""
function temporal_motif_counts(g::OrderedEdgeList, δ::Real; strict::Bool=true)
    δ >= 0 || throw(ArgumentError("δ must be non-negative"))
    n = num_nodes(g)
    es = edges(g)
    T = time_type(g)
    M = zeros(Int, 6, 6)
    # stars and 2-node motifs
    off, idx, nbr, dir = _incidence(g)
    stars = tmapreduce(.+, chunks(1:n; n=min(max(n, 1), 8 * Threads.nthreads())); init=zeros(Int, 32)) do us
        c = _StarCounter()
        loc = zeros(Int, n)
        ts = T[]
        for u in us
            _star_motifs!(c, ts, loc, off, idx, nbr, dir, es, u, δ, strict)
        end
        vcat(c.cnt_pre, c.cnt_post, c.cnt_mid, c.cnt_two)
    end
    for d1 in 1:2, d2 in 1:2, d3 in 1:2
        x = _i3(d1, d2, d3)
        M[_motif_cell(_motif_edge(d1, 1, 2), _motif_edge(d2, 1, 2), _motif_edge(d3, 1, 3))...] += stars[x]
        M[_motif_cell(_motif_edge(d1, 1, 3), _motif_edge(d2, 1, 2), _motif_edge(d3, 1, 2))...] += stars[8+x]
        M[_motif_cell(_motif_edge(d1, 1, 2), _motif_edge(d2, 1, 3), _motif_edge(d3, 1, 2))...] += stars[16+x]
        M[_motif_cell(_motif_edge(d1, 1, 2), _motif_edge(d2, 1, 2), _motif_edge(d3, 1, 2))...] += stars[24+x]
    end
    # triangles
    pairs, pair_off, edge_pair, adj_off, adj, adj_pair = _static_pairs(g)
    gpair, goff, _, qa, qb = _triangle_groups(n, pairs, pair_off, adj_off, adj, adj_pair)
    soff, sidx, scode = _triangle_streams(edge_pair, pair_off, gpair, goff, qa, qb)
    G = length(gpair)
    tris = tmapreduce(.+, chunks(1:G; n=min(max(G, 1), 8 * Threads.nthreads())); init=zeros(Int, 48)) do gs
        c = _TriangleCounter()
        ts = T[]
        for x in gs
            a, b = pairs[gpair[x]]
            _triangle_motifs!(c, ts, es, a, b, sidx, scode, soff[x]:soff[x+1]-1, goff[x+1] - goff[x], δ, strict)
        end
        vcat(c.cnt_pre, c.cnt_post, c.cnt_mid)
    end
    for d in 1:2, s in 1:2, d1 in 1:2, d2 in 1:2
        y = _i4(d, s, d1, d2)
        ab = d == 1 ? (1, 2) : (2, 1)
        x = _motif_edge(d1, s, 3)          # side s: node 1 (a) or 2 (b)
        z = _motif_edge(d2, 3 - s, 3)
        M[_motif_cell(x, z, ab)...] += tris[y]
        M[_motif_cell(ab, x, z)...] += tris[16+y]
        M[_motif_cell(x, ab, z)...] += tris[32+y]
    end
    return M
end

# ----------------------------------------------------------------------------------
# General motifs (Section 4.1 of the paper): enumerate the maps of the motif nodes to
# the graph nodes that realize the static graph H of the motif, and count with a
# sliding window the subsequences of the edges on the mapped pairs that spell the
# motif. counts[a, b] is the number of subsequences of the window matching edges a to
# b of the motif; when a block leaves the window, the subsequences starting in it are
# removed, when it enters, those ending in it are added.

mutable struct _SequenceCounter
    pattern::Vector{Int}    # static edge (label) of every motif edge
    counts::Matrix{Int}
    g::Vector{Int}          # edges of the current block per label
    total::Int
end

function _count_sequence!(c::_SequenceCounter, labels::Vector{Int}, ts::Vector, δ, strict::Bool)
    p = c.pattern
    l = length(p)
    C = c.counts
    fill!(C, 0)
    L = length(labels)
    start = 1
    i = 1
    @inbounds while i <= L
        j = _block_last(ts, i, L, strict)
        τ = ts[i]
        while ts[start] + δ < τ  # remove the earliest block
            s2 = _block_last(ts, start, L, strict)
            fill!(c.g, 0)
            for x in start:s2
                c.g[labels[x]] += 1
            end
            for a in l:-1:1
                ga = c.g[p[a]]
                ga == 0 && continue
                C[a, a] -= ga
                for b in a+1:l
                    (a == 1 && b == l) && continue
                    C[a, b] -= ga * C[a+1, b]
                end
            end
            start = s2 + 1
        end
        fill!(c.g, 0)
        for x in i:j
            c.g[labels[x]] += 1
        end
        for len in l-1:-1:1, a in 1:l-len  # subsequences a..b ending in the block, longest first
            b = a + len
            gb = c.g[p[b]]
            gb == 0 && continue
            if a == 1 && b == l
                c.total += gb * C[1, l-1]
            else
                C[a, b] += gb * C[a, b-1]
            end
        end
        if l == 1
            c.total += c.g[p[1]]
        else
            for b in 1:l
                C[b, b] += c.g[p[b]]
            end
        end
        i = j + 1
    end
    return c
end

"""
    temporal_motif_count(g::OrderedEdgeList, motif, δ; strict = true)

Number of instances of the δ-temporal motif `motif` (Paranjape, Benson and Leskovec,
2017), given as a vector of directed edges `(x, y)` on the nodes `1:k`, in their
temporal order: for example `[(1, 2), (2, 3), (3, 1)]` is a cyclic triangle and
`[(1, 2), (2, 1), (1, 2), (2, 1)]` a 4-edge exchange between two nodes. The static
graph of the motif must be connected. Instances and `strict` are as in
[`temporal_motif_counts`](@ref).

This is the general algorithm of the paper: it enumerates the embeddings of the
static graph of the motif in the static graph of `g` and counts the matching
subsequences of their edges with a sliding window, in `O(l² m')` time per embedding
(`l` motif edges, `m'` temporal edges on the embedding). The embeddings are
enumerated in parallel. It works for any motif but can be slow for motifs with
high-degree centers; use [`temporal_motif_counts`](@ref) for the motifs with 3 edges.

# References

- A. Paranjape, A. R. Benson, and J. Leskovec. *Motifs in temporal networks.* WSDM, 2017. [DOI](https://doi.org/10.1145/3018661.3018731), [arXiv](https://arxiv.org/abs/1612.09259)
"""
function temporal_motif_count(g::OrderedEdgeList, motif::AbstractVector{<:Tuple{Integer,Integer}}, δ::Real;
                              strict::Bool=true)
    δ >= 0 || throw(ArgumentError("δ must be non-negative"))
    l = length(motif)
    l >= 1 || throw(ArgumentError("the motif needs at least one edge"))
    mot = [(Int(x), Int(y)) for (x, y) in motif]
    k = maximum(max(x, y) for (x, y) in mot)
    all(x != y for (x, y) in mot) || throw(ArgumentError("motifs cannot have self loops"))
    Set(Iterators.flatten(mot)) == Set(1:k) || throw(ArgumentError("the motif nodes must be 1:$k"))
    H = unique(mot)
    pattern = [findfirst(==(e), H) for e in mot]
    # order of the motif nodes such that every node after the first is adjacent to an earlier one
    order = [1]
    while length(order) < k
        x = findfirst(y -> !(y in order) && any(e -> (e[1] == y && e[2] in order) || (e[2] == y && e[1] in order), H), 1:k)
        x === nothing && throw(ArgumentError("the static graph of the motif must be connected"))
        push!(order, x)
    end
    n = num_nodes(g)
    es = edges(g)
    T = time_type(g)
    # directed static graph: sorted out- and in-neighbors and the edges of every directed pair
    keys = [(Int(e.u), Int(e.v)) for e in es]
    perm = [i for i in sortperm(keys; alg=Base.Sort.DEFAULT_STABLE) if keys[i][1] != keys[i][2]]
    dpairs = Tuple{Int,Int}[]
    doff = Int[]
    for (r, i) in enumerate(perm)
        if r == 1 || keys[i] != keys[perm[r-1]]
            push!(dpairs, keys[i])
            push!(doff, r)
        end
    end
    push!(doff, length(perm) + 1)
    out = [Int[] for _ in 1:n]
    inn = [Int[] for _ in 1:n]
    aoff = zeros(Int, n + 1)  # the pairs (a, ⋅) are dpairs[aoff[a]+1:aoff[a+1]], sorted
    for (a, b) in dpairs
        push!(out[a], b)
        push!(inn[b], a)
        aoff[a+1] += 1
    end
    for a in 1:n
        aoff[a+1] += aoff[a]
    end
    foreach(sort!, inn)
    function pid(a, b)
        o = out[a]
        i = searchsortedfirst(o, b)
        return i <= length(o) && o[i] == b ? aoff[a] + i : 0
    end
    total = tmapreduce(+, chunks(1:n; n=min(max(n, 1), 8 * Threads.nthreads())); init=0) do roots
        c = _SequenceCounter(pattern, zeros(Int, l, l), zeros(Int, length(H)), 0)
        f = zeros(Int, k)
        used = falses(n)
        buf = Tuple{Int,Int}[]
        labels = Int[]
        ts = T[]
        function count_embedding()
            empty!(buf)
            for (h, (x, y)) in enumerate(H)
                q = pid(f[x], f[y])
                for r in doff[q]:doff[q+1]-1
                    push!(buf, (perm[r], h))
                end
            end
            sort!(buf)
            resize!(labels, length(buf))
            resize!(ts, length(buf))
            for (r, (i, h)) in enumerate(buf)
                labels[r] = h
                ts[r] = es[i].t
            end
            _count_sequence!(c, labels, ts, δ, strict)
        end
        function extend(level)
            if level > k
                count_embedding()
                return
            end
            x = order[level]
            # candidates: neighbors of the image of an earlier adjacent node
            e = H[findfirst(e -> (e[1] == x && e[2] in view(order, 1:level-1)) ||
                                 (e[2] == x && e[1] in view(order, 1:level-1)), H)]
            cands = e[1] == x ? inn[f[e[2]]] : out[f[e[1]]]
            for z in cands
                used[z] && continue
                f[x] = z
                ok = true
                for (a, b) in H
                    if (a == x || b == x) && (a == x || f[a] != 0) && (b == x || f[b] != 0) && pid(f[a], f[b]) == 0
                        ok = false
                        break
                    end
                end
                if ok
                    used[z] = true
                    extend(level + 1)
                    used[z] = false
                end
                f[x] = 0
            end
        end
        for r in roots
            f[order[1]] = r
            used[r] = true
            extend(2)
            used[r] = false
            f[order[1]] = 0
        end
        c.total
    end
    return total
end
