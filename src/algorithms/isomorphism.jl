# Isomorphisms of temporal graphs, following Heeg, Sauer, Mutzel and Scholtes, "Weisfeiler
# and Leman follow the arrow of time: expressive power of message passing in temporal
# event graphs" (arXiv:2505.24438).
#
# Three notions, from the weakest to the strongest:
# - time-aggregated isomorphism: isomorphism of the static graphs with the number of
#   temporal edges of every pair of nodes as edge label;
# - time-respecting path isomorphism (consistent event graph isomorphism): bijections of
#   the nodes and of the temporal edges that preserve all time-respecting paths, i.e.,
#   an isomorphism of the augmented event graphs (Theorem 1 of the paper); timestamps
#   can change as long as the time-respecting paths are the same;
# - time-concatenated isomorphism (timewise isomorphism of Wałęga and Rawson): the
#   static graphs with the sets of timestamps (relative to the first one) as labels.
# Isomorphisms are found exactly by individualization and refinement on the two graphs
# together (exponential in the worst case, fast in practice); the directed
# Weisfeiler–Leman refinement (D-WL) alone is the one-sided test of the paper.

# Directed graph with node and arc labels (integers), stored as arc lists.
struct _LabeledDigraph
    n::Int
    label::Vector{Int}
    out::Vector{Vector{Tuple{Int,Int}}}  # (neighbor, arc label)
    inn::Vector{Vector{Tuple{Int,Int}}}
end

function _labeled_digraph(n, label, arcs)
    out = [Tuple{Int,Int}[] for _ in 1:n]
    inn = [Tuple{Int,Int}[] for _ in 1:n]
    for (u, v, l) in arcs
        push!(out[u], (v, l))
        push!(inn[v], (u, l))
    end
    return _LabeledDigraph(n, label, out, inn)
end

"""
    augmented_event_graph(g::OrderedEdgeList, ti = time_interval(g); β = Inf)

The *augmented event graph* of `g` (Heeg, Sauer, Mutzel and Scholtes, 2025): the
temporal event graph (the [`edge_expansion`](@ref), with a node per temporal edge and
an arc `e → f` when `f` can follow `e` within the waiting time `β`), augmented with the
nodes of `g` and the arcs `u → e` and `e → v` for every temporal edge `e = (u, v, t, tt)`.
Its paths are the time-respecting paths of `g`, so two temporal graphs have the same
time-respecting paths up to renaming iff their augmented event graphs are isomorphic.

Returns `(graph, labels, edges)`: `graph` is a `Graphs.SimpleDiGraph` whose nodes
`1:n` are the nodes of `g` (label 0) and whose node `n + i` is the temporal edge
`edges(g)[edges[i]]` (label 1). The paper uses unit transition times and a maximum
time difference `δ` between consecutive edges; this is `tt = 1` and `β = δ - 1`.

# References

- F. Heeg, J. Sauer, P. Mutzel, and I. Scholtes. *Weisfeiler and Leman follow the arrow of time: expressive power of message passing in temporal event graphs.* 2025. [arXiv](https://arxiv.org/abs/2505.24438)
"""
function augmented_event_graph(g::OrderedEdgeList, ti=time_interval(g); β::Real=Inf)
    E = edge_expansion(g, ti; β=β)
    n, k = num_nodes(g), length(E.edges)
    h = Graphs.SimpleDiGraph(n + k)
    for e in Graphs.edges(E.graph)
        Graphs.add_edge!(h, n + Graphs.src(e), n + Graphs.dst(e))
    end
    es = edges(g)
    for (x, i) in enumerate(E.edges)
        Graphs.add_edge!(h, Int(es[i].u), n + x)
        Graphs.add_edge!(h, n + x, Int(es[i].v))
    end
    return (graph=h, labels=[fill(0, n); fill(1, k)], edges=E.edges)
end

function _augmented_digraph(g::OrderedEdgeList, ti, β)
    A = augmented_event_graph(g, ti; β=β)
    arcs = [(Graphs.src(e), Graphs.dst(e), 0) for e in Graphs.edges(A.graph)]
    return _labeled_digraph(Graphs.nv(A.graph), A.labels, arcs), A.edges
end

# Static graph of the pairs (u, v) with a label computed from their timestamps.
function _static_digraph(g::OrderedEdgeList{V,T}, ti, labelof) where {V,T}
    a, b = _interval(T, ti)
    times = Dict{Tuple{Int,Int},Vector{T}}()
    es = [e for e in edges(g) if a <= e.t && e.t + e.tt <= b]
    t0 = isempty(es) ? zero(T) : minimum(e.t for e in es)
    for e in es
        push!(get!(times, (Int(e.u), Int(e.v)), T[]), e.t - t0)
    end
    arcs = [(u, v, labelof(sort!(ts))) for ((u, v), ts) in times]
    return _labeled_digraph(num_nodes(g), zeros(Int, num_nodes(g)), arcs)
end

# 64-bit mixing (splitmix64).
@inline function _mix(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

# Directed Weisfeiler–Leman refinement of the disjoint union of graphs, with colors
# shared by all graphs: the new color of a node is determined by its color and the
# multisets of (color, arc label) of its in- and out-neighbors, combined with
# commutative sums of 64-bit hashes (no sorting; a collision, of probability about
# (n / 2^32)², could only make the refinement coarser). Returns the stable coloring
# (or the coloring after `rounds` rounds).
function _dwl_colors(gs::Vector{_LabeledDigraph}, init::Vector{Vector{Int}}; rounds=typemax(Int),
                     history::Union{Nothing,Vector}=nothing)
    colors = [copy(c) for c in init]
    ncol = length(unique(reduce(vcat, colors; init=Int[])))
    sig = [Vector{UInt64}(undef, g.n) for g in gs]
    r = 0
    @inbounds while r < rounds
        r += 1
        for (k, g) in enumerate(gs)
            c = colors[k]
            for v in 1:g.n
                so = UInt64(0)
                for (w, l) in g.out[v]
                    so += _mix(UInt64(c[w]) << 20 ⊻ UInt64(l))
                end
                si = UInt64(0)
                for (w, l) in g.inn[v]
                    si += _mix(UInt64(c[w]) << 20 ⊻ UInt64(l) ⊻ 0x5555555555555555)
                end
                sig[k][v] = _mix(_mix(UInt64(c[v])) ⊻ so) ⊻ _mix(si + 0x2545f4914f6cdd1d)
            end
        end
        # colors numbered by the sorted signatures: the same in every graph
        distinct = sort!(unique!(reduce(vcat, sig; init=UInt64[])))
        colors = [[searchsortedfirst(distinct, x) for x in sig[k]] for k in eachindex(gs)]
        history === nothing || push!(history, deepcopy(colors))
        length(distinct) == ncol && break
        ncol = length(distinct)
    end
    return colors
end

# Signature of node v of graph g under the coloring c (see _dwl_colors).
@inline function _signature(g::_LabeledDigraph, c::Vector{Int}, v::Int)
    so = UInt64(0)
    @inbounds for (w, l) in g.out[v]
        so += _mix(UInt64(c[w]) << 20 ⊻ UInt64(l))
    end
    si = UInt64(0)
    @inbounds for (w, l) in g.inn[v]
        si += _mix(UInt64(c[w]) << 20 ⊻ UInt64(l) ⊻ 0x5555555555555555)
    end
    return _mix(so) ⊻ _mix(si + 0x2545f4914f6cdd1d)
end

# Coloring of the two graphs during the search: colors, number of nodes of every
# color in both graphs (`size`) and in the first graph (`size1`), and a trail of the
# recolorings (graph, node, old color) to undo them when backtracking.
mutable struct _IsoState
    cs::NTuple{2,Vector{Int}}
    size::Vector{Int}
    size1::Vector{Int}
    trail::Vector{NTuple{3,Int}}
    unbalanced::Int   # number of colors x with 2 size1[x] != size[x]
end

@inline _isbalanced(st::_IsoState, x::Int) = 2 * st.size1[x] == st.size[x]

function _newcolor!(st::_IsoState)
    push!(st.size, 0)
    push!(st.size1, 0)
    return length(st.size)
end

# Give node w of graph k the color y, recording the change.
@inline function _recolor!(st::_IsoState, k::Int, w::Int, y::Int)
    x = st.cs[k][w]
    before = _isbalanced(st, x) + _isbalanced(st, y)
    push!(st.trail, (k, w, x))
    st.cs[k][w] = y
    st.size[x] -= 1
    st.size[y] += 1
    if k == 1
        st.size1[x] -= 1
        st.size1[y] += 1
    end
    st.unbalanced += before - (_isbalanced(st, x) + _isbalanced(st, y))
    return st
end

# Undo the recolorings after position `mark` of the trail and drop the colors after `ncolors`.
function _undo!(st::_IsoState, mark::Int, ncolors::Int)
    while length(st.trail) > mark
        k, w, x = pop!(st.trail)
        y = st.cs[k][w]
        before = _isbalanced(st, x) + _isbalanced(st, y)
        st.cs[k][w] = x
        st.size[y] -= 1
        st.size[x] += 1
        if k == 1
            st.size1[y] -= 1
            st.size1[x] += 1
        end
        st.unbalanced += before - (_isbalanced(st, x) + _isbalanced(st, y))
    end
    resize!(st.size, ncolors)
    resize!(st.size1, ncolors)
    return st
end

# Refine the coloring after the colors of the nodes in `changed` changed, until it is
# stable: only the neighbors of recolored nodes are examined, and a class is split by
# the signatures of its affected members (the others keep the old color; if the whole
# class is affected, its largest group does). Colors are assigned identically in both
# graphs, so the result is the same partition as a full refinement.
function _refine!(gs::NTuple{2,_LabeledDigraph}, st::_IsoState, changed::Vector{Tuple{Int,Int}})
    affected = Set{Tuple{Int,Int}}()
    groups = Dict{Int,Vector{Tuple{UInt64,Int,Int}}}()
    while !isempty(changed)
        empty!(affected)
        for (k, v) in changed, list in (gs[k].out[v], gs[k].inn[v]), (w, _) in list
            push!(affected, (k, w))
        end
        empty!(changed)
        empty!(groups)
        for (k, w) in affected
            push!(get!(groups, st.cs[k][w], Tuple{UInt64,Int,Int}[]), (_signature(gs[k], st.cs[k], w), k, w))
        end
        for x in sort!(collect(keys(groups)))
            members = sort!(groups[x])
            sigs = unique(first.(members))
            whole = length(members) == st.size[x]
            (whole && length(sigs) == 1) && continue  # the class does not split
            keep = whole ? argmax(sig -> count(m -> m[1] == sig, members), sigs) : nothing
            for sig in sigs
                sig == keep && continue
                y = _newcolor!(st)
                for (h, k, w) in members
                    h == sig || continue
                    _recolor!(st, k, w, y)
                    push!(changed, (k, w))
                end
            end
        end
    end
    return st
end

# An isomorphism from g1 to g2 (vector of images), or nothing: individualization and
# refinement (depth-first, with an explicit stack and a trail of recolorings), with an
# incremental refinement after every individualization.
function _find_isomorphism(g1::_LabeledDigraph, g2::_LabeledDigraph)
    g1.n == g2.n || return nothing
    n = g1.n
    sum(length, g1.out; init=0) == sum(length, g2.out; init=0) || return nothing
    arcs2 = Set((u, v, l) for u in 1:n for (v, l) in g2.out[u])
    function verify(π)
        for u in 1:n, (v, l) in g1.out[u]
            (π[u], π[v], l) in arcs2 || return false
        end
        return true
    end
    gs = (g1, g2)
    c1, c2 = _dwl_colors([g1, g2], [copy(g1.label), copy(g2.label)])
    ncol = maximum(vcat(c1, c2); init=0)
    size, size1 = zeros(Int, ncol), zeros(Int, ncol)
    for x in c1
        size[x] += 1
        size1[x] += 1
    end
    for x in c2
        size[x] += 1
    end
    st = _IsoState((c1, c2), size, size1, NTuple{3,Int}[], count(x -> 2 * size1[x] != size[x], 1:ncol))
    # frames: (node v of g1, candidates in g2, next candidate, trail mark, number of colors)
    stack = Tuple{Int,Vector{Int},Int,Int,Int}[]
    descend = true
    while true
        if descend && st.unbalanced == 0
            best = 0  # the smallest class with more than one node per graph
            for x in eachindex(st.size)
                st.size[x] > 2 && (best == 0 || st.size[x] < st.size[best]) && (best = x)
            end
            if best == 0
                pos = zeros(Int, length(st.size))
                for (v, x) in enumerate(c2)
                    pos[x] = v
                end
                π = [pos[x] for x in c1]
                verify(π) && return π
            else
                push!(stack, (findfirst(==(best), c1), findall(==(best), c2), 1, length(st.trail), length(st.size)))
            end
        end
        # try the next candidate of the deepest frame
        isempty(stack) && return nothing
        v, ws, i, mark, nc = stack[end]
        _undo!(st, mark, nc)
        if i > length(ws)
            pop!(stack)
            descend = false
            continue
        end
        stack[end] = (v, ws, i + 1, mark, nc)
        y = _newcolor!(st)
        _recolor!(st, 1, v, y)
        _recolor!(st, 2, ws[i], y)
        _refine!(gs, st, [(1, v), (2, ws[i])])
        descend = true
    end
end

const _ISO_KINDS = (:event, :aggregated, :concatenated)

function _iso_graphs(g1::OrderedEdgeList, g2::OrderedEdgeList, kind::Symbol, β, ti1, ti2)
    kind in _ISO_KINDS || throw(ArgumentError("kind must be :event, :aggregated or :concatenated"))
    if kind == :event
        h1, e1 = _augmented_digraph(g1, ti1, β)
        h2, e2 = _augmented_digraph(g2, ti2, β)
        return h1, h2, e1, e2
    end
    ids = Dict{Any,Int}()  # labels shared by both graphs
    labelof = kind == :aggregated ? (ts -> length(ts)) : (ts -> get!(ids, ts, length(ids) + 1))
    return _static_digraph(g1, ti1, labelof), _static_digraph(g2, ti2, labelof), nothing, nothing
end

"""
    temporal_isomorphism(g1::OrderedEdgeList, g2::OrderedEdgeList; kind = :event, β = Inf,
                         ti1 = time_interval(g1), ti2 = time_interval(g2))

An isomorphism between the temporal graphs `g1` and `g2` of the given `kind`
(Heeg, Sauer, Mutzel and Scholtes, 2025), or `nothing` if there is none:

- `:event` — *time-respecting path isomorphism* (equivalently, consistent event graph
  isomorphism): bijections of the nodes and of the temporal edges that map every
  temporal edge `(u, v, t)` to an edge between the images of `u` and `v` and preserve
  all time-respecting paths with waiting time at most `β`. Timestamps may change as
  long as the time-respecting paths do not. Computed as an isomorphism of the
  [`augmented_event_graph`](@ref)s (Theorem 1 of the paper). Returns
  `(nodes, edges)`: node `v` of `g1` maps to `nodes[v]`, and edge `i` of `edges(g1)`
  to edge `edges[i]` of `edges(g2)` (0 for edges outside the time interval).
- `:aggregated` — *time-aggregated isomorphism*: the static graphs with the number of
  temporal edges of every pair of nodes are isomorphic. Returns `(nodes,)`.
- `:concatenated` — *time-concatenated isomorphism* (timewise isomorphism of Wałęga and
  Rawson): the static graphs with the sets of timestamps relative to the first one are
  isomorphic. Returns `(nodes,)`.

Time-concatenated isomorphic graphs are event isomorphic, and event isomorphic graphs
are time-aggregated isomorphic (Theorem 2). Graph isomorphism has no known polynomial
algorithm; this exact test refines the colors of both graphs together with the
directed Weisfeiler–Leman algorithm and individualizes nodes when the refinement
stops, which is fast unless the graphs are very symmetric.

# References

- F. Heeg, J. Sauer, P. Mutzel, and I. Scholtes. *Weisfeiler and Leman follow the arrow of time: expressive power of message passing in temporal event graphs.* 2025. [arXiv](https://arxiv.org/abs/2505.24438)
- P. A. Wałęga and M. Rawson. *Expressive power of temporal message passing.* AAAI, 2025. [DOI](https://doi.org/10.1609/aaai.v39i20.35396)
"""
function temporal_isomorphism(g1::OrderedEdgeList, g2::OrderedEdgeList; kind::Symbol=:event, β::Real=Inf,
                              ti1=time_interval(g1), ti2=time_interval(g2))
    h1, h2, e1, e2 = _iso_graphs(g1, g2, kind, β, ti1, ti2)
    π = _find_isomorphism(h1, h2)
    π === nothing && return nothing
    n1, n2 = num_nodes(g1), num_nodes(g2)
    kind == :event || return (nodes=π,)
    emap = zeros(Int, num_edges(g1))
    for (x, i) in enumerate(e1)
        emap[i] = e2[π[n1+x]-n2]
    end
    return (nodes=π[1:n1], edges=emap)
end

"""
    is_temporally_isomorphic(g1::OrderedEdgeList, g2::OrderedEdgeList; kind = :event, β = Inf)

Whether [`temporal_isomorphism`](@ref) finds an isomorphism of the given `kind`.
"""
is_temporally_isomorphic(g1::OrderedEdgeList, g2::OrderedEdgeList; kwargs...) =
    temporal_isomorphism(g1, g2; kwargs...) !== nothing

"""
    temporal_wl_equivalent(g1::OrderedEdgeList, g2::OrderedEdgeList; β = Inf, kind = :event)

The temporal Weisfeiler–Leman test of Heeg, Sauer, Mutzel and Scholtes (2025): run the
directed Weisfeiler–Leman color refinement (D-WL) on the [`augmented_event_graph`](@ref)s
of both graphs until the colors are stable and compare the color histograms. It is a
one-sided test: `false` proves that the graphs are not time-respecting path
isomorphic, `true` means that D-WL (and message passing neural networks on the
augmented event graph) cannot distinguish them. With `kind = :aggregated` or
`:concatenated` the refinement runs on the corresponding labeled static graphs.
Runs in `O(n + a)` time per round (plus sorting the colors) for `a` arcs; neighbor
multisets are compared through 64-bit hashes, so a collision (of negligible
probability) could only make the test weaker.

# References

- F. Heeg, J. Sauer, P. Mutzel, and I. Scholtes. *Weisfeiler and Leman follow the arrow of time: expressive power of message passing in temporal event graphs.* 2025. [arXiv](https://arxiv.org/abs/2505.24438)
- E. Rossi, B. Charpentier, F. Di Giovanni, F. Frasca, S. Günnemann, and M. M. Bronstein. *Edge directionality improves learning on heterophilic graphs.* Learning on Graphs Conference (LoG), 2023. [PMLR](https://proceedings.mlr.press/v231/rossi24a.html)
"""
function temporal_wl_equivalent(g1::OrderedEdgeList, g2::OrderedEdgeList; β::Real=Inf, kind::Symbol=:event,
                                ti1=time_interval(g1), ti2=time_interval(g2))
    h1, h2, _, _ = _iso_graphs(g1, g2, kind, β, ti1, ti2)
    h1.n == h2.n || return false
    c1, c2 = _dwl_colors([h1, h2], [h1.label, h2.label])
    return sort(c1) == sort(c2)
end

"""
    temporal_wl_kernel(graphs; iterations = 3, β = Inf, augmented = true)

The Weisfeiler–Leman subtree kernel of a collection of temporal graphs on their
[`augmented_event_graph`](@ref)s (Heeg, Sauer, Mutzel and Scholtes, 2025): every graph
is represented by the numbers of nodes of every color after `0, 1, …, iterations`
rounds of directed color refinement, with colors shared by all graphs, and `K[i, j]`
is the dot product of these vectors. With `augmented = false` the refinement runs on
the temporal event graphs alone, as in the temporal graph kernel of Oettershagen,
Kriege, Morris and Mutzel (2020), which ignores which nodes the temporal edges join.
Returns `(kernel, features)` where `features[i]` maps colors to counts.

# References

- F. Heeg, J. Sauer, P. Mutzel, and I. Scholtes. *Weisfeiler and Leman follow the arrow of time: expressive power of message passing in temporal event graphs.* 2025. [arXiv](https://arxiv.org/abs/2505.24438)
- L. Oettershagen, N. M. Kriege, C. Morris, and P. Mutzel. *Temporal graph kernels for classifying dissemination processes.* SDM, 2020. [DOI](https://doi.org/10.1137/1.9781611976236.56), [arXiv](https://arxiv.org/abs/1911.05496)
"""
function temporal_wl_kernel(graphs::AbstractVector{<:OrderedEdgeList}; iterations::Integer=3, β::Real=Inf,
                            augmented::Bool=true)
    hs = _LabeledDigraph[]
    for g in graphs
        if augmented
            push!(hs, first(_augmented_digraph(g, time_interval(g), β)))
        else
            E = edge_expansion(g; β=β)
            arcs = [(Graphs.src(e), Graphs.dst(e), 0) for e in Graphs.edges(E.graph)]
            push!(hs, _labeled_digraph(Graphs.nv(E.graph), zeros(Int, Graphs.nv(E.graph)), arcs))
        end
    end
    history = Vector{Vector{Vector{Int}}}()
    init = [h.label for h in hs]
    push!(history, init)
    _dwl_colors(hs, init; rounds=iterations, history=history)
    while length(history) < iterations + 1  # the refinement stabilized early
        push!(history, history[end])
    end
    features = [Dict{Tuple{Int,Int},Int}() for _ in hs]
    for (r, colors) in enumerate(history), (k, c) in enumerate(colors), x in c
        features[k][(r, x)] = get(features[k], (r, x), 0) + 1
    end
    N = length(hs)
    K = zeros(Int, N, N)
    for i in 1:N, j in i:N
        s = 0
        for (key, x) in features[i]
            s += x * get(features[j], key, 0)
        end
        K[i, j] = K[j, i] = s
    end
    return (kernel=K, features=features)
end
