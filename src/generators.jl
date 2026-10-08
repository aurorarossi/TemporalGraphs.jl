# Generators of random and structured temporal graphs: the classes and constructions
# listed in the Temporal Graph Wiki (temporalgraph.notion.site). Undirected graphs are
# stored with both directions of every edge, as by `make_undirected`; all times are
# integers.

const _GenEdge = TemporalEdge{Int32,Int}

# Positions in 1:N chosen independently with probability p, in increasing order, by
# geometric jumps: O(1 + pN) expected time.
function _bernoulli_positions!(out::Vector{Int}, rng::AbstractRNG, N::Integer, p::Real)
    empty!(out)
    (p == 0 || N <= 0) && return out
    p == 1 && return append!(out, 1:N)
    c = log1p(-p)
    i = 0
    while true
        j = log(1 - rand(rng)) / c  # failures before the next success
        j >= N - i && return out
        i += 1 + floor(Int, j)
        push!(out, i)
    end
end

function _check_probability(p)
    0 <= p <= 1 || throw(ArgumentError("p must be in [0, 1]"))
    return nothing
end

function _check_transition_time(tt)
    tt >= 0 || throw(ArgumentError("the transition time must be non-negative"))
    return nothing
end

# The temporal edge (u, v, t, tt), and (v, u, t, tt) if undirected.
function _push_edge!(es::Vector{_GenEdge}, u, v, t, tt, directed::Bool)
    push!(es, _GenEdge(u, v, t, tt))
    directed || u == v || push!(es, _GenEdge(v, u, t, tt))
    return es
end

"""
    random_temporal_graph(n, p, T; directed = false, transition_time = 1,
                          rng = Random.default_rng())

A random temporal graph with `n` nodes whose snapshots at the times `1, ..., T` are
independent Erdős–Rényi graphs: every pair of distinct nodes (ordered pair if
`directed`) has an edge at every time with probability `p`, independently. Undirected
edges are stored in both directions. The default transition time 1 gives strict paths
(at most one edge per snapshot); use `transition_time = 0` for non-strict paths.

The pairs are drawn by geometric jumps, in `O(nT + m)` expected time for `m` edges.
The time interval is `(1, T + transition_time)`.
"""
function random_temporal_graph(n::Integer, p::Real, T::Integer; directed::Bool=false, transition_time::Integer=1,
                               rng::AbstractRNG=Random.default_rng())
    n >= 0 || throw(ArgumentError("number of nodes must be non-negative"))
    T >= 0 || throw(ArgumentError("T must be non-negative"))
    _check_probability(p)
    _check_transition_time(transition_time)
    es = _GenEdge[]
    pos = Int[]
    for t in 1:T, u in 1:n
        # row u: the nodes v ≠ u (directed) or v > u (undirected)
        _bernoulli_positions!(pos, rng, directed ? n - 1 : n - u, p)
        for k in pos
            v = directed ? (k < u ? k : k + 1) : u + k
            _push_edge!(es, u, v, t, transition_time, directed)
        end
    end
    return OrderedEdgeList(n, es, (1, T + transition_time))
end

"""
    random_simple_temporal_graph(n, p; directed = false, rng = Random.default_rng())

A *random simple temporal graph* (Casteigts, Raskin, Renken and Zamaraev, 2024): an
Erdős–Rényi graph `G(n, p)` (directed if `directed`) whose `m` edges get the times
`1, ..., m` in uniformly random order, one per edge. The labeling is *simple* (one
time per edge) and *proper* (adjacent edges have different times), so the temporal
paths are the paths with increasing times; transition times are 1. With `p = 1` this
is a random temporal clique.

Casteigts et al. show sharp thresholds for temporal connectivity at `p = log n / n`
(some node reaches all others) and `p = 3 log n / n` (all pairs reach each other).

# References

- A. Casteigts, M. Raskin, M. Renken, and V. Zamaraev. *Sharp thresholds in random simple temporal graphs.* SIAM Journal on Computing, 2024. [DOI](https://doi.org/10.1137/22M1511916), [arXiv](https://arxiv.org/abs/2011.03738)
"""
function random_simple_temporal_graph(n::Integer, p::Real; directed::Bool=false, rng::AbstractRNG=Random.default_rng())
    n >= 0 || throw(ArgumentError("number of nodes must be non-negative"))
    _check_probability(p)
    pairs = Tuple{Int,Int}[]
    pos = Int[]
    for u in 1:n
        _bernoulli_positions!(pos, rng, directed ? n - 1 : n - u, p)
        for k in pos
            push!(pairs, (u, directed ? (k < u ? k : k + 1) : u + k))
        end
    end
    labels = Random.randperm(rng, length(pairs))
    es = _GenEdge[]
    for ((u, v), t) in zip(pairs, labels)
        _push_edge!(es, u, v, t, 1, directed)
    end
    return OrderedEdgeList(n, es)
end

# k distinct integers of 1:T, uniformly at random (Floyd's algorithm).
function _distinct_times(rng::AbstractRNG, T::Int, k::Int)
    S = Set{Int}()
    for j in T-k+1:T
        t = rand(rng, 1:j)
        push!(S, t in S ? j : t)
    end
    return S
end

"""
    random_temporal_labeling(G::Graphs.AbstractGraph, T; labels = 1, transition_time = 1,
                             rng = Random.default_rng())

The temporal graph of the static graph `G` in which every edge gets `labels` distinct
times drawn uniformly at random from `1, ..., T`, in both directions if `G` is
undirected. Together with the generators of Graphs.jl it gives random temporal paths,
stars, trees, grids and cliques, e.g. `random_temporal_labeling(Graphs.star_graph(10), 5)`
or, for a random spanning tree with random times,
`random_temporal_labeling(Graphs.uniform_tree(n), T)`. With `labels = 1` the result is
simple (one time per edge). The time interval is `(1, T + transition_time)`.
"""
function random_temporal_labeling(G::Graphs.AbstractGraph, T::Integer; labels::Integer=1, transition_time::Integer=1,
                                  rng::AbstractRNG=Random.default_rng())
    1 <= labels <= T || throw(ArgumentError("labels must be between 1 and T"))
    _check_transition_time(transition_time)
    directed = Graphs.is_directed(G)
    es = _GenEdge[]
    for e in Graphs.edges(G), t in sort!(collect(_distinct_times(rng, Int(T), Int(labels))))
        _push_edge!(es, Graphs.src(e), Graphs.dst(e), t, transition_time, directed)
    end
    return OrderedEdgeList(Graphs.nv(G), es, (1, T + transition_time))
end

"""
    round_robin_temporal_clique(n)

The *round-robin temporal clique* on `n` nodes: the undirected complete graph whose
edges are scheduled by the circle method of round-robin tournaments. Round `r` (time
`r`) is a matching, and every pair of nodes meets exactly once, in `n - 1` rounds for
even `n` and `n` rounds for odd `n` (one node rests in each round). The labeling is
simple and proper; transition times are 1.
"""
function round_robin_temporal_clique(n::Integer)
    n >= 0 || throw(ArgumentError("number of nodes must be non-negative"))
    N = isodd(n) ? n + 1 : n  # node N is fixed, the others rotate; with odd n it is a rest
    es = _GenEdge[]
    for r in 0:N-2
        pairs = [(N, r + 1); [(mod(r + i, N - 1) + 1, mod(r - i, N - 1) + 1) for i in 1:N÷2-1]]
        for (u, v) in pairs
            (u <= n && v <= n) && _push_edge!(es, u, v, r + 1, 1, false)
        end
    end
    return OrderedEdgeList(n, es)
end

"""
    temporal_hypercube(d)

The `d`-dimensional hypercube on `2^d` nodes whose edges along dimension `k` have
time `k` (`k = 1, ..., d`), undirected, with transition times 1. Every pair of nodes is
joined by exactly one temporal path (fix the differing bits in increasing order), so
the graph is temporally connected and every one of its `d 2^(d-1)` edges is needed: no
temporal spanner is smaller (Kempe, Kleinberg and Kumar, 2002). Node `x + 1` is the
bit string `x`.

# References

- D. Kempe, J. Kleinberg, and A. Kumar. *Connectivity and inference problems for temporal networks.* Journal of Computer and System Sciences 64(4), 2002. [DOI](https://doi.org/10.1006/jcss.2002.1829)
"""
function temporal_hypercube(d::Integer)
    0 <= d <= 30 || throw(ArgumentError("d must be between 0 and 30"))
    es = _GenEdge[]
    for k in 0:d-1, x in 0:2^d-1
        x & (1 << k) == 0 && _push_edge!(es, x + 1, (x | (1 << k)) + 1, k + 1, 1, false)
    end
    return OrderedEdgeList(2^d, es)
end

"""
    temporal_knodel_graph(d)

The Knödel graph `W(d, 2^d)` on `2^d` nodes with the edges of dimension `k` at time
`k` (`k = 1, ..., d`), undirected, with transition times 1. Nodes `j + 1` and
`2^(d-1) + j + 1` form the two sides `j = 0, ..., 2^(d-1) - 1`; dimension `k` joins
`j` on the first side to `j + 2^(k-1) - 1 (mod 2^(d-1))` on the second. If in round
`k` every node exchanges what it knows along dimension `k`, all nodes know everything
after `d` rounds, the minimum possible for `2^d` nodes (Knödel, 1975). So the graph is
temporally connected; its labeling is proper, and since the knowledge of every node
doubles in every round, every pair of nodes is joined by exactly one temporal path and
every edge is needed.

# References

- W. Knödel. *New gossips and telephones.* Discrete Mathematics 13(1), 1975. [DOI](https://doi.org/10.1016/0012-365X(75)90090-4)
- G. Fertin and A. Raspaud. *A survey on Knödel graphs.* Discrete Applied Mathematics 137(2), 2004. [DOI](https://doi.org/10.1016/S0166-218X(03)00260-9)
"""
function temporal_knodel_graph(d::Integer)
    0 <= d <= 30 || throw(ArgumentError("d must be between 0 and 30"))
    d == 0 && return OrderedEdgeList(1, _GenEdge[])
    h = 2^(d - 1)
    es = _GenEdge[]
    for k in 0:d-1, j in 0:h-1
        _push_edge!(es, j + 1, h + mod(j + 2^k - 1, h) + 1, k + 1, 1, false)
    end
    return OrderedEdgeList(2^d, es)
end
