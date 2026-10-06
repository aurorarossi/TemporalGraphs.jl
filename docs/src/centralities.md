# Centralities and statistics

In the formulas below ``d(u, v)`` is a temporal distance of the chosen
[`DistanceType`](@ref) and ``n`` the number of nodes.

## Closeness

The temporal closeness of ``u`` is

```math
C(u) = \sum_{v :\, 0 < d(u, v) < ∞} \frac{1}{d(u, v)}.
```

For [`EarliestArrival`](@ref) ``d(u, v)`` is the arrival time at ``v``; for
[`LatestDeparture`](@ref) it is the end of the time window minus the latest time one
can leave ``u`` towards ``v``.

- [`temporal_closeness`](@ref)`(g, s, dt)` computes it for one node and
  `temporal_closeness(g, dt)` for all nodes, in parallel.
- [`compute_topk_closeness`](@ref) returns the ``k`` nodes with the highest
  closeness. For fastest and shortest paths it uses the pruning algorithm of
  Oettershagen and Mutzel (2020), which stops the computation for a node as soon as an
  upper bound of its closeness drops below the ``k``-th best value found so far.
- [`temporal_closeness_approximation`](@ref) estimates the fastest path closeness of
  all nodes from ``h`` random samples on the reversed graph; the estimate is unbiased.

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(4, [(1, 2, 1, 1), (2, 3, 2, 1), (3, 4, 3, 1), (1, 4, 9, 1)])) # hide
```

```jldoctest
julia> using TemporalGraphs

julia> g = OrderedEdgeList(4, [(1, 2, 1, 1), (2, 3, 2, 1), (3, 4, 3, 1), (1, 4, 9, 1)]);

julia> temporal_closeness(g, MinimumHops())
4-element Vector{Float64}:
 2.5
 1.5
 1.0
 0.0

julia> compute_topk_closeness(g, 1, Fastest())
1-element Vector{Tuple{Int64, Float64}}:
 (1, 2.5)
```

### Closeness integrated over time

[`temporal_harmonic_closeness`](@ref) follows Crescenzi, Magnien and Marino (2020):
instead of fixing a departure time, it averages over all departure times ``τ`` in
the time window ``[a, b]``:

```math
C(s) = \frac{1}{(n-1)(b-a)} \sum_{v \neq s} \int_a^b \frac{dτ}{ℓ_τ(s, v)},
```

where the *latency* ``ℓ_τ(s, v)`` is the earliest arrival time at ``v`` when leaving
``s`` at time ``τ`` or later, minus ``τ`` (unreachable nodes contribute 0).
Transition times must be positive. The latency is a step function of ``τ`` whose
steps are the Pareto optimal (earliest arrival, latest departure) pairs of the paths
from ``s`` to ``v``; these are computed for all ``v`` in one scan of the edges, so
the integral is exact and the closeness of all nodes takes ``O(nM)`` time (in
parallel). Real-valued times are supported.

[`temporal_harmonic_closeness_approximation`](@ref) estimates the closeness of all
nodes from random target nodes: the contributions of a target to all sources come
from one scan of the reversed temporal graph. The estimates are unbiased and most
accurate for the nodes with high closeness, which makes them suited to find the
top-k nodes.

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1), (1, 3, 5, 1)], (0, 6))) # hide
```

```jldoctest
julia> using TemporalGraphs

julia> g = OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1), (1, 3, 5, 1)], (0, 6));

julia> temporal_harmonic_closeness(g)
3-element Vector{Float64}:
 0.22567085009185084
 0.09155102405567582
 0.0
```

For node 1 the latency to node 2 is ``2 - τ`` for ``τ ≤ 1``, and the latency to
node 3 is ``3 - τ`` for ``τ ≤ 1`` (through node 2) and ``6 - τ`` for ``1 < τ ≤ 5``
(direct edge), so ``C(1) = (\log 2 + \log \frac{3}{2} + \log 5) / 12``.

The program of the paper uses unit transition times and normalizes by the span of
the time stamps `t_max - t_min`; with the default window `(t_min, t_max + 1)` its
values are those of [`temporal_harmonic_closeness`](@ref) multiplied by
`(t_max + 1 - t_min) / (t_max - t_min)`.

## Eccentricity, diameter and efficiency

- [`temporal_eccentricity`](@ref): the largest finite distance from a node.
- [`temporal_diameter`](@ref): the largest eccentricity.
- [`temporal_efficiency`](@ref): ``\frac{1}{n(n-1)} \sum_u C(u)``.

## Temporal betweenness

[`temporal_betweenness`](@ref) computes the betweenness of all nodes for optimal
temporal *walks* (Brunelli, Crescenzi and Viennot, 2024):

```math
B(u) = \sum_{s \neq u,\; t \neq u,\; s \neq t} \frac{σ_{s,t}(u)}{σ_{s,t}},
```

where ``σ_{s,t}`` is the number of optimal walks from ``s`` to ``t`` and
``σ_{s,t}(u)`` the number of those passing through ``u``, counted with multiplicity
(a walk through ``u`` twice counts twice). Walks can be *restless*: with a waiting
constraint `β`, the edge after an edge ``e`` must depart within `β` time units of the
arrival of ``e``.

| Criterion | Optimal walks |
|---|---|
| [`MinimumHops`](@ref) | fewest edges (shortest) |
| [`EarliestArrival`](@ref) | earliest arrival (foremost) |
| [`LatestDeparture`](@ref) | latest departure (latest) |
| [`Fastest`](@ref) | minimum duration |
| [`ShortestForemost`](@ref) | foremost, then fewest edges |
| [`ShortestLatest`](@ref) | latest, then fewest edges |
| [`ShortestFastest`](@ref) | fastest, then fewest edges |
| [`PrefixForemost`](@ref) | prefix foremost paths (no waiting constraint) |

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(4, [(1, 2, 1, 1), (2, 3, 3, 1), (3, 4, 4, 1), (1, 4, 9, 1), (2, 4, 7, 1)])) # hide
```

```jldoctest
julia> using TemporalGraphs

julia> g = OrderedEdgeList(4, [(1, 2, 1, 1), (2, 3, 3, 1), (3, 4, 4, 1), (1, 4, 9, 1), (2, 4, 7, 1)]);

julia> temporal_betweenness(g)                       # shortest walks
4-element Vector{Float64}:
 0.0
 1.0
 0.0
 0.0

julia> temporal_betweenness(g, EarliestArrival())    # foremost walks
4-element Vector{Float64}:
 0.0
 2.0
 2.0
 0.0

julia> temporal_betweenness(g, EarliestArrival(); β=0)   # no waiting allowed
4-element Vector{Float64}:
 0.0
 0.0
 1.0
 0.0
```

With [`MinimumHops`](@ref) and no waiting constraint this is the *shortest temporal
betweenness* of Buß et al. (2020) used in TSBProxy, since shortest strict walks are
paths.

**Algorithm.** For every source the algorithm makes one forward pass over the edges
by departure time, counting the optimal walks that end with each edge, and one
backward pass accumulating the dependencies of the edges. Without waiting constraint
the best predecessors of the outgoing edges of a node improve monotonically, so the
successors of an edge form a contiguous block and each edge is processed in constant
time: ``O(nM)`` in total. With a waiting constraint the predecessors form a sliding
window of arrival times (a monotone deque), and ``O(nM \log M)`` time is needed.
Sources are processed in parallel.

**Numbers of walks.** The number of optimal walks can be astronomically large, in
particular for foremost, latest and fastest walks, which are rarely unique. Counts
are kept in `Float64`; if they overflow, the computation is repeated with `BigFloat`
and a warning is issued. Pass `count_type = BigFloat` to skip the first attempt, or
`count_type = Rational{BigInt}` to count exactly; the betweenness values are returned
as `Float64` in every case. The accumulations only add
non-negative numbers (no differences of prefix sums), so the results are accurate
even when the counts span hundreds of orders of magnitude.

**Strict and non-strict walks.** With positive transition times walks are *strict*:
on integer times, transition time 1 gives the strict walks of Buß et al. and of
Brunelli et al., where consecutive edges have increasing time stamps. Transition time
0 gives the *non-strict* walks of Zhang et al. (2024), where consecutive edges can
have the same time stamp; both kinds of edges can be mixed. The edges with
transition time 0 and the same time stamp are processed as a block: by Dijkstra's
algorithm on the numbers of edges for the shortest criteria, in topological order for
the others. If these edges form a cycle that some walk reaches, a walk can go around
the cycle any number of times without changing its arrival or departure time, so
there are infinitely many foremost, latest or fastest walks and an `ArgumentError` is
thrown; the shortest criteria are always defined.

The non-strict graph, with transition times 0, and the strict one:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 5, 0), (2, 3, 5, 0)])) # hide
```

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 5, 1), (2, 3, 5, 1)])) # hide
```

```jldoctest
julia> using TemporalGraphs

julia> g = OrderedEdgeList(3, [(1, 2, 5, 0), (2, 3, 5, 0)]);   # non-strict: 1 → 2 → 3 at time 5

julia> temporal_betweenness(g)
3-element Vector{Float64}:
 0.0
 1.0
 0.0

julia> temporal_betweenness(OrderedEdgeList(3, [(1, 2, 5, 1), (2, 3, 5, 1)]))   # strict
3-element Vector{Float64}:
 0.0
 0.0
 0.0
```

**Latest walks.** Latest walks are the only optimal walks that can visit their target
more than once. As in Fact 3 of Brunelli et al., these earlier visits count as
passages through the target. Note that on the example of Figure 1 of the paper the
program TWBC and Table 2 give values for the latest and shortest latest criteria
that differ from the definition (and from a brute-force enumeration of the walks);
all other criteria agree exactly with TWBC.

### Proxies and approximation

Exact betweenness needs one computation per node. Becker, Crescenzi, Cruciani and
Kodric (2023) compare cheaper proxies for the shortest temporal betweenness ranking:

- [`temporal_betweenness`](@ref)`(g, PrefixForemost())`: the prefix foremost
  betweenness of Buß et al., for paths whose every prefix is foremost. The prefix
  foremost paths from a source form a DAG ordered by earliest arrival times, so this
  takes ``O(M)`` time per source.
- [`temporal_ego_betweenness`](@ref): the betweenness of each node in its *ego
  network*, the subgraph induced by the node and its neighbors (any criterion).
- [`temporal_pass_through_degree`](@ref): the square root of the number of pairs of
  neighbors that are temporally connected through the node, in ``O(M \log M)`` time.
- [`temporal_betweenness_approximation`](@ref): the ONBRA estimator of Santoro and
  Sarpe (2022), which samples pairs of nodes and returns unbiased estimates of the
  normalized betweenness ``B(u) / (n(n-1))`` together with an error bound that holds
  with probability ``1 - η`` (empirical Bernstein bound).

- [`temporal_betweenness_mantra`](@ref): the MANTRA progressive sampling algorithm of
  Cruciani (2024), for shortest, shortest foremost and prefix foremost paths. It
  returns estimates of the normalized betweenness that are within a target error
  `ε` of the exact values with probability ``1 - δ``. The number of samples is not
  fixed in advance: after a bootstrap phase that bounds the variance of the
  estimators and the average path length, the sample grows until a bound on the
  maximum error computed from Monte Carlo empirical Rademacher averages drops below
  `ε`. Pairs of nodes (as in ONBRA) or single sources can be sampled.

- [`temporal_betweenness_atbc`](@ref): the ATBC progressive sampling algorithm of
  Zhang et al. (2024), for shortest, foremost, fastest and shortest foremost walks,
  strict or non-strict. It samples pairs of nodes until the bound of Riondato and
  Upfal (ABRA) on the maximum error, computed from the empirical Rademacher average
  of the sample, is at most `ε`, and returns the same named tuple as MANTRA.

```julia
r = temporal_betweenness_mantra(g, MinimumHops(); ε = 0.005, δ = 0.1)
r.estimates      # normalized betweenness of all nodes
r.error_bound    # every estimate is within this of the exact value w.p. 1 - δ
r.samples        # number of samples used
```

The static betweenness of the underlying graph, another proxy, is available through
Graphs.jl: `Graphs.betweenness_centrality(static_graph(g))`.

### Distance statistics

[`temporal_distance_statistics`](@ref) computes the diameter, the effective diameter
(the `τ` quantile of the lengths), the temporal connectivity rate (the fraction of
connected ordered pairs) and the average length of the shortest, shortest foremost
or prefix foremost paths, either exactly or, as in Algorithm 2 of Cruciani (2024),
from a sample of source nodes.

```jldoctest
julia> using TemporalGraphs

julia> g = OrderedEdgeList(4, [(1, 2, 1, 1), (2, 3, 2, 1), (3, 4, 3, 1), (1, 4, 9, 1)]);

julia> temporal_distance_statistics(g)
(diameter = 2, effective_diameter = 2, connectivity_rate = 0.5, average_distance = 1.3333333333333333)
```

## Edge betweenness

[`temporal_edge_betweenness`](@ref) computes, for every temporal edge, the sum over
all pairs of edges ``(e, f)`` of the fraction of minimum hop temporal paths from
``e`` to ``f`` that pass through it. It runs Brandes' algorithm on the
[`DirectedLineGraph`](@ref), whose nodes are the temporal edges.

## Walk based centralities

A *temporal walk* is a sequence of edges ``e_1, \dots, e_k`` with ``e_{i+1}``
leaving the head of ``e_i`` and strictly increasing time stamps.

- [`temporal_katz_centrality`](@ref)`(g, β)`: ``\sum_{W} β^{|W|}`` over all walks
  ending at the node (Béres et al., 2018).
- [`temporal_walk_centrality`](@ref)`(g, α, β)`: the weighted number of walks that
  pass through a node, with incoming walks weighted by ``α`` and outgoing walks by
  ``β`` per edge (Oettershagen, Mutzel and Kriege, 2022). The implementation runs in
  ``O(n + m)`` time.
- [`temporal_pagerank`](@ref)`(g, α, β, γ)`: temporal PageRank (Rozenshtein and
  Gionis, 2016).

## Local statistics

- [`edge_burstiness`](@ref) and [`node_burstiness`](@ref): ``(σ - μ) / (σ + μ)`` of
  the inter-contact times, from ``-1`` (periodic) to ``1`` (bursty).
- [`temporal_clustering_coefficient`](@ref): temporal edges among the
  out-neighbors of a node, normalized by the number of pairs and the length of the
  time window.
- [`topological_overlap`](@ref): mean overlap ``|A ∩ B| / \sqrt{|A||B|}`` of the
  out-neighborhoods at consecutive time stamps.
- [`get_statistics`](@ref): numbers of nodes, edges, time stamps and degree ranges.

## Cores

- [`kcores`](@ref): core numbers of a static graph (via `Graphs.core_number`).
- [`temporal_khcores`](@ref): ``(k, h)``-cores, the cores of the graph of the node
  pairs with at least ``h`` contacts (Wu et al., 2015).
- [`temporal_lkcores`](@ref): the largest ``(L, k)``-lasting core, a set of nodes
  that forms a ``k``-core during ``L`` consecutive time stamps (Hung and Tseng, 2021).

## References

- F. Brunelli, P. Crescenzi, and L. Viennot. *Making temporal betweenness computation faster and restless.* KDD, 2024. [DOI](https://doi.org/10.1145/3637528.3671825), [arXiv](https://arxiv.org/abs/2501.12708)
- T. Zhang, Y. Gao, J. Zhao, L. Chen, L. Jin, Z. Yang, B. Cao, and J. Fan. *Efficient exact and approximate betweenness centrality computation for temporal graphs.* The Web Conference (WWW), 2024. [DOI](https://doi.org/10.1145/3589334.3645438)
- M. Riondato and E. Upfal. *ABRA: approximating betweenness centrality in static and dynamic graphs with Rademacher averages.* ACM Transactions on Knowledge Discovery from Data 12(5), 2018. [DOI](https://doi.org/10.1145/3208351)
- S. Buß, H. Molter, R. Niedermeier, and M. Rymar. *Algorithmic aspects of temporal betweenness.* KDD, 2020. [DOI](https://doi.org/10.1145/3394486.3403259), [arXiv](https://arxiv.org/abs/2006.08668)
- R. Becker, P. Crescenzi, A. Cruciani, and B. Kodric. *Proxying betweenness centrality rankings in temporal networks.* SEA, 2023. [DOI](https://doi.org/10.4230/LIPIcs.SEA.2023.6)
- A. Cruciani. *MANTRA: temporal betweenness centrality approximation through sampling.* ECML PKDD, 2024. [DOI](https://doi.org/10.1007/978-3-031-70341-6_8), [arXiv](https://arxiv.org/abs/2304.08356)
- D. Santoro and I. Sarpe. *ONBRA: rigorous estimation of the temporal betweenness centrality in temporal networks.* The Web Conference (WWW), 2022. [DOI](https://doi.org/10.1145/3485447.3512204), [arXiv](https://arxiv.org/abs/2203.00653)
- L. Oettershagen and P. Mutzel. *Efficient top-k temporal closeness calculation in temporal networks.* ICDM, 2020. [DOI](https://doi.org/10.1109/ICDM50108.2020.00049)
- P. Crescenzi, C. Magnien, and A. Marino. *Finding top-k nodes for temporal closeness in large temporal graphs.* Algorithms 13(9), 2020. [DOI](https://doi.org/10.3390/a13090211)
- F. Béres, R. Pálovics, A. Oláh, and A. A. Benczúr. *Temporal walk based centrality metric for graph streams.* Applied Network Science 3, 2018. [DOI](https://doi.org/10.1007/s41109-018-0080-5)
- L. Oettershagen, P. Mutzel, and N. M. Kriege. *Temporal walk centrality: ranking nodes in evolving networks.* The Web Conference (WWW), 2022. [DOI](https://doi.org/10.1145/3485447.3512210), [arXiv](https://arxiv.org/abs/2202.03706)
- P. Rozenshtein and A. Gionis. *Temporal PageRank.* ECML PKDD, 2016. [DOI](https://doi.org/10.1007/978-3-319-46227-1_42)
- H. Wu, J. Cheng, Y. Lu, Y. Ke, Y. Huang, D. Yan, and H. Wu. *Core decomposition in large temporal graphs.* IEEE International Conference on Big Data, 2015. [DOI](https://doi.org/10.1109/BigData.2015.7363809)
- W.-C. Hung and C.-Y. Tseng. *Maximum (L, K)-lasting cores in temporal social networks.* DASFAA, 2021. [DOI](https://doi.org/10.1007/978-3-030-73216-5_23)
