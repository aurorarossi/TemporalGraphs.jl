# Connectivity, flows and spanners

This page covers classical problems of temporal graph theory. In temporal graphs,
several notions that are easy for static graphs behave differently:

- reachability is not transitive;
- connected components can overlap;
- Menger's theorem fails for vertex-disjoint paths;
- spanning trees and spanners can be hard to compute.

For each function, the docstring gives the complexity of the problem it solves. Some
problems are NP-hard; for those the package provides exact exponential-time
algorithms that are practical on small and medium instances, or approximation
algorithms. Every function was validated against brute-force enumeration and against
the examples of the original papers.

!!! info "Conventions"
    An edge `(u, v, t, tt)` can be followed by an edge leaving `v` at time `t + tt` or
    later. Positive transition times give *strict* paths; transition times 0 give
    *non-strict* paths, which can use several edges with the same time stamp.

## Reachability and temporal components

[`temporal_reachability`](@ref) computes the reachability matrix between all pairs of
nodes. It processes 64 sources at a time with bitsets, in one scan of the edges per
block of sources, so it is much faster than one search per source.
[`is_temporally_connected`](@ref) tests whether every node reaches every other node.

A *temporal connected component* (Bhadra and Ferreira, 2003) is a maximal set of
nodes that reach each other.

- **Overlapping components.** Since reachability is not transitive, components can
  overlap, and finding the largest one is NP-hard.
- **Open components.** These are the maximal cliques of the graph that joins the
  pairs of nodes reaching each other (Costa, Lopes, Marino and Silva, 2024);
  [`temporal_connected_components`](@ref) enumerates them, and
  [`largest_temporal_connected_component`](@ref) finds a largest one with the maximum
  clique algorithm of Tomita and Seki (2003).
- **Closed components.** In a *closed* component the paths must also stay inside the
  set. Even deciding whether a set is a closed component is NP-hard.
  [`largest_temporal_connected_component`](@ref) finds a largest closed component
  exactly, by branching inside the open components.
- **Unilateral components.** Both functions accept `unilateral = true`: then it is
  enough that one node of every pair reaches the other.

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(4, [(1, 3, 1, 0), (3, 2, 2, 0), (2, 4, 3, 0), (4, 1, 4, 0)])) # hide
```

```jldoctest theory
julia> using TemporalGraphs

julia> g = OrderedEdgeList(4, [(1, 3, 1, 0), (3, 2, 2, 0), (2, 4, 3, 0), (4, 1, 4, 0)]);

julia> temporal_connected_components(g)   # overlapping components
3-element Vector{Vector{Int64}}:
 [1, 2]
 [1, 3]
 [1, 4]

julia> largest_temporal_connected_component(g; closed = true)   # pairs need other nodes
1-element Vector{Int64}:
 1
```

## More temporal components

The [Temporal Graph Wiki](https://temporalgraph.notion.site/) lists several other
component notions. All of them are
computed in polynomial time, except the Δ-temporal connected components, which can
overlap like the temporal connected components.

- **Source and sink components.** [`source_component`](@ref)`(g, S)` is the largest
  set of nodes reached by every node of `S`; [`sink_component`](@ref) is the largest
  set of nodes that reach every node of `S`. With `closed = true` the paths must stay
  inside the set. Unlike for temporal connected components, the largest set is
  unique in both cases. The closed one is the greatest fixed point of restricting the
  set to the subgraph it induces.
- **Window components.** [`window_components`](@ref) returns the connected components
  of the static graph of the edges in a time window, ignoring their order.
- **Persistent components** (Vernet, Pigné and Sanlaville, 2023).
  [`persistent_components`](@ref) finds the sets that are a connected component of
  every snapshot during an interval of consecutive time steps, with the interval as
  long as possible. It takes one sweep over the snapshots.
- **T-interval connected components** (after Kuhn, Lynch and Oshman, 2010, and
  Casteigts, Klasing, Neggaz and Peters, 2015, for the whole graph).
  [`interval_connected_components`](@ref) finds the maximal sets of nodes that are
  connected, in every window of `T` consecutive time steps, by the edges present in
  all the snapshots of that window. These sets partition the nodes, and are found by
  refining a partition until it no longer changes.
- **Δ-temporal connected components.**
  [`delta_temporal_connected_components`](@ref)`(g, Δ)` finds the maximal sets whose
  nodes reach each other inside every time window of length `Δ`, e.g. groups that stay
  in touch every week and not only once over the whole interval. Only the windows that
  start right after a time stamp need to be checked; the components are then
  enumerated as the temporal connected components.
- **Stream components** (Latapy, Viard and Magnien, 2018). [`stream_components`](@ref)
  sees `g` as a stream graph: a node is present while it has an edge, and the elements
  of the components are temporal nodes, a node together with a presence interval.
  Presence intervals at most `gap` apart are merged. Time order and directions are
  ignored, so a component is a group that keeps interacting without long
  interruptions.

The time steps are the grid of the time interval with step the
[resolution](@ref "Snapshot graphs") of the time stamps, here 1:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(4, [(1, 2, 0, 0), (1, 2, 1, 0), (2, 3, 1, 0), (1, 2, 2, 0), (2, 3, 2, 0), (3, 4, 5, 0)])) # hide
```

```jldoctest theory
julia> c = OrderedEdgeList(4, [(1, 2, 0, 0), (1, 2, 1, 0), (2, 3, 1, 0),
                               (1, 2, 2, 0), (2, 3, 2, 0), (3, 4, 5, 0)]);

julia> persistent_components(c)
3-element Vector{@NamedTuple{nodes::Vector{Int64}, interval::Tuple{Int64, Int64}, length::Int64}}:
 (nodes = [1, 2, 3], interval = (1, 2), length = 2)
 (nodes = [1, 2], interval = (0, 0), length = 1)
 (nodes = [3, 4], interval = (5, 5), length = 1)

julia> interval_connected_components(c, 2, (0, 2))   # windows {0, 1} and {1, 2}
3-element Vector{Vector{Int64}}:
 [1, 2]
 [3]
 [4]
```

In `d` the nodes 1 and 2 meet every 3 time units, while 3 only meets 2 at the start;
in `s` the nodes 1 and 2 meet twice, at times 0 and 5, so each of them has two
temporal nodes:

```jldoctest theory
julia> d = make_undirected(OrderedEdgeList(3, [(1, 2, 1, 0), (2, 3, 1, 0),
                                               (1, 2, 4, 0), (1, 2, 7, 0)]));

julia> temporal_connected_components(d)
1-element Vector{Vector{Int64}}:
 [1, 2, 3]

julia> delta_temporal_connected_components(d, 3)
2-element Vector{Vector{Int64}}:
 [1, 2]
 [3]

julia> s = OrderedEdgeList(3, [(1, 2, 0, 1), (2, 3, 1, 1), (1, 2, 5, 1)]);

julia> stream_components(s)
2-element Vector{Vector{Tuple{Int64, Tuple{Int64, Int64}}}}:
 [(1, (0, 1)), (2, (0, 2)), (3, (1, 2))]
 [(1, (5, 6)), (2, (5, 6))]
```

## Dissemination times

If every node forwards what it knows on every edge it can use (*flooding*, see
Clementi et al., 2010), then:

- information starting at `s` at the beginning of the time interval reaches everybody
  after the [`temporal_flooding_time`](@ref) of `s`;
- everybody knows everything after the [`temporal_gossip_time`](@ref), the largest
  flooding time.

Both are computed from [`earliest_arrival_matrix`](@ref), the earliest arrival times
between all pairs of nodes. It uses the bitset scan of the reachability matrix and
records when every source first reaches every node, in `O(m n / 64 + n²)` time.

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1), (3, 1, 4, 1), (1, 2, 6, 1)])) # hide
```

```jldoctest theory
julia> g = OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1), (3, 1, 4, 1), (1, 2, 6, 1)]);

julia> earliest_arrival_matrix(g)
3×3 TemporalDistances{Int64}:
 0  2  3
 5  0  3
 5  7  0

julia> temporal_flooding_times(g), temporal_gossip_time(g)   # the interval starts at time 1
([2, 4, 6], 6)
```

## Seed selection

Which nodes should be informed, or infected, at the start of the time interval so
that the information reaches everybody, or as many nodes as possible, if every
contact transmits? Both are covering problems on the rows of the reachability matrix,
and both are NP-hard.

- **Reachability dominating sets.** [`reachability_dominating_set`](@ref) returns a
  small set of nodes from which every node is reached (a TaRDiS, Kutner and
  Larios-Jones, 2026). The greedy algorithm is at most `ln n + 1` times larger than
  optimal; `exact = true` finds a smallest set by branch and bound.
- **Maximum reach.** [`max_reach_seeds`](@ref)`(g, k)` returns `k` seeds that reach
  the most nodes: influence maximization (Kempe, Kleinberg and Tardos, 2003) when every
  contact transmits. The greedy algorithm reaches at least `1 - 1/e` of the optimum;
  `exact = true` finds an optimal set.

```jldoctest theory
julia> p = OrderedEdgeList(5, [(1, 2, 1, 1), (2, 3, 2, 1), (4, 5, 1, 1), (5, 1, 3, 1)]);

julia> reachability_dominating_set(p)   # 1 reaches 2 and 3, 4 reaches 5 and 1
2-element Vector{Int64}:
 1
 4

julia> max_reach_seeds(p, 1)
(seeds = [1], reached = 3)
```

## Static expansions

Static expansions turn temporal walks into directed walks of an ordinary static
graph. They are the standard tool for reducing temporal questions (reachability,
optimal paths, flows, cuts) to static ones, which can then be answered with Graphs.jl.

- **Vertex expansion.** [`vertex_expansion`](@ref) is the time-expanded graph. It has
  a node `(v, τ)` per node and time, waiting arcs between consecutive times of a
  node, and an arc `(u, t) → (v, t + tt)` per temporal edge. By default only the
  times at which the node has an event are used.
- **Edge expansion.** [`edge_expansion`](@ref) is the temporal line graph, with one
  node per temporal edge.
  - With a waiting constraint `β`, its walks are exactly the restless walks.
  - With `hubs = true`, one hub per departure time replaces the `Θ(m²)` arcs between
    consecutive edges by `O(m)` arcs.

With the graph `g` of [Dissemination times](@ref):

```jldoctest theory
julia> using Graphs

julia> x = vertex_expansion(g);

julia> x.nodes
7-element Vector{Tuple{Int64, Int64}}:
 (1, 1)
 (1, 5)
 (1, 6)
 (2, 2)
 (2, 7)
 (3, 3)
 (3, 4)

julia> has_path(x.graph, 6, 5)   # (3, 3) reaches (2, 7): 3 → 1 → 2
true
```

## Waiting constraints: restless walks and paths

Some processes cannot wait indefinitely at a node: a packet in a network, or an
infection whose carrier recovers. A *β-restless* walk (Casteigts, Himmel, Molter and
Zschoche, 2021) waits at most `β` at every intermediate node.

- **Walks.** Optimal restless walks can be computed in polynomial time.
  [`temporal_distances`](@ref) with `β` gives earliest arrival, latest departure,
  fastest, minimum hop and minimum transition time walks, in `O(m log m)` time per
  source.
- **Paths.** A restless *path* visits every node at most once, and finding one is
  NP-hard for every `β ≥ 1`. [`restless_path`](@ref) is an exact search, pruned with
  the polynomial walk computation.

In the example below, the only 1-restless walk from 1 to 5 goes around the cycle
2 → 3 → 4 → 2, so it is not a path:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(5, [(1, 2, 1, 0), (2, 3, 2, 0), (3, 4, 3, 0), (4, 2, 4, 0), (2, 5, 5, 0)])) # hide
```

```jldoctest theory
julia> g = OrderedEdgeList(5, [(1, 2, 1, 0), (2, 3, 2, 0), (3, 4, 3, 0),
                               (4, 2, 4, 0), (2, 5, 5, 0)]);

julia> temporal_distances(g, 1, EarliestArrival(); β = 1)
5-element TemporalDistances{Int64}:
 0
 1
 2
 3
 5

julia> restless_path(g, 1, 5, 1) === nothing
true

julia> restless_path(g, 1, 5, 4)   # waiting 4 at node 2 is allowed
2-element Vector{TemporalEdge{Int32, Int64}}:
 (1 2 1 0)
 (2 5 5 0)
```

## Flows, cuts and disjoint paths

[`temporal_max_flow`](@ref) computes the maximum amount that can be sent from `s` to
`z` (Akrida, Czyzowicz, Gąsieniec, Kuszner and Spirakis, 2019). Each temporal edge
carries at most its capacity at its time, and nodes store flow between edges. It is
solved as a static maximum flow (Dinic's algorithm) in the time-expanded network,
which has one node per time an edge leaves or reaches a node. Optional arguments
model:

- **capacities:** `capacity[i]` for the `i`-th edge of `edges(g)` (1 by default);
- **bounded storage:** the `buffers` of the nodes;
- **no storage beyond a waiting time:** `β`, with which the flow follows restless
  walks.

**Menger's theorem.** The disjointness notion matters:

| Disjoint paths | Min-max theorem | Complexity | Function |
|:---|:---|:---|:---|
| edge-disjoint | = minimum temporal edge cut (Berman, 1996) | polynomial | [`temporal_edge_disjoint_paths`](@ref), [`temporal_min_cut`](@ref) |
| out-disjoint (never leaving a node at the same time) | = minimum set of node departure times (Mertzios, Michail and Spirakis, 2019) | polynomial | [`temporal_out_disjoint_paths`](@ref) |
| vertex-disjoint | fails (Kempe, Kleinberg and Kumar, 2002) | NP-hard | — |
| vertex separator | — | NP-hard (Zschoche, Fluschnik, Molter and Niedermeier, 2020) | [`temporal_vertex_separator`](@ref) |

The example of Kempe, Kleinberg and Kumar has no two vertex-disjoint paths from 1 to
5, but a minimum vertex separator needs two nodes:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(5, [(1, 2, 1, 1), (2, 1, 1, 1), (2, 3, 2, 1), (3, 2, 2, 1), (3, 5, 3, 1), (5, 3, 3, 1), (2, 4, 4, 1), (4, 2, 4, 1), (1, 3, 5, 1), (3, 1, 5, 1), (3, 4, 6, 1), (4, 3, 6, 1), (4, 5, 7, 1), (5, 4, 7, 1)])) # hide
```

```jldoctest theory
julia> und(es) = reduce(vcat, [[(u, v, t, tt), (v, u, t, tt)] for (u, v, t, tt) in es]);

julia> g = OrderedEdgeList(5, und([(1, 2, 1, 1), (2, 3, 2, 1), (3, 5, 3, 1), (2, 4, 4, 1),
                                   (1, 3, 5, 1), (3, 4, 6, 1), (4, 5, 7, 1)]));

julia> temporal_vertex_separator(g, 1, 5)
2-element Vector{Int64}:
 2
 3

julia> temporal_edge_disjoint_paths(g, 1, 5)
2-element Vector{Vector{TemporalEdge{Int32, Int64}}}:
 [(1 3 5 1), (3 4 6 1), (4 5 7 1)]
 [(1 2 1 1), (2 3 2 1), (3 5 3 1)]
```

The separator is exact. It branches on the inner nodes of a minimum hop path, which
is fixed-parameter tractable in the size of the separator and the number of time
stamps for strict paths.

## Spanning trees and spanners

A *temporal spanning tree* rooted at `r` (Huang, Fu and Liu, 2015) contains one
incoming edge for each node reachable from `r`, such that every tree path is a
temporal path. There are two versions:

- [`earliest_arrival_tree`](@ref) (*MST_a*) reaches every node at its earliest arrival
  time, in linear time. It uses the stack-based algorithm of the paper, which is also
  correct for transition times 0.
- [`minimum_weight_spanning_tree`](@ref) (*MST_w*) minimizes the total weight. This is
  NP-hard, so it is approximated by the paper's reduction to the Directed Steiner Tree
  problem, solved with the greedy algorithm of Charikar et al.

With the default `level = 2`, on random graphs the approximation was within 0.5% of
the optimum on average (computed by brute force), and never worse than 1.5 times it.

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(6, [(1, 2, 1, 2), (1, 3, 1, 4), (1, 3, 3, 3), (1, 2, 4, 1), (2, 4, 4, 2), (2, 5, 5, 3), (3, 6, 6, 2), (3, 5, 7, 2)]); order=[1, 3, 2, 4, 5, 6]) # hide
```

```jldoctest theory
julia> g = OrderedEdgeList(6, [(1, 2, 1, 2), (1, 3, 1, 4), (1, 3, 3, 3), (1, 2, 4, 1),
                               (2, 4, 4, 2), (2, 5, 5, 3), (3, 6, 6, 2), (3, 5, 7, 2)]);

julia> earliest_arrival_tree(g, 1)
5-element Vector{TemporalEdge{Int32, Int64}}:
 (1 2 1 2)
 (1 3 1 4)
 (2 4 4 2)
 (2 5 5 3)
 (3 6 6 2)

julia> minimum_weight_spanning_tree(g, 1).weight   # weights = transition times
11.0
```

A *temporal spanner* is a subset of the edges that keeps the reachability between all
pairs of nodes. Minimum spanners are NP-hard to compute and to approximate (Axiotis
and Fotakis, 2016), and some graphs need ``Ω(n^2)`` edges.

- [`temporal_spanner`](@ref) returns the union of the earliest arrival trees of all
  nodes. With `minimal = true` it then removes edges as long as the reachability does
  not change, so the result is inclusion-minimal.
- [`temporal_clique_spanner`](@ref) handles *temporal cliques*: undirected graphs
  where every pair of nodes has exactly one time label, with transition times 0. These always admit spanners of size ``O(n \log n)``
  (Casteigts, Peters and Schoeters, 2021), and the function computes one with the
  algorithm of Angrick et al. (2024). On random temporal cliques with 128 nodes the
  spanners keep about 480 of the 8128 pairs.

## References

- S. Bhadra and A. Ferreira. *Complexity of connected components in evolving graphs and the computation of multicast trees in wireless networks.* ADHOC-NOW, 2003. [DOI](https://doi.org/10.1007/978-3-540-39611-6_23)
- I. L. Costa, R. Lopes, A. Marino, and A. Silva. *On computing large temporal (unilateral) connected components.* Journal of Computer and System Sciences 144, 2024. [DOI](https://doi.org/10.1016/j.jcss.2024.103548), [arXiv](https://arxiv.org/abs/2302.12068)
- E. Tomita and T. Seki. *An efficient branch-and-bound algorithm for finding a maximum clique.* DMTCS, 2003. [DOI](https://doi.org/10.1007/3-540-45066-1_22)
- A. Casteigts, A.-S. Himmel, H. Molter, and P. Zschoche. *Finding temporal paths under waiting time constraints.* Algorithmica 83(9), 2021. [DOI](https://doi.org/10.1007/s00453-021-00831-w), [arXiv](https://arxiv.org/abs/1909.06437)
- E. C. Akrida, J. Czyzowicz, L. Gąsieniec, Ł. Kuszner, and P. G. Spirakis. *Temporal flows in temporal networks.* Journal of Computer and System Sciences 103, 2019. [DOI](https://doi.org/10.1016/j.jcss.2019.02.003), [arXiv](https://arxiv.org/abs/1606.01091)
- K. A. Berman. *Vulnerability of scheduled networks and a generalization of Menger's theorem.* Networks 28(3), 1996. [DOI](https://doi.org/10.1002/%28SICI%291097-0037%28199610%2928%3A3%3C125%3A%3AAID-NET1%3E3.0.CO%3B2-P)
- D. Kempe, J. Kleinberg, and A. Kumar. *Connectivity and inference problems for temporal networks.* Journal of Computer and System Sciences 64(4), 2002. [DOI](https://doi.org/10.1006/jcss.2002.1829)
- G. B. Mertzios, O. Michail, and P. G. Spirakis. *Temporal network optimization subject to connectivity constraints.* Algorithmica 81(4), 2019. [DOI](https://doi.org/10.1007/s00453-018-0478-6), [arXiv](https://arxiv.org/abs/1502.04382)
- P. Zschoche, T. Fluschnik, H. Molter, and R. Niedermeier. *The complexity of finding small separators in temporal graphs.* Journal of Computer and System Sciences 107, 2020. [DOI](https://doi.org/10.1016/j.jcss.2019.07.006), [arXiv](https://arxiv.org/abs/1711.00963)
- S. Huang, A. W.-C. Fu, and R. Liu. *Minimum spanning trees in temporal graphs.* SIGMOD, 2015. [DOI](https://doi.org/10.1145/2723372.2723717)
- M. Charikar, C. Chekuri, T. Cheung, Z. Dai, A. Goel, S. Guha, and M. Li. *Approximation algorithms for directed Steiner problems.* Journal of Algorithms 33(1), 1999. [DOI](https://doi.org/10.1006/jagm.1999.1042)
- K. Axiotis and D. Fotakis. *On the size and the approximability of minimum temporally connected subgraphs.* ICALP, 2016. [DOI](https://doi.org/10.4230/LIPIcs.ICALP.2016.149), [arXiv](https://arxiv.org/abs/1602.06411)
- A. Casteigts, J. G. Peters, and J. Schoeters. *Temporal cliques admit sparse spanners.* Journal of Computer and System Sciences 121, 2021. [DOI](https://doi.org/10.1016/j.jcss.2021.04.004), [arXiv](https://arxiv.org/abs/1810.00104)
- S. Angrick, B. Bals, T. Friedrich, H. Gawendowicz, N. Hastrich, N. Klodt, P. Lenzner, J. Schmidt, G. Skretas, and A. Wells. *How to reduce temporal cliques to find sparse spanners.* ESA, 2024. [DOI](https://doi.org/10.4230/LIPIcs.ESA.2024.11), [arXiv](https://arxiv.org/abs/2402.13624)
- M. Vernet, Y. Pigné, and É. Sanlaville. *A study of connectivity on dynamic graphs: computing persistent connected components.* 4OR 21, 2023. [DOI](https://doi.org/10.1007/s10288-022-00507-3)
- F. Kuhn, N. Lynch, and R. Oshman. *Distributed computation in dynamic networks.* STOC, 2010. [DOI](https://doi.org/10.1145/1806689.1806760)
- A. Casteigts, R. Klasing, Y. M. Neggaz, and J. G. Peters. *Efficiently testing T-interval connectivity in dynamic graphs.* CIAC, 2015. [DOI](https://doi.org/10.1007/978-3-319-18173-8_6)
- A. E. F. Clementi, C. Macci, A. Monti, F. Pasquale, and R. Silvestri. *Flooding time of edge-Markovian evolving graphs.* SIAM Journal on Discrete Mathematics 24, 2010. [DOI](https://doi.org/10.1137/090756053)
- M. Latapy, T. Viard, and C. Magnien. *Stream graphs and link streams for the modeling of interactions over time.* Social Network Analysis and Mining 8, 2018. [DOI](https://doi.org/10.1007/s13278-018-0537-7), [arXiv](https://arxiv.org/abs/1710.04073)
- D. C. Kutner and L. Larios-Jones. *Temporal reachability dominating sets: contagion in temporal graphs.* Journal of Computer and System Sciences 155, 2026. [DOI](https://doi.org/10.1016/j.jcss.2025.103701), [arXiv](https://arxiv.org/abs/2306.06999)
- D. Kempe, J. Kleinberg, and É. Tardos. *Maximizing the spread of influence through a social network.* KDD, 2003. [DOI](https://doi.org/10.1145/956750.956769)
- M. Döring. *Temporal Graph Wiki.* [temporalgraph.notion.site](https://temporalgraph.notion.site/)
