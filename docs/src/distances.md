# Temporal paths and distances

## Temporal paths

A **temporal path** (or journey) from `s` to `z` is a sequence of temporal edges
``e_1 = (s, v_1, t_1, λ_1), e_2 = (v_1, v_2, t_2, λ_2), \dots, e_k = (v_{k-1}, z, t_k, λ_k)``
in which every edge starts after the previous one has arrived:

```math
t_{i+1} \geq t_i + λ_i \quad \text{for } 1 \leq i < k.
```

It departs at ``t_1``, arrives at ``t_k + λ_k``, has **duration**
``t_k + λ_k - t_1``, **length** ``k`` and **transition time** ``\sum_i λ_i``.

Waiting at a node is allowed for any amount of time. Edges with transition time 0
are supported: several of them can be traversed at the same time stamp.

## Distance types

The distance type is chosen with a singleton [`DistanceType`](@ref), which lets the
compiler generate a specialized algorithm for each of them:

| Distance type | Optimal path | Value for node `v` | Function |
|---|---|---|---|
| [`EarliestArrival`](@ref) | arrives first | arrival time | [`earliest_arrival_times`](@ref) |
| [`LatestDeparture`](@ref) | leaves the source last | departure time | [`latest_departure_times`](@ref) |
| [`Fastest`](@ref) | shortest duration | duration | [`minimum_durations`](@ref) |
| [`MinimumTransitionTimes`](@ref) | smallest sum of transition times ("shortest") | sum of transition times | [`minimum_transition_times`](@ref) |
| [`MinimumHops`](@ref) | fewest edges | number of edges | [`minimum_hops`](@ref) |

For the source itself the value is 0 (for latest departure the end of the time
interval). Nodes that cannot be reached get `typemax` of the element type (for
latest departure `typemin`), which is [`INF`](@ref) for `Int64` times. The results
are [`TemporalDistances`](@ref) vectors, which print these values as `∞` (`-∞`) and
otherwise behave as plain vectors. Minimum hops are always counted as `Int`; all
other values have the time type of the graph.

## Time windows

Every algorithm takes an optional time interval `(a, b)`, by default the interval
spanned by the graph. Only edges that depart and arrive inside the window are used,
i.e. edges with `a ≤ t` and `t + tt ≤ b`:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 4, 1), (1, 3, 6, 1)])) # hide
```

```jldoctest
julia> using TemporalGraphs

julia> g = OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 4, 1), (1, 3, 6, 1)]);

julia> earliest_arrival_times(g, 1)
3-element TemporalDistances{Int64}:
 0
 2
 5

julia> earliest_arrival_times(g, 1, (0, 4))       # (2 3 4 1) arrives at 5 > 4
3-element TemporalDistances{Int64}:
 0
 2
 ∞
```

## Algorithms

| Representation | Earliest arrival | Fastest, latest departure, shortest, minimum hops |
|---|---|---|
| [`OrderedEdgeList`](@ref) | one scan over the edges | one scan, Pareto fronts of labels per node |
| [`IncidentLists`](@ref) | Dijkstra's algorithm | label setting (Dijkstra on labels) |
| [`TRSGraph`](@ref) | depth-first search | depth-first search (fastest), Dijkstra (shortest, hops) |

The edge stream algorithms (Wu et al., 2016) are usually the fastest. They keep for
every node the non-dominated pairs *(arrival time, cost)* in a sorted vector; since
the edges are scanned in chronological order, labels older than the one used by the
last query are discarded. Edges with transition time 0 at the same time stamp are
scanned repeatedly until nothing changes, so they are handled correctly (a single
pass, as in TGLib, misses paths through several such edges).

The label-setting algorithms on incident lists keep the same Pareto fronts, search
them with binary search and stop as soon as every node is settled. They also return
the optimal paths, see [`minimum_duration_path`](@ref) and the other path functions.

## Paths

| Distance type | Function |
|---|---|
| earliest arrival | [`earliest_arrival_path`](@ref) |
| fastest | [`minimum_duration_path`](@ref) |
| shortest | [`minimum_transition_time_path`](@ref) |
| minimum hops | [`minimum_hops_path`](@ref) |

Each returns the vector of edges of one optimal path from `s` to `target`, or an
empty vector if there is none.

### Refined optimality criteria

When there are several optimal paths, as for counting them in a betweenness
centrality, four more criteria select some of them. For a path ``P = e_1, \dots, e_k``
from `s` to `z`, let ``\mathrm{dep}(P) = t_1``, ``\mathrm{arr}(P) = t_k + λ_k`` and
``|P| = k``:

| Criterion | Optimal paths from `s` to `z` |
|---|---|
| [`ShortestForemost`](@ref) | among the paths with the earliest ``\mathrm{arr}(P)``, those with the fewest edges |
| [`ShortestFastest`](@ref) | among the paths with the smallest ``\mathrm{arr}(P) - \mathrm{dep}(P)``, those with the fewest edges |
| [`ShortestLatest`](@ref) | among the paths with the latest ``\mathrm{dep}(P)``, those with the fewest edges |
| [`PrefixForemost`](@ref) | the foremost paths whose every prefix ``e_1, \dots, e_i`` is foremost too, i.e. reaches ``v_i`` at its earliest arrival time (Buß et al., 2020) |

The distance values are those of [`EarliestArrival`](@ref), [`Fastest`](@ref) and
[`LatestDeparture`](@ref): the criteria change only *which* optimal paths count. For
this reason they are accepted only by the betweenness functions
([`temporal_betweenness`](@ref), [`temporal_ego_betweenness`](@ref) and the
approximations, see [Centralities and statistics](centralities.md)), not by
[`temporal_distances`](@ref).

In the following graph there are three foremost paths from 1 to 3, all arriving at
time 6: ``P_1 = (1\,2\,1\,1), (2\,3\,5\,1)``, ``P_2 = (1\,2\,3\,1), (2\,3\,5\,1)`` and
``P_3 = (1\,3\,5\,1)``. ``P_2`` is not prefix foremost, since its prefix ``(1\,2\,3\,1)``
reaches 2 at time 4 instead of 2, and only the direct ``P_3`` is shortest foremost:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 1, 1), (1, 2, 3, 1), (2, 3, 5, 1), (1, 3, 5, 1)])) # hide
```

The betweenness of node 2 is the fraction of the optimal paths from 1 to 3 that pass
through it: 2/3 (foremost), 1/2 (prefix foremost) and 0 (shortest foremost):

```jldoctest
julia> using TemporalGraphs

julia> h = OrderedEdgeList(3, [(1, 2, 1, 1), (1, 2, 3, 1), (2, 3, 5, 1), (1, 3, 5, 1)]);

julia> round(temporal_betweenness(h, EarliestArrival())[2]; digits=3)
0.667

julia> temporal_betweenness(h, PrefixForemost())[2]
0.5

julia> temporal_betweenness(h, ShortestForemost())[2]
0.0
```

## Reachability

[`number_of_reachable_nodes`](@ref) counts the nodes reachable from a source, the
source included.

## References

- S. Buß, H. Molter, R. Niedermeier, and M. Rymar. *Algorithmic aspects of temporal betweenness.* KDD, 2020. [DOI](https://doi.org/10.1145/3394486.3403259), [arXiv](https://arxiv.org/abs/2006.08668)
- H. Wu, J. Cheng, Y. Ke, S. Huang, Y. Huang, and H. Wu. *Efficient algorithms for temporal path computation.* IEEE Transactions on Knowledge and Data Engineering 28(11), 2016. [DOI](https://doi.org/10.1109/TKDE.2016.2594065)
- L. Oettershagen and P. Mutzel. *TGLib: an open-source library for temporal graph analysis.* IEEE International Conference on Data Mining Workshops (ICDMW), 2022. [DOI](https://doi.org/10.1109/ICDMW58026.2022.00160), [arXiv](https://arxiv.org/abs/2209.12587)
