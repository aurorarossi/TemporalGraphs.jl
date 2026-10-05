# Graph representations

All representations are parametric in the node id type `V` and the time type `T`
and store their data in a few contiguous vectors.

| Type | Layout | Best for |
|---|---|---|
| [`OrderedEdgeList`](@ref) | edges sorted by time | distances, closeness, walk centralities, statistics |
| [`IncidentLists`](@ref) | outgoing edges of each node, sorted by time (CSR) | paths, top-k closeness, clustering, overlap |
| [`TRSGraph`](@ref) | time-respecting static DAG of time-nodes (CSR) | repeated reachability and fastest path queries |
| [`DirectedLineGraph`](@ref) | edges as nodes (CSR) | betweenness |

Convert between them with [`to_incident_lists`](@ref), [`to_ordered_edge_list`](@ref),
[`to_trs_graph`](@ref) and [`to_directed_line_graph`](@ref). Distances work on the
first three representations with the same functions:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 4, 1), (1, 3, 6, 1)])) # hide
```

```jldoctest
julia> using TemporalGraphs

julia> g = OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 4, 1), (1, 3, 6, 1)]);

julia> il = to_incident_lists(g); trs = to_trs_graph(g);

julia> minimum_durations(g, 1) == minimum_durations(il, 1) == minimum_durations(trs, 1)
true

julia> out_edges(il, 1)
2-element view(::Vector{TemporalEdge{Int32, Int64}}, 1:2) with eltype TemporalEdge{Int32, Int64}:
 (1 2 1 1)
 (1 3 6 1)
```

## Edge stream

[`OrderedEdgeList`](@ref) keeps the edges sorted by time in one vector: by
`(t, tt, u, v)`, or in the given order if they were already sorted by time (with the
edges of transition time 0 first at every time stamp). The nodes are `1:n` at every
time (see [Nodes over time](@ref)). It also stores the ids that the nodes had in the input file ([`original_id`](@ref),
[`node_map`](@ref)). This is the representation returned by the loaders.

## Incident lists

[`IncidentLists`](@ref) groups the edges by tail in compressed sparse row format:
`out_edges(g, u)` is a view of a contiguous block of the edge vector, sorted by time,
so the edges leaving `u` after a given time are found by binary search.

## Time-respecting static graph

A [`TRSGraph`](@ref) has one *time-node* `(v, τ)` for each time `τ` at which `v` has
an outgoing edge, plus one for the latest arrival at `v`. Consecutive time-nodes of
the same node are linked by waiting edges; a temporal edge `(u, v, t, tt)` links
`(u, t)` to the first time-node of `v` not before `t + tt`. Use
[`num_trs_nodes`](@ref), [`trs_node`](@ref), [`trs_time`](@ref) and
[`trs_neighbors`](@ref) to inspect it.

## Directed line graph

The [`DirectedLineGraph`](@ref) has one node per temporal edge and an edge from `e`
to `f` whenever `f` can follow `e` in a temporal path. It can be much larger than the
temporal graph.

## Snapshot graphs

A discrete-time temporal graph is often given, or best analyzed, as a sequence of
static *snapshot graphs*, one per time step.

- [`snapshots`](@ref) cuts the time interval into steps of length `resolution`. Each
  step becomes a Graphs.jl graph of the edges whose time stamp falls in it. By
  default the step is the resolution of the time stamps, so every time stamp has its
  own snapshot.
- [`OrderedEdgeList`](@ref)`(graphs, times)` converts a sequence of snapshots back
  into a temporal graph. Undirected snapshots give both directions. With
  `transition_time = 1` and consecutive times the paths are strict, using at most
  one edge per snapshot; with `transition_time = 0` they are non-strict.
- [`aggregate_time`](@ref) reduces the time resolution, for example from seconds to
  hours, by moving every edge to the start of its step.
- [`static_graph`](@ref)`(g, ti)` returns the window graph of the edges in a time
  window.

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 0, 1), (2, 3, 0, 1), (1, 3, 5, 1), (2, 1, 7, 1)])) # hide
```

```jldoctest snapshots
julia> using TemporalGraphs, Graphs

julia> g = OrderedEdgeList(3, [(1, 2, 0, 1), (2, 3, 0, 1), (1, 3, 5, 1), (2, 1, 7, 1)]);

julia> S = snapshots(g; resolution = 5);

julia> S.times, ne.(S.graphs)
([0, 5], [2, 2])

julia> edges(OrderedEdgeList(S.graphs, S.times; transition_time = 5))
4-element Vector{TemporalEdge{Int32, Int64}}:
 (1 2 0 5)
 (2 3 0 5)
 (1 3 5 5)
 (2 1 5 5)

julia> edges(aggregate_time(g, 5))   # the edge at time 7 moves to the step starting at 5
4-element Vector{TemporalEdge{Int32, Int64}}:
 (1 2 0 1)
 (2 3 0 1)
 (1 3 5 1)
 (2 1 5 1)
```

The snapshot-based measures of the package ([`persistent_components`](@ref),
[`interval_connected_components`](@ref), [`window_components`](@ref) and the
snapshot shufflings of [`randomize`](@ref)) use the same time steps.

## Transformations

- [`normalize_graph`](@ref): remove identical edges (and self-loops).
- [`make_undirected`](@ref): add the reverse of every edge.
- [`scale_timestamps`](@ref), [`unit_transition_times`](@ref).
- `reverse(g, ti)`: the temporal transpose, which turns paths from `u` to `w` into
  paths from `w` to `u` with the same duration.
- [`to_aggregated_edge_list`](@ref): the static graph with contact counts;
  [`static_graph`](@ref) returns it as a Graphs.jl `SimpleDiGraph`.

## Graphs.jl

The package extends `Graphs.edges`, `Graphs.nv` and `Graphs.ne` for its graph types,
so `using Graphs, TemporalGraphs` does not cause name clashes.
