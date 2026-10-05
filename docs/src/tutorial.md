# Getting started

This page walks through a typical analysis: building a temporal graph, computing
distances and paths, and ranking nodes by centrality.

```@meta
DocTestSetup = quote
    using TemporalGraphs
    g = OrderedEdgeList(4, [(1, 4, 1, 5), (1, 2, 2, 1), (1, 2, 5, 2), (3, 2, 6, 1),
                            (4, 3, 6, 2), (2, 4, 7, 2), (4, 3, 8, 1)])
end
```

## Temporal edges

A [`TemporalEdge`](@ref) `(u, v, t, tt)` is a directed contact from node `u` to node
`v` that starts at time `t` and needs the *transition time* `tt`, so it arrives at
`v` at time `t + tt`. If you do not care about transition times, leave them at the
default 1.

```jldoctest
julia> e = TemporalEdge(1, 2, 10)        # tt = 1
(1 2 10 1)

julia> e.t + e.tt                        # arrival time
11
```

Node ids are 1-based integers. By default they are stored as `Int32`, which keeps an
edge with integer times at 24 bytes; the time type is inferred from `t` and `tt` and
can be any real type, e.g. `Float64` for continuous time:

```jldoctest
julia> TemporalEdge(1, 2, 0.5, 0.25)
(1 2 0.5 0.25)

julia> typeof(ans)
TemporalEdge{Int32, Float64}
```

## Building a temporal graph

The main representation is the [`OrderedEdgeList`](@ref), a chronologically sorted
edge stream. You can build it from a vector of tuples `(u, v, t, tt)` or of
`TemporalEdge`s, for nodes `1:n`:

```jldoctest
julia> g = OrderedEdgeList(3, [(1, 2, 5, 1), (2, 3, 7, 1), (1, 3, 1, 3)])
OrderedEdgeList{Int32, Int64} with 3 nodes, 3 temporal edges, time interval (1, 8)

julia> edges(g)                          # sorted by (t, tt, u, v)
3-element Vector{TemporalEdge{Int32, Int64}}:
 (1 3 1 3)
 (1 2 5 1)
 (2 3 7 1)
```

The time interval defaults to `(earliest departure, latest arrival)`; you can pass
another one as third argument.

### Nodes over time

The nodes `1:n` exist during the whole time interval: only the edges depend on time.
There is no separate notion of a node joining or leaving the network. A node that has
no edges at some time is *isolated* then: it is still a node of the graph, but it
cannot reach or be reached by the others at that time. To model a node that is active
only during a period, give it edges only during that period.

Results over all nodes therefore include the inactive ones: a graph with a node that
never has an edge is not temporally connected ([`is_temporally_connected`](@ref)),
isolated nodes form components of size 1, and every snapshot ([`snapshots`](@ref))
has all `n` nodes. The loaders ([`load_ordered_edge_list`](@ref),
[`load_dataset`](@ref)) number only the nodes that appear in some edge. To study a
subset of the nodes, for example those active in a time window, build a new graph
from the edges between them.

### From files

[`load_ordered_edge_list`](@ref) reads text files with one edge `u v t [tt]` per
line, the format of TGLib and of many public datasets (e.g. the
[SNAP temporal networks](https://snap.stanford.edu/data/#temporal)). Spaces, tabs and
commas are accepted as separators and lines starting with `#` or `%` are skipped:

```julia
g = load_ordered_edge_list("CollegeMsg.txt")                        # integer times
g = load_ordered_edge_list("contacts.csv"; time_type=Float64)       # real times
g = load_ordered_edge_list("contacts.txt"; directed=false)          # both directions
```

The node ids of the file can be arbitrary integers; they are renumbered to `1:n` in
order of first appearance, and [`original_id`](@ref) maps them back:

```jldoctest
julia> path = tempname(); write(path, "100 7 1\n7 42 3\n");

julia> h = load_ordered_edge_list(path);

julia> edges(h)
2-element Vector{TemporalEdge{Int32, Int64}}:
 (1 2 1 1)
 (2 3 3 1)

julia> original_id(h, 3)
42
```

[`save_ordered_edge_list`](@ref) writes a graph in the same format.

### From tables

Any [Tables.jl](https://github.com/JuliaData/Tables.jl) source, such as a
`DataFrame` or a `CSV.File`, can be turned into a temporal graph by naming its
columns:

```jldoctest
julia> tbl = (from=[10, 20, 10], to=[20, 30, 30], time=[1.0, 2.5, 4.0]);

julia> OrderedEdgeList(tbl; u=:from, v=:to, t=:time)
OrderedEdgeList{Int32, Float64} with 3 nodes, 3 temporal edges, time interval (1.0, 5.0)
```

## Distances

Temporal distances from a source node come in five flavors (see
[Temporal paths and distances](distances.md) for the definitions). Each has its own function,
and [`temporal_distances`](@ref) selects one with a [`DistanceType`](@ref):

```jldoctest
julia> earliest_arrival_times(g, 1)
4-element Vector{Int64}:
 0
 3
 8
 6

julia> minimum_hops(g, 1)
4-element Vector{Int64}:
 0
 1
 2
 1

julia> temporal_distances(g, 1, MinimumTransitionTimes())
4-element Vector{Int64}:
 0
 1
 6
 3
```

Unreachable nodes get `typemax` of the element type (`INF` for `Int64`, `Inf` for
`Float64`):

```jldoctest
julia> minimum_durations(g, 3)
4-element Vector{Int64}:
 9223372036854775807
                   1
                   0
                   3
```

All functions accept a time window as last argument; only edges that depart and
arrive inside the window are used:

```jldoctest
julia> earliest_arrival_times(g, 1, (2, 9))
4-element Vector{Int64}:
                   0
                   3
 9223372036854775807
                   9
```

## Paths

The path functions return the edges of one optimal path (empty if the target cannot
be reached):

```jldoctest
julia> earliest_arrival_path(g, 1, 3)
2-element Vector{TemporalEdge{Int32, Int64}}:
 (1 4 1 5)
 (4 3 6 2)

julia> minimum_transition_time_path(g, 1, 4)
2-element Vector{TemporalEdge{Int32, Int64}}:
 (1 2 2 1)
 (2 4 7 2)
```

## Centralities

Closeness, eccentricity, diameter and efficiency take a distance type:

```jldoctest
julia> temporal_closeness(g, Fastest())
4-element Vector{Float64}:
 1.3928571428571428
 0.5
 1.3333333333333333
 1.0

julia> compute_topk_closeness(g, 2, Fastest())       # (node, closeness) pairs
2-element Vector{Tuple{Int64, Float64}}:
 (1, 1.3928571428571428)
 (3, 1.3333333333333333)

julia> temporal_diameter(g, MinimumHops())
2
```

Walk based centralities work on the edge stream directly:

```jldoctest
julia> temporal_katz_centrality(g, 1.0)
4-element Vector{Float64}:
 0.0
 3.0
 8.0
 5.0

julia> temporal_walk_centrality(g, 1.0, 1.0)
4-element Vector{Float64}:
 0.0
 6.0
 0.0
 6.0
```

[`get_statistics`](@ref) summarizes a graph:

```jldoctest
julia> get_statistics(g)
number of nodes: 4
number of edges: 7
number of static edges: 5
number of time stamps: 6
number of transition times: 3
min. time stamp: 1
max. time stamp: 8
min. transition time: 1
max. transition time: 5
min. temporal in-degree: 0
max. temporal in-degree: 3
min. temporal out-degree: 1
max. temporal out-degree: 3
```

See [Centralities and statistics](centralities.md) for all measures.

## Drawing

[`draw_graph`](@ref) draws the static graph, with every edge labelled by its
`(t, tt)`, and [`draw_timelines`](@ref) draws one timeline per node, with every edge
going from `u` at time `t` to `v` at time `t + tt`:

```@example drawing
using TemporalGraphs
g = OrderedEdgeList(4, [(1, 4, 1, 5), (1, 2, 2, 1), (1, 2, 5, 2), (3, 2, 6, 1),
                        (4, 3, 6, 2), (2, 4, 7, 2), (4, 3, 8, 1)])
draw_graph(g)
```

```@example drawing
draw_timelines(g; order=[1, 2, 4, 3])   # the timelines from top to bottom
```

Both return a [`TemporalGraphDrawing`](@ref), shown as an image by notebooks, the VS
Code plot pane and this documentation; `write("graph.svg", draw_graph(g))` saves it.
They draw every edge, so they are meant for small graphs.

## Next steps

- [Graph representations](representations.md): when to convert to incident lists or to the
  time-respecting static graph.
- [Performance tips](performance.md): in-place computations, multithreading and the choice of representation.

```@meta
DocTestSetup = nothing
```
