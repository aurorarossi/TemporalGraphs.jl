# Performance tips

TemporalGraphs.jl follows the [Julia performance tips](https://docs.julialang.org/en/v1/manual/performance-tips/):
concrete parametric types, type-stable kernels specialized for every distance type,
contiguous memory layouts, in-place variants of all hot loops and task-based
parallelism. The test suite checks type stability (`@inferred`, JET.jl) and that the
in-place kernels do not allocate.

## Choosing types

The graph types are parametric in the node id type `V` and the time type `T`. The
defaults, `Int32` node ids and the time type of the input, make an edge with `Int64`
times 24 bytes long. Use `TemporalEdge{Int64,T}` (or `node_type=Int64` when loading)
only for graphs with more than two billion nodes or time-nodes, and floating point
times only when the data is really continuous.

## Reusing buffers

Every distance computation has an in-place form. [`distance_workspace`](@ref) creates
the buffers once, and [`temporal_distances!`](@ref) reuses them; after the buffers
have grown to their final size, no memory is allocated:

```jldoctest
julia> using TemporalGraphs

julia> g = OrderedEdgeList(4, [(1, 2, 1, 1), (2, 3, 2, 1), (3, 4, 3, 1), (1, 4, 9, 1)]);

julia> ws = distance_workspace(g, Fastest());

julia> dist = Vector{distance_eltype(g, Fastest())}(undef, num_nodes(g));

julia> for s in 1:num_nodes(g)
           temporal_distances!(dist, ws, g, s, Fastest())
           # ... use dist ...
       end

julia> dist
4-element Vector{Int64}:
 9223372036854775807
 9223372036854775807
 9223372036854775807
                   0
```

Workspaces reset only the part of the buffers that the previous computation used, so
a query that reaches few nodes is cheap even in a large graph. A workspace must not
be shared by tasks that run at the same time.

## Multithreading

Measures that need one computation per node — [`temporal_closeness`](@ref) of all
nodes, [`temporal_diameter`](@ref), [`temporal_efficiency`](@ref) and
[`temporal_edge_betweenness`](@ref) — run in parallel with
[OhMyThreads.jl](https://github.com/JuliaFolds2/OhMyThreads.jl), with one workspace
per task. Start Julia with several threads to use them:

```
julia --threads=auto
```

## Choosing a representation

- Single-source distances and closeness: the edge stream ([`OrderedEdgeList`](@ref))
  is almost always the fastest, since it scans one contiguous vector.
- Paths and top-k closeness: [`IncidentLists`](@ref); the search stops as soon as the
  target (or all nodes) is settled.
- Many earliest arrival or fastest path queries on the same graph: the
  [`TRSGraph`](@ref) is an option, but its construction costs about as much as a
  few dozen edge stream queries.

## Latency

A precompilation workload (PrecompileTools.jl) compiles the main code paths for
integer and floating point times when the package is installed, so the first call
in a new session is fast.
