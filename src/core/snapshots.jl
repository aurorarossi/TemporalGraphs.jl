# Snapshot representations: a temporal graph as a sequence of static graphs, one per
# time step, and back; reducing the time resolution.

# Time grid (a, Δ, number of steps) of ti; Δ === nothing means one step per distinct
# time stamp (floating point times without resolution).
function _snapshot_grid(g::OrderedEdgeList{V,T}, ti, resolution) where {V,T}
    a, b = _interval(T, ti)
    if resolution === nothing
        _, _, Δ = _time_grid(g, nothing)
    else
        Δ = convert(T, resolution)
        Δ > 0 || throw(ArgumentError("the resolution must be positive"))
    end
    return a, b, Δ
end

"""
    snapshots(g::OrderedEdgeList, ti = nothing; resolution = nothing, directed = true)

The *snapshot graphs* of `g`: the time interval `ti = (a, b)` (by default from the start
of the time interval of `g` to its last time stamp) is cut into steps
`[a, a + Δ), [a + Δ, a + 2Δ), …` of length `Δ = resolution`, and the snapshot of a
step is the static graph (a `Graphs.SimpleDiGraph`, or `SimpleGraph` if
`directed = false`) of the edges whose time stamp falls in it. By default `Δ` is the
resolution of the time stamps (the gcd of the times minus `a`), so that every time
stamp has its own snapshot and empty steps give empty snapshots; a larger `Δ`
aggregates the edges of every step (see also [`aggregate_time`](@ref)). For floating
point times without `resolution` there is one snapshot per distinct time stamp.

Returns `(times, graphs)`, where `times[k]` is the start of step `k`. Multiple edges
of a step become one edge of the snapshot. [`OrderedEdgeList`](@ref)`(graphs, times)`
converts a sequence of snapshots back into a temporal graph.
"""
function snapshots(g::OrderedEdgeList{V,T}, ti=nothing; resolution=nothing, directed::Bool=true) where {V,T}
    if ti === nothing  # from the start of the graph to its last time stamp
        a = time_interval(g)[1]
        ti = (a, isempty(edges(g)) ? a : max(a, last(edges(g)).t))
    end
    a, b, Δ = _snapshot_grid(g, ti, resolution)
    n = num_nodes(g)
    es = [e for e in edges(g) if a <= e.t <= b]
    if Δ === nothing
        times = unique!([e.t for e in es])  # the edges are sorted by time
        step = Dict(t => k for (k, t) in enumerate(times))
        index = e -> step[e.t]
    else
        K = fld(b - a, Δ) + 1
        times = [a + Δ * k for k in 0:K-1]
        index = e -> fld(e.t - a, Δ) + 1
    end
    graphs = [directed ? Graphs.SimpleDiGraph(n) : Graphs.SimpleGraph(n) for _ in times]
    for e in es
        Graphs.add_edge!(graphs[index(e)], Int(e.u), Int(e.v))
    end
    return (times=times, graphs=graphs)
end

"""
    OrderedEdgeList(graphs::AbstractVector{<:Graphs.AbstractGraph}, times = 1:length(graphs);
                    transition_time = 1, node_type = Int32)

The temporal graph of a sequence of snapshot graphs: every edge `(u, v)` of
`graphs[k]` becomes the temporal edge `(u, v, times[k], transition_time)`, in both
directions for undirected graphs. The default transition time 1 with consecutive
integer times gives *strict* paths, which use at most one edge per snapshot; use
`transition_time = 0` for non-strict paths through several edges of a snapshot. All
graphs must have the same number of nodes. The time interval is
`(times[1], times[end] + transition_time)`.
"""
function OrderedEdgeList(graphs::AbstractVector{<:Graphs.AbstractGraph}, times::AbstractVector{<:Real}=1:length(graphs);
                         transition_time::Real=1, node_type::Type{<:Integer}=Int32)
    length(graphs) == length(times) || throw(ArgumentError("one time is needed per snapshot"))
    isempty(graphs) && return OrderedEdgeList(0, TemporalEdge{node_type,eltype(times)}[])
    issorted(times; lt=<=) || throw(ArgumentError("the times must be increasing"))
    n = Graphs.nv(first(graphs))
    all(h -> Graphs.nv(h) == n, graphs) || throw(ArgumentError("all snapshots must have the same nodes"))
    T = promote_type(eltype(times), typeof(transition_time))
    tt = convert(T, transition_time)
    es = TemporalEdge{node_type,T}[]
    for (h, t) in zip(graphs, times), e in Graphs.edges(h)
        u, v = Graphs.src(e), Graphs.dst(e)
        push!(es, TemporalEdge{node_type,T}(u, v, t, tt))
        (!Graphs.is_directed(h) && u != v) && push!(es, TemporalEdge{node_type,T}(v, u, t, tt))
    end
    return OrderedEdgeList(n, sort!(es), (convert(T, times[1]), convert(T, times[end]) + tt))
end

"""
    aggregate_time(g::OrderedEdgeList, Δ, ti = time_interval(g); transition_time = nothing,
                   multiedges = false)

Reduce the time resolution of `g` to `Δ`: every edge inside `ti = (a, b)` (by its time
stamp) is moved to the start `a + kΔ` of its step `[a + kΔ, a + (k+1)Δ)`, for example
from seconds to hours. Edges that become identical are merged unless `multiedges`.
The transition times are kept, or all set to `transition_time` (`Δ` makes the paths
strict at the new resolution, `0` non-strict). The result has the same nodes; its time
interval covers the steps.
"""
function aggregate_time(g::OrderedEdgeList{V,T}, Δ::Real, ti=time_interval(g); transition_time=nothing,
                        multiedges::Bool=false) where {V,T}
    a, b = _interval(T, ti)
    Δ > 0 || throw(ArgumentError("Δ must be positive"))
    D = convert(T, Δ)
    es = TemporalEdge{V,T}[]
    for e in edges(g)
        a <= e.t <= b || continue
        t = a + D * fld(e.t - a, D)
        push!(es, TemporalEdge{V,T}(e.u, e.v, t, transition_time === nothing ? e.tt : convert(T, transition_time)))
    end
    sort!(es)
    multiedges || unique!(es)
    hi = isempty(es) ? a : max(maximum(e.t + e.tt for e in es), a + D * fld(b - a, D) + D)
    return OrderedEdgeList(num_nodes(g), es, (a, hi); original_ids=g.original_ids)
end
