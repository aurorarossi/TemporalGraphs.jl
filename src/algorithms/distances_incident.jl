# Label-setting algorithms on incident lists.
#
# A label is a path ending at a node, arriving at time `arr` with cost `cost` (as for
# the edge stream: c = -(start time) for fastest and latest departure paths). Labels
# live in a pool of plain vectors (struct of arrays) and are identified by their
# index; label 1 is the source, which can use every edge leaving it inside the
# interval. Every node keeps its non-dominated labels in a Pareto front sorted by
# arrival time, so dominance checks are binary searches. Labels are settled in order
# of increasing distance with a binary heap.

"""
    IncidentListsWorkspace

Reusable buffers for temporal distance computations on [`IncidentLists`](@ref),
created by [`distance_workspace`](@ref).
"""
struct IncidentListsWorkspace{T,C}
    node::Vector{Int}
    arr::Vector{T}
    cost::Vector{C}
    parent::Vector{Int}
    edge::Vector{Int}          # index into edges(g) of the last edge, 0 for the source
    deleted::Vector{Bool}
    fronts::LazyLists{Tuple{T,C,Int}}
    visited::Vector{Bool}
    heap::MinHeap{Tuple{C,Int}}
end

function IncidentListsWorkspace{T,C}(n::Integer) where {T,C}
    return IncidentListsWorkspace{T,C}(Int[], T[], C[], Int[], Int[], Bool[],
                                       LazyLists{Tuple{T,C,Int}}(n), zeros(Bool, n),
                                       MinHeap{Tuple{C,Int}}())
end

struct EarliestArrivalWorkspace{T}
    heap::MinHeap{Tuple{T,Int}}
    done::Vector{Bool}
    pred::Vector{Int}
end
EarliestArrivalWorkspace{T}(n::Integer) where {T} = EarliestArrivalWorkspace{T}(MinHeap{Tuple{T,Int}}(), zeros(Bool, n), zeros(Int, n))

distance_workspace(g::IncidentLists{V,T}, dt::DistanceType) where {V,T} =
    IncidentListsWorkspace{T,_cost_type(dt, T)}(num_nodes(g))
distance_workspace(g::IncidentLists{V,T}, ::EarliestArrival) where {V,T} =
    EarliestArrivalWorkspace{T}(num_nodes(g))

function _reset!(ws::IncidentListsWorkspace)
    _reset!(ws.fronts)
    for x in (ws.node, ws.arr, ws.cost, ws.parent, ws.edge, ws.deleted)
        empty!(x)
    end
    fill!(ws.visited, false)
    empty!(ws.heap)
    return ws
end

@inline function _push_label!(ws::IncidentListsWorkspace, v, a, c, parent, edge)
    push!(ws.node, v)
    push!(ws.arr, a)
    push!(ws.cost, c)
    push!(ws.parent, parent)
    push!(ws.edge, edge)
    push!(ws.deleted, false)
    return length(ws.node)
end

# priority of a label in the heap: the duration for fastest and latest departure
# paths, the cost otherwise
@inline _priority(::Union{Fastest,LatestDeparture}, a, c) = a + c
@inline _priority(::DistanceType, a, c) = c
# the search can stop as soon as all nodes are settled (not for latest departure,
# whose value is not monotone in the settling order)
_can_stop_early(::LatestDeparture) = false
_can_stop_early(::DistanceType) = true

# Runs the label-setting algorithm from s and writes the distances into `dist`.
# Stops when `target` is settled and returns the label that settled it (0 if none).
function _il_labels!(dist::AbstractVector{C}, ws::IncidentListsWorkspace{T,C}, g::IncidentLists{V,T}, s::Int,
                     ti::TimeInterval{T}, dt::DistanceType, target::Int=0) where {V,T,C}
    a0, b0 = ti
    n = num_nodes(g)
    fill!(dist, _unreachable(dt, C))
    dist[s] = dt isa LatestDeparture ? b0 : zero(C)
    _reset!(ws)
    _push_label!(ws, s, typemin(T), zero(C), 0, 0)
    push!(ws.heap, (zero(C), 1))
    deleted = ws.deleted
    mark_deleted = y -> (@inbounds deleted[y[3]] = true; nothing)
    nvisited = 0
    @inbounds while !isempty(ws.heap) && !(_can_stop_early(dt) && nvisited == n)
        _, l = pop!(ws.heap)
        deleted[l] && continue
        u = ws.node[l]
        if !ws.visited[u]
            ws.visited[u] = true
            nvisited += 1
        end
        u == target && return l
        src = l == 1
        la, lc = ws.arr[l], ws.cost[l]
        for k in _first_edge_from(g, u, src ? a0 : max(la, a0)):g.offsets[u+1]-1
            e = g.edges[k]
            e.t > b0 && break
            arr = e.t + e.tt
            arr > b0 && continue
            c = _next_cost(dt, src ? _init_cost(dt, e.t, C) : lc, e)
            v = Int(e.v)
            id = length(ws.node) + 1
            if _front_insert!(_list!(ws.fronts, v), 1, (arr, c, id), mark_deleted)
                _push_label!(ws, v, arr, c, l, k)
                val = _value(dt, arr, c)
                _better(dt, val, dist[v]) && (dist[v] = val)
                push!(ws.heap, (_priority(dt, arr, c), id))
            end
        end
    end
    return 0
end

# Earliest arrival times with Dijkstra's algorithm; `pred[v]` is the edge with which
# v is reached first.
function _il_ea!(arr::AbstractVector{T}, ws::EarliestArrivalWorkspace{T}, g::IncidentLists{V,T}, s::Int,
                 ti::TimeInterval{T}, target::Int=0) where {V,T}
    a0, b0 = ti
    fill!(arr, typemax(T))
    arr[s] = typemin(T)
    fill!(ws.done, false)
    heap = empty!(ws.heap)
    push!(heap, (typemin(T), s))
    @inbounds while !isempty(heap)
        au, u = pop!(heap)
        ws.done[u] && continue
        ws.done[u] = true
        u == target && break
        for k in _first_edge_from(g, u, max(au, a0)):g.offsets[u+1]-1
            e = g.edges[k]
            e.t > b0 && break
            av = e.t + e.tt
            v = Int(e.v)
            if av <= b0 && av < arr[v]
                arr[v] = av
                ws.pred[v] = k
                push!(heap, (av, v))
            end
        end
    end
    arr[s] = zero(T)
    return arr
end

function temporal_distances!(dist::AbstractVector, ws, g::IncidentLists{V,T}, s::Integer, dt::DistanceType,
                             ti=time_interval(g)) where {V,T}
    length(dist) == num_nodes(g) || throw(DimensionMismatch("dist must have length $(num_nodes(g))"))
    ti = _interval(T, ti)
    s = _check_node(g, s)
    if dt isa EarliestArrival
        _il_ea!(dist, ws, g, s, ti)
    else
        _il_labels!(dist, ws, g, s, ti, dt)
    end
    return dist
end
