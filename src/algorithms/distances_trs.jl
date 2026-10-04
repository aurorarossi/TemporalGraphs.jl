# Temporal distances on the time-respecting static graph. An edge leaving the
# time-node x departs at time(x); edges between time-nodes of the same temporal
# graph node (waiting edges and self-loops) cost nothing for hops and transition times.

"""
    TRSWorkspace

Reusable buffers for temporal distance computations on a [`TRSGraph`](@ref), created
by [`distance_workspace`](@ref).
"""
struct TRSWorkspace{C}
    visited::BitVector
    stack::Vector{Int}
    touched::Vector{Int}
    dist::Vector{C}
    heap::MinHeap{Tuple{C,Int}}
end

TRSWorkspace{C}(N::Integer) where {C} =
    TRSWorkspace{C}(falses(N), Int[], Int[], fill(typemax(C), N), MinHeap{Tuple{C,Int}}())

distance_workspace(g::TRSGraph{V,T}, dt::DistanceType) where {V,T} = TRSWorkspace{_cost_type(dt, T)}(num_trs_nodes(g))
# depth-first searches do not need the distance buffer
distance_workspace(g::TRSGraph{V,T}, ::Union{EarliestArrival,Fastest}) where {V,T} =
    TRSWorkspace{T}(falses(num_trs_nodes(g)), Int[], Int[], T[], MinHeap{Tuple{T,Int}}())
distance_workspace(::TRSGraph, ::LatestDeparture) =
    throw(ArgumentError("latest departure times are not implemented for TRSGraph"))

# Nodes reached by a search are recorded in `touched` (to reset only them) up to
# N/16 entries; beyond that the whole buffers are reset, which is cheaper.
@inline _touch!(ws::TRSWorkspace, x) = length(ws.touched) <= length(ws.visited) >> 4 && push!(ws.touched, x)

function _reset!(ws::TRSWorkspace{C}) where {C}
    reset_dist = !isempty(ws.dist)
    if length(ws.touched) > length(ws.visited) >> 4
        fill!(ws.visited, false)
        reset_dist && fill!(ws.dist, typemax(C))
    else
        @inbounds for x in ws.touched
            ws.visited[x] = false
            reset_dist && (ws.dist[x] = typemax(C))
        end
    end
    empty!(ws.touched)
    empty!(ws.stack)
    empty!(ws.heap)
    return ws
end

@inline function _visit!(ws::TRSWorkspace, x)
    @inbounds ws.visited[x] = true
    _touch!(ws, x)
    push!(ws.stack, x)
    return nothing
end

# first time-node of s with time ≥ a0 (0 if there is none inside the interval)
function _trs_start(g::TRSGraph, s, a0, b0)
    g.first_pos[s] == 0 && return 0
    p = g.first_pos[s]
    @inbounds while p <= g.last_pos[s] && g.time[p] < a0
        p += 1
    end
    return (p > g.last_pos[s] || @inbounds(g.time[p]) > b0) ? 0 : p
end

function _trs_ea!(arr::AbstractVector{T}, ws::TRSWorkspace, g::TRSGraph{V,T}, s::Int, ti::TimeInterval{T}) where {V,T}
    a0, b0 = ti
    fill!(arr, typemax(T))
    arr[s] = zero(T)
    _reset!(ws)
    start = _trs_start(g, s, a0, b0)
    start == 0 && return arr
    _visit!(ws, start)
    stack = ws.stack
    @inbounds while !isempty(stack)
        x = pop!(stack)
        tx = g.time[x]
        tx > b0 && continue
        for k in g.offsets[x]:g.offsets[x+1]-1
            arc = g.arcs[k]
            av = tx + arc.tt
            av > b0 && continue
            v = arc.node
            v != s && av < arr[v] && (arr[v] = av)
            ws.visited[arc.head] || _visit!(ws, Int(arc.head))
        end
    end
    return arr
end

function _trs_fastest!(durs::AbstractVector{T}, ws::TRSWorkspace, g::TRSGraph{V,T}, s::Int, ti::TimeInterval{T}) where {V,T}
    a0, b0 = ti
    fill!(durs, typemax(T))
    durs[s] = zero(T)
    _reset!(ws)
    g.first_pos[s] == 0 && return durs
    stack = ws.stack
    # later departures first: nodes reached from a later start are not revisited
    @inbounds for p in g.last_pos[s]:-1:g.first_pos[s]
        τ = g.time[p]
        τ < a0 && break
        (τ > b0 || ws.visited[p]) && continue
        _visit!(ws, p)
        while !isempty(stack)
            x = pop!(stack)
            tx = g.time[x]
            tx > b0 && continue
            for k in g.offsets[x]:g.offsets[x+1]-1
                arc = g.arcs[k]
                av = tx + arc.tt
                av > b0 && continue
                v = arc.node
                d = av - τ
                d < durs[v] && (durs[v] = d)
                ws.visited[arc.head] || _visit!(ws, Int(arc.head))
            end
        end
    end
    return durs
end

# minimum hops (dt = MinimumHops()) or transition times with Dijkstra's algorithm
function _trs_cost!(result::AbstractVector{C}, ws::TRSWorkspace{C}, g::TRSGraph{V,T}, s::Int,
                    ti::TimeInterval{T}, dt::DistanceType) where {V,T,C}
    a0, b0 = ti
    fill!(result, typemax(C))
    result[s] = zero(C)
    _reset!(ws)
    start = _trs_start(g, s, a0, b0)
    start == 0 && return result
    dist, heap = ws.dist, ws.heap
    dist[start] = zero(C)
    _touch!(ws, start)
    push!(heap, (zero(C), start))
    @inbounds while !isempty(heap)
        d, x = pop!(heap)
        d > dist[x] && continue
        tx = g.time[x]
        tx > b0 && continue
        u = g.tg_node[x]
        for k in g.offsets[x]:g.offsets[x+1]-1
            arc = g.arcs[k]
            tx + arc.tt > b0 && continue
            y = Int(arc.head)
            dy = arc.node == u ? d : d + (dt isa MinimumHops ? one(C) : C(arc.tt))
            if dy < dist[y]
                dist[y] == typemax(C) && _touch!(ws, y)
                dist[y] = dy
                arc.node != u && dy < result[arc.node] && (result[arc.node] = dy)
                push!(heap, (dy, y))
            end
        end
    end
    return result
end

function temporal_distances!(dist::AbstractVector, ws::TRSWorkspace, g::TRSGraph{V,T}, s::Integer, dt::DistanceType,
                             ti=time_interval(g)) where {V,T}
    length(dist) == num_nodes(g) || throw(DimensionMismatch("dist must have length $(num_nodes(g))"))
    ti = _interval(T, ti)
    s = _check_node(g, s)
    if dt isa EarliestArrival
        _trs_ea!(dist, ws, g, s, ti)
    elseif dt isa Fastest
        _trs_fastest!(dist, ws, g, s, ti)
    elseif dt isa LatestDeparture
        throw(ArgumentError("latest departure times are not implemented for TRSGraph"))
    else
        _trs_cost!(dist, ws, g, s, ti, dt)
    end
    return dist
end
