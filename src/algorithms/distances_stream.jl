# Temporal distances on the edge stream (OrderedEdgeList) in one pass over the edges.
#
# Fastest, latest departure, minimum hops and minimum transition times paths keep at
# every node a Pareto front of labels (a, c): a path arriving at time `a` with cost
# `c`. The front is sorted by increasing `a` and strictly decreasing `c` (for fastest
# and latest departure paths c = -(start time)). As edges are scanned
# chronologically, labels older than the one used by the last query at a node can
# never be the best one again, so each front keeps a moving start index (`head`).

# Type of the label costs and of the reported distances.
_cost_type(::MinimumHops, ::Type{T}) where {T} = Int
_cost_type(::DistanceType, ::Type{T}) where {T} = T

"""
    distance_eltype(g, dt::DistanceType)

Element type of the distance vectors computed for `g` and `dt`: `Int` for
[`MinimumHops`](@ref), the time type of `g` otherwise.
"""
distance_eltype(g, dt::DistanceType) = _cost_type(dt, time_type(g))

"""
    StreamWorkspace

Reusable buffers for temporal distance computations on an [`OrderedEdgeList`](@ref),
created by [`distance_workspace`](@ref).
"""
struct StreamWorkspace{T,C}
    fronts::LazyLists{Tuple{T,C}}
    head::Vector{Int}
end
StreamWorkspace{T,C}(n::Integer) where {T,C} = StreamWorkspace{T,C}(LazyLists{Tuple{T,C}}(n), ones(Int, n))

"""
    distance_workspace(g, dt::DistanceType)

Buffers for [`temporal_distances!`](@ref). Repeated single-source computations with
the same workspace do not allocate once the buffers have grown to their final size.
A workspace must not be shared between tasks running in parallel.
"""
distance_workspace(g::OrderedEdgeList{V,T}, dt::DistanceType) where {V,T} =
    StreamWorkspace{T,_cost_type(dt, T)}(num_nodes(g))
distance_workspace(::OrderedEdgeList, ::EarliestArrival) = nothing

function _reset!(ws::StreamWorkspace)
    @inbounds for v in ws.fronts.touched
        ws.head[v] = 1
    end
    _reset!(ws.fronts)
    return ws
end

# largest index i ≥ h with f[i][1] ≤ a, or h - 1
@inline function _last_leq(f::Vector{<:Tuple}, h::Int, a)
    lo, hi = h, length(f)
    @inbounds while lo <= hi
        mid = (lo + hi) >>> 1
        if f[mid][1] <= a
            lo = mid + 1
        else
            hi = mid - 1
        end
    end
    return hi
end

# Insert the label `x = (a, c, ...)` into the front `f` (active part f[h:end]) unless
# it is dominated. `onremove(y)` is called for every label `y` that the new label
# dominates. Returns true iff the label was inserted.
@inline function _front_insert!(f::Vector{L}, h::Int, x::L, onremove::F=nothing) where {L<:Tuple,F}
    a, c = x[1], x[2]
    j = _last_leq(f, h, a)
    @inbounds if j >= h && f[j][2] <= c
        return false  # dominated
    end
    first = (j >= h && @inbounds(f[j][1]) == a) ? j : j + 1
    k = j + 1
    @inbounds while k <= length(f) && f[k][2] >= c
        k += 1
    end
    # labels first:k-1 are dominated by x
    if onremove !== nothing
        @inbounds for i in first:k-1
            onremove(f[i])
        end
    end
    _replace!(f, first, k - 1, x)
    return true
end

# Replace f[first:last] (empty if last == first - 1) by the single element x. Unlike
# deleteat!/insert!, which may move the beginning of the vector inside its memory
# and later force a reallocation, the vector only grows and shrinks at its end, so
# a front that is reused does not allocate.
@inline function _replace!(f::Vector, first::Int, last::Int, x)
    n = length(f)
    if last < first
        push!(f, x)
        first <= n && copyto!(f, first + 1, f, first, n - first + 1)
    elseif last > first
        last < n && copyto!(f, first + 1, f, last + 1, n - last)
        resize!(f, n - (last - first))
    end
    @inbounds f[first] = x
    return f
end

# remove f[1:k] keeping the memory layout (see _replace!)
@inline function _drop_first!(f::Vector, k::Int)
    n = length(f)
    k < n && copyto!(f, 1, f, k + 1, n - k)
    resize!(f, n - k)
    return f
end

# Edges with the same time stamp and transition time 0 can be traversed one after
# the other, so a single pass in stream order is not enough. They form a contiguous
# block at the beginning of each time stamp (edges are sorted by (t, tt, u, v)), which
# is scanned until nothing changes. `relax!(e)` processes one edge and returns true
# iff it changed the state. Graphs without such edges use a plain loop.
@inline function _scan_stream(relax!::F, g::OrderedEdgeList{V,T}, ti::TimeInterval{T}) where {F,V,T}
    a0, b0 = ti
    es = g.edges
    if !g.has_zero_tt && g.max_arrival <= b0  # no edge can end after the interval
        @inbounds for e in es
            t = e.t
            t > b0 && break
            t < a0 && continue
            relax!(e)
        end
        return nothing
    elseif !g.has_zero_tt
        @inbounds for e in es
            t = e.t
            t > b0 && break
            (t < a0 || t + e.tt > b0) && continue
            relax!(e)
        end
        return nothing
    end
    i = 1
    m = length(es)
    @inbounds while i <= m
        e = es[i]
        t = e.t
        t > b0 && break
        if t < a0 || t + e.tt > b0
            i += 1
            continue
        end
        if iszero(e.tt) && i < m && es[i+1].t == t && iszero(es[i+1].tt)
            j = i + 1
            while j < m && es[j+1].t == t && iszero(es[j+1].tt)
                j += 1
            end
            changed = true
            while changed
                changed = false
                for k in i:j
                    changed |= relax!(es[k])
                end
            end
            i = j + 1
        else
            relax!(e)
            i += 1
        end
    end
    return nothing
end

@inline _init_cost(::Union{Fastest,LatestDeparture}, t, ::Type{C}) where {C} = -t
@inline _init_cost(::Union{MinimumHops,MinimumTransitionTimes}, t, ::Type{C}) where {C} = zero(C)
@inline _next_cost(::Union{Fastest,LatestDeparture}, c, e) = c
@inline _next_cost(::MinimumHops, c, e) = c + 1
@inline _next_cost(::MinimumTransitionTimes, c, e) = c + e.tt
# value reported for a path arriving at time a with cost c
@inline _value(::Fastest, a, c) = a + c
@inline _value(::LatestDeparture, a, c) = -c
@inline _value(::Union{MinimumHops,MinimumTransitionTimes}, a, c) = c
@inline _better(::LatestDeparture, x, y) = x > y
@inline _better(::DistanceType, x, y) = x < y
@inline _unreachable(::LatestDeparture, ::Type{C}) where {C} = typemin(C)
@inline _unreachable(::DistanceType, ::Type{C}) where {C} = typemax(C)

function _stream_labels!(dist::AbstractVector{C}, ws::StreamWorkspace{T,C}, g::OrderedEdgeList{V,T},
                         s::Int, ti::TimeInterval{T}, dt::DistanceType) where {V,T,C}
    fill!(dist, _unreachable(dt, C))
    dist[s] = dt isa LatestDeparture ? ti[2] : zero(C)
    _reset!(ws)
    fronts, head = ws.fronts, ws.head
    _scan_stream(g, ti) do e
        @inline  # the scan is the hot loop: avoid a call per edge
        @inbounds begin
            t = e.t
            u = Int(e.u)
            if u == s
                fu = _list!(fronts, u)
                changed = _front_insert!(fu, head[u], (t, _init_cost(dt, t, C)))
            else
                _has_list(fronts, u) || return false
                fu = fronts.lists[u]
                changed = false
            end
            h = head[u]
            i = _last_leq(fu, h, t)
            i < h && return changed
            if i > 32 && 2i > length(fu)  # drop labels that can no longer be used
                _drop_first!(fu, i - 1)
                i = 1
            end
            head[u] = i
            c = _next_cost(dt, fu[i][2], e)
            arr = t + e.tt
            v = Int(e.v)
            val = _value(dt, arr, c)
            _better(dt, val, dist[v]) && (dist[v] = val)
            return _front_insert!(_list!(fronts, v), head[v], (arr, c)) | changed
        end
    end
    return dist
end

function _stream_ea!(arr::AbstractVector{T}, g::OrderedEdgeList{V,T}, s::Int, ti::TimeInterval{T}) where {V,T}
    fill!(arr, typemax(T))
    arr[s] = typemin(T)
    reached = Ref(1)
    _scan_stream(g, ti) do e
        @inline
        @inbounds if arr[e.u] <= e.t && arr[e.v] > e.t + e.tt
            arr[e.v] == typemax(T) && (reached[] += 1)
            arr[e.v] = e.t + e.tt
            return true
        end
        return false
    end
    arr[s] = zero(T)
    return reached[]
end

"""
    temporal_distances!(dist, ws, g, s, dt::DistanceType, ti = time_interval(g))

In-place version of [`temporal_distances`](@ref): writes the distances from `s` into
`dist` (a vector of length `num_nodes(g)` and element type
[`distance_eltype`](@ref)`(g, dt)`), using the buffers `ws =`
[`distance_workspace`](@ref)`(g, dt)`. Returns `dist`.
"""
function temporal_distances!(dist::AbstractVector, ws, g::OrderedEdgeList{V,T}, s::Integer, dt::DistanceType,
                             ti=time_interval(g)) where {V,T}
    length(dist) == num_nodes(g) || throw(DimensionMismatch("dist must have length $(num_nodes(g))"))
    ti = _interval(T, ti)
    s = _check_node(g, s)
    if dt isa EarliestArrival
        _stream_ea!(dist, g, s, ti)
    else
        _stream_labels!(dist, ws, g, s, ti, dt)
    end
    return dist
end

"""
    number_of_reachable_nodes(g::OrderedEdgeList, s, ti = time_interval(g))

Number of nodes reachable from `s` inside `ti`, including `s` itself.
"""
function number_of_reachable_nodes(g::OrderedEdgeList{V,T}, s::Integer, ti=time_interval(g)) where {V,T}
    arr = Vector{T}(undef, num_nodes(g))
    return _stream_ea!(arr, g, _check_node(g, s), _interval(T, ti))
end
