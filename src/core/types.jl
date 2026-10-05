"""
    TemporalEdge{V<:Integer,T<:Real}(u, v, t, tt = 1)
    TemporalEdge(u, v, t, tt = 1)

A directed temporal edge from `u` to `v` that departs at time `t` and needs the
transition time `tt`, i.e., it arrives at `v` at time `t + tt`.

`V` is the node id type and `T` the time type (integer or floating point times).
Without explicit parameters node ids are stored as `Int32` (an edge with `Int64`
times takes 24 bytes) and `T` is the promoted type of `t` and `tt`.
Edges are ordered lexicographically by `(t, tt, u, v)`.
"""
struct TemporalEdge{V<:Integer,T<:Real}
    u::V
    v::V
    t::T
    tt::T
    TemporalEdge{V,T}(u, v, t, tt) where {V<:Integer,T<:Real} = new{V,T}(u, v, t, tt)
end

TemporalEdge{V}(u::Integer, v::Integer, t::Real, tt::Real=one(t)) where {V<:Integer} =
    TemporalEdge{V,promote_type(typeof(t), typeof(tt))}(u, v, t, tt)
TemporalEdge(u::Integer, v::Integer, t::Real, tt::Real=one(t)) = TemporalEdge{Int32}(u, v, t, tt)

Base.convert(::Type{TemporalEdge{V,T}}, e::TemporalEdge) where {V,T} = TemporalEdge{V,T}(e.u, e.v, e.t, e.tt)

@inline _key(e::TemporalEdge) = (e.t, e.tt, e.u, e.v)
Base.isless(a::TemporalEdge, b::TemporalEdge) = isless(_key(a), _key(b))
Base.show(io::IO, e::TemporalEdge) = print(io, "(", e.u, " ", e.v, " ", e.t, " ", e.tt, ")")

"""
    INF

`typemax(Int64)`, the distance of unreachable nodes in graphs with `Int64` times.
In general unreachable nodes get `typemax(T)` (`Inf` for floating point times).
"""
const INF = typemax(Int64)

"""
    TemporalDistances(A::AbstractArray)

The vector or matrix `A` of distances or times, wrapped without copying, that prints
the entries of unreachable nodes (`typemax`, i.e. [`INF`](@ref) for `Int64` times, and
`typemin` for latest departure times) as `∞` and `-∞`. Otherwise it behaves as `A`:
the entries keep their values, so `d[v] == INF` tests if `v` is unreachable, and
`parent(d)` returns `A`. Returned by [`temporal_distances`](@ref) and the functions
based on it, [`earliest_arrival_matrix`](@ref) and [`temporal_flooding_times`](@ref).
"""
struct TemporalDistances{T,N,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
    data::A
end

Base.parent(d::TemporalDistances) = d.data
Base.size(d::TemporalDistances) = size(d.data)
Base.IndexStyle(::Type{<:TemporalDistances{T,N,A}}) where {T,N,A} = IndexStyle(A)
Base.@propagate_inbounds Base.getindex(d::TemporalDistances, i::Int...) = d.data[i...]
Base.@propagate_inbounds Base.setindex!(d::TemporalDistances, x, i::Int...) = (d.data[i...] = x; d)
Base.showarg(io::IO, ::TemporalDistances{T}, toplevel) where {T} = print(io, "TemporalDistances{", T, "}")

# An entry as it is printed, and a lazy array of them for printing a TemporalDistances
struct _ShownDistance{T}
    x::T
end
_unreachable_sign(x::T) where {T} = x == typemax(T) ? 1 : (typemin(T) < zero(T) && x == typemin(T)) ? -1 : 0
function Base.show(io::IO, s::_ShownDistance)
    k = _unreachable_sign(s.x)
    k == 0 ? show(io, s.x) : print(io, k > 0 ? "∞" : "-∞")
end
Base.alignment(io::IO, s::_ShownDistance) =
    _unreachable_sign(s.x) == 0 ? Base.alignment(io, s.x) : (_unreachable_sign(s.x) > 0 ? (1, 0) : (2, 0))

struct _ShownDistances{T,N,A<:AbstractArray{T,N}} <: AbstractArray{_ShownDistance{T},N}
    data::A
end
Base.size(s::_ShownDistances) = size(s.data)
Base.IndexStyle(::Type{<:_ShownDistances{T,N,A}}) where {T,N,A} = IndexStyle(A)
Base.getindex(s::_ShownDistances, i::Int...) = _ShownDistance(s.data[i...])

function Base.show(io::IO, ::MIME"text/plain", d::TemporalDistances{T}) where {T}
    summary(io, d)
    isempty(d) && return
    println(io, ":")
    Base.print_array(IOContext(io, :typeinfo => T), _ShownDistances(d.data))
end
function Base.show(io::IO, d::TemporalDistances)
    s = _ShownDistances(d.data)
    show(IOContext(io, :typeinfo => typeof(s)), s)
end

"""
    TimeInterval{T}

A closed time interval `(a, b)`. Algorithms restricted to `(a, b)` only use edges
`e` with `a ≤ e.t` and `e.t + e.tt ≤ b`.
"""
const TimeInterval{T} = Tuple{T,T}

"""
    StaticWeightedEdge(u, v, weight)

A static (non-temporal) edge with an integer weight, e.g., the number of temporal
contacts between `u` and `v`.
"""
struct StaticWeightedEdge
    u::Int
    v::Int
    weight::Int
end

"""
    DistanceType

Abstract supertype of the notions of optimal temporal paths: [`EarliestArrival`](@ref),
[`Fastest`](@ref), [`LatestDeparture`](@ref), [`MinimumTransitionTimes`](@ref) and
[`MinimumHops`](@ref). They are singleton types, so the algorithms are specialized
at compile time for each of them.
"""
abstract type DistanceType end

"""Earliest arrival paths: the distance is the arrival time."""
struct EarliestArrival <: DistanceType end
"""Fastest paths: the distance is the duration (arrival time minus departure time)."""
struct Fastest <: DistanceType end
"""Latest departure paths: the latest time one can leave the source."""
struct LatestDeparture <: DistanceType end
"""Shortest paths: the distance is the sum of the transition times."""
struct MinimumTransitionTimes <: DistanceType end
"""Minimum hop paths: the distance is the number of edges."""
struct MinimumHops <: DistanceType end

"""
Shortest foremost walks: earliest arrival, ties broken by the number of edges.
Only used by [`temporal_betweenness`](@ref) and [`temporal_ego_betweenness`](@ref).
"""
struct ShortestForemost <: DistanceType end
"""
Shortest fastest walks: minimum duration, ties broken by the number of edges.
Only used by [`temporal_betweenness`](@ref) and [`temporal_ego_betweenness`](@ref).
"""
struct ShortestFastest <: DistanceType end
"""
Shortest latest walks: latest departure, ties broken by the number of edges.
Only used by [`temporal_betweenness`](@ref) and [`temporal_ego_betweenness`](@ref).
"""
struct ShortestLatest <: DistanceType end
"""
Prefix foremost paths (Buß et al., 2020): foremost paths whose every prefix is also
foremost, i.e., every node is reached at its earliest arrival time. Only used by
[`temporal_betweenness`](@ref) and [`temporal_ego_betweenness`](@ref).
"""
struct PrefixForemost <: DistanceType end

"""
    num_nodes(g)

Number of nodes of the temporal graph `g`.
"""
function num_nodes end

"""
    num_edges(g)

Number of (temporal) edges of `g`.
"""
function num_edges end

"""
    time_interval(g)

The time interval `(a, b)` spanned by `g`.
"""
function time_interval end

"""
    time_type(g)

The type of the time stamps of `g`.
"""
function time_type end

"""
    node_type(g)

The type of the node ids stored in the edges of `g`.
"""
function node_type end

# (min t, max t + tt) of a collection of edges
function _span(es::AbstractVector{TemporalEdge{V,T}}) where {V,T}
    isempty(es) && return (zero(T), zero(T))
    a = typemax(T)
    b = typemin(T)
    @inbounds for e in es
        a = min(a, e.t)
        b = max(b, e.t + e.tt)
    end
    return (a, b)
end

_interval(::Type{T}, ti) where {T} = (convert(T, ti[1]), convert(T, ti[2]))

function _check_node(g, u::Integer)
    1 <= u <= num_nodes(g) || throw(BoundsError(1:num_nodes(g), u))
    return Int(u)
end

function _check_edges(n, es)
    for e in es
        (1 <= e.u <= n && 1 <= e.v <= n) ||
            throw(ArgumentError("edge $e has a node id outside 1:$n"))
    end
    return nothing
end

# time-sorted copy of `es`: sorted by (t, tt, u, v) if it was not sorted by time.
# Otherwise the order is kept, except that the edges with transition time 0 are moved
# (stably) to the beginning of their time stamp, which the stream algorithms rely on.
_zero_first(e::TemporalEdge) = (e.t, !iszero(e.tt))
function _chronological(::Type{E}, es) where {E}
    v = collect(E, es)
    if !issorted(v; by=e -> e.t)
        sort!(v)
    elseif !issorted(v; by=_zero_first)
        sort!(v; by=_zero_first, alg=Base.Sort.DEFAULT_STABLE)
    end
    return v
end
