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

# time-sorted copy of `es` (sorted by (t, tt, u, v) if it was not sorted by time)
function _chronological(::Type{E}, es) where {E}
    v = collect(E, es)
    issorted(v; by=e -> e.t) || sort!(v)
    return v
end
