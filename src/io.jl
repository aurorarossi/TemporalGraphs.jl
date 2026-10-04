@inline _isdelim(c::UInt8) = c == UInt8(' ') || c == UInt8('\t') || c == UInt8(',') || c == UInt8('\r')

# integer token str[i:j] parsed without allocations
@inline function _parse_token(::Type{T}, bytes, str, i, j, line) where {T<:Integer}
    neg = @inbounds bytes[i] == UInt8('-')
    k = neg ? i + 1 : i
    k > j && throw(ArgumentError("invalid integer on line $line"))
    x = zero(T)
    @inbounds for p in k:j
        d = bytes[p] - UInt8('0')
        d > 9 && throw(ArgumentError("invalid integer on line $line"))
        x = T(10) * x + T(d)
    end
    return neg ? -x : x
end

@inline function _parse_token(::Type{T}, bytes, str, i, j, line) where {T<:AbstractFloat}
    x = tryparse(T, SubString(str, i, j))
    x === nothing && throw(ArgumentError("invalid number on line $line"))
    return x
end

# Parse the edge file into (u, v, t, tt) tuples with the ids of the file.
function _parse_edges(str::String, ::Type{T}, default_tt::T) where {T<:Real}
    bytes = codeunits(str)
    out = Tuple{Int,Int,T,T}[]
    n = length(bytes)
    i = 1
    line = 1
    @inbounds while i <= n
        k = 0
        u = v = 0
        t = tt = default_tt
        while i <= n && _isdelim(bytes[i])
            i += 1
        end
        if i <= n && (bytes[i] == UInt8('#') || bytes[i] == UInt8('%'))
            while i <= n && bytes[i] != UInt8('\n')
                i += 1
            end
        end
        while i <= n && bytes[i] != UInt8('\n')
            if _isdelim(bytes[i])
                i += 1
                continue
            end
            j = i
            while j < n && !_isdelim(bytes[j+1]) && bytes[j+1] != UInt8('\n')
                j += 1
            end
            k += 1
            if k <= 2
                x = _parse_token(Int, bytes, str, i, j, line)
                k == 1 ? (u = x) : (v = x)
            elseif k == 3
                t = _parse_token(T, bytes, str, i, j, line)
            elseif k == 4
                tt = _parse_token(T, bytes, str, i, j, line)
            end
            i = j + 1
        end
        i += 1
        k >= 3 && push!(out, (u, v, t, tt))
        line += 1
    end
    return out
end

"""
    load_ordered_edge_list(filename; directed = true, default_tt = 1,
                           time_type = Int64, node_type = Int32)

Load a temporal graph from a text file in which each line `u v t [tt]` is a temporal
edge (separated by spaces, tabs or commas; lines starting with `#` or `%` and lines
with fewer than three values are skipped; missing transition times are `default_tt`;
further columns are ignored). Use `time_type = Float64` for real-valued times.

Node ids are mapped to `1:n` in order of first appearance, the ids of the file are
available via [`original_id`](@ref). If `directed` is `false`, every edge is inserted
in both directions.
"""
function load_ordered_edge_list(filename::AbstractString; directed::Bool=true, default_tt::Real=1,
                                time_type::Type{T}=Int64, node_type::Type{V}=Int32) where {T<:Real,V<:Integer}
    raw = _parse_edges(read(filename, String), T, T(default_tt))
    return _from_contacts(V, T, raw, directed)
end

# Build an OrderedEdgeList from (u, v, t, tt) tuples with arbitrary node ids, mapped
# to 1:n in order of first appearance.
function _from_contacts(::Type{V}, ::Type{T}, raw, directed::Bool) where {V,T}
    ids = Dict{Int,Int}()
    original = Int[]
    es = Vector{TemporalEdge{V,T}}(undef, 0)
    sizehint!(es, directed ? length(raw) : 2 * length(raw))
    id!(x) = get!(ids, x) do
        push!(original, x)
        length(original)
    end
    for (a, b, t, tt) in raw
        u = id!(a)
        v = id!(b)
        push!(es, TemporalEdge{V,T}(u, v, t, tt))
        directed || push!(es, TemporalEdge{V,T}(v, u, t, tt))
    end
    return OrderedEdgeList(length(original), sort!(es); original_ids=original)
end

"""
    OrderedEdgeList(table; u = 1, v = 2, t = 3, tt = nothing, directed = true,
                    default_tt = 1, node_type = Int32)

Build a temporal graph from any table supported by Tables.jl (e.g., a `DataFrame`
or a `CSV.File`). `u`, `v`, `t` and `tt` are the column names or indices of the
tails, heads, time stamps and transition times; without a `tt` column all
transition times are `default_tt`. Node ids (integers) are mapped to `1:n` in order
of first appearance, as in [`load_ordered_edge_list`](@ref).
"""
function OrderedEdgeList(table; u=1, v=2, t=3, tt=nothing, directed::Bool=true, default_tt::Real=1,
                         node_type::Type{V}=Int32) where {V<:Integer}
    Tables.istable(table) || throw(ArgumentError("expected a Tables.jl compatible table"))
    cols = Tables.columns(table)
    col(c) = c isa Integer ? Tables.getcolumn(cols, Int(c)) : Tables.getcolumn(cols, Symbol(c))
    us, vs, ts = col(u), col(v), col(t)
    tts = tt === nothing ? nothing : col(tt)
    T = tts === nothing ? promote_type(eltype(ts), typeof(default_tt)) : promote_type(eltype(ts), eltype(tts))
    raw = [(Int(us[i]), Int(vs[i]), T(ts[i]), T(tts === nothing ? default_tt : tts[i])) for i in eachindex(us)]
    return _from_contacts(V, T, raw, directed)
end

"""
    load_incident_lists(filename; kwargs...)

Load a temporal graph in [`IncidentLists`](@ref) representation, see
[`load_ordered_edge_list`](@ref).
"""
load_incident_lists(filename::AbstractString; kwargs...) =
    to_incident_lists(load_ordered_edge_list(filename; kwargs...))

"""
    load_trs_graph(filename; kwargs...)

Load a temporal graph in [`TRSGraph`](@ref) representation, see
[`load_ordered_edge_list`](@ref).
"""
load_trs_graph(filename::AbstractString; kwargs...) =
    to_trs_graph(load_ordered_edge_list(filename; kwargs...))

"""
    save_ordered_edge_list(filename, g::OrderedEdgeList; original_ids = true)

Write `g` in the format read by [`load_ordered_edge_list`](@ref), one line `u v t tt`
per edge. With `original_ids = true` the node ids of the input file are used.
"""
function save_ordered_edge_list(filename::AbstractString, g::OrderedEdgeList; original_ids::Bool=true)
    open(filename, "w") do io
        for e in edges(g)
            u = original_ids ? g.original_ids[e.u] : e.u
            v = original_ids ? g.original_ids[e.v] : e.v
            println(io, u, ' ', v, ' ', e.t, ' ', e.tt)
        end
    end
    return filename
end
