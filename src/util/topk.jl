"""
    TopkResult(k)

Collects the `k` largest values inserted with `insert!(r, id, value)`. Ids with equal
values are grouped, so more than `k` ids can be returned when there are ties.
"""
mutable struct TopkResult
    k::Int
    bound::Float64
    values::Vector{Float64}      # distinct values in decreasing order
    ids::Vector{Vector{Int}}     # ids for each value, in insertion order
end
TopkResult(k::Integer) = TopkResult(Int(k), 0.0, Float64[], Vector{Int}[])

# true iff there are already k distinct values
_is_full(r::TopkResult) = length(r.values) >= r.k
# the smallest of the stored values
_min_topk(r::TopkResult) = r.bound

Base.length(r::TopkResult) = length(r.values)

function Base.insert!(r::TopkResult, id::Integer, value::Real)
    value = Float64(value)
    _is_full(r) && r.bound > value && return r
    i = searchsortedfirst(r.values, value; rev=true)
    if i <= length(r.values) && r.values[i] == value
        push!(r.ids[i], id)
    else
        insert!(r.values, i, value)
        insert!(r.ids, i, [Int(id)])
    end
    if r.bound < value && length(r.values) > r.k && r.values[end] == r.bound
        pop!(r.values)
        pop!(r.ids)
    end
    r.bound = r.values[end]
    return r
end

"""
    results(r::TopkResult)

The stored `(id, value)` pairs by decreasing value.
"""
function results(r::TopkResult)
    out = Tuple{Int,Float64}[]
    for (v, ids) in zip(r.values, r.ids), id in ids
        push!(out, (id, v))
    end
    return out
end
