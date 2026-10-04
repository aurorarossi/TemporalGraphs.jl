# (σ - μ) / (σ + μ) of the inter-contact times of sorted time stamps `ts`, with the
# mean μ and the (uncorrected) standard deviation σ
function _burstiness(ts::AbstractVector{<:Real})
    ict = diff(ts)
    μ = mean(ict)
    σ = std(ict; corrected=false, mean=μ)
    return (σ - μ) / (σ + μ)
end

"""
    edge_burstiness(g::OrderedEdgeList, ti = time_interval(g))

Burstiness `(σ - μ) / (σ + μ)` of the inter-contact times of every static edge
`(u, v)` with at least two contacts with time stamps in `ti`, where `μ` and `σ` are
the mean and standard deviation. Returns a `Dict` from `(u, v)` to the burstiness.
"""
function edge_burstiness(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    contacts = [(Int(e.u), Int(e.v), e.t) for e in edges(g) if a <= e.t <= b]
    sort!(contacts)  # stable grouping by edge, time stamps stay sorted
    result = Dict{Tuple{Int,Int},Float64}()
    i = 1
    while i <= length(contacts)
        j = i
        while j < length(contacts) && contacts[j+1][1:2] == contacts[i][1:2]
            j += 1
        end
        if j > i
            result[contacts[i][1:2]] = _burstiness([c[3] for c in view(contacts, i:j)])
        end
        i = j + 1
    end
    return result
end

"""
    node_burstiness(g::OrderedEdgeList, ti = time_interval(g))

Burstiness of the inter-contact times of the outgoing edges of every node with time
stamps in `ti` (0 for nodes with fewer than two contacts).
"""
function node_burstiness(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    n = num_nodes(g)
    times = [T[] for _ in 1:n]
    for e in edges(g)
        a <= e.t <= b && push!(times[e.u], e.t)
    end
    return [length(ts) < 2 ? 0.0 : _burstiness(ts) for ts in times]
end
