# Top-k temporal closeness, Oettershagen and Mutzel, ICDM 2020.
#
# The closeness of a candidate is computed with the label-setting algorithm, which
# is stopped as soon as an upper bound of the closeness falls below the k-th best
# value found so far. Labels are settled by increasing distance, so with `dmin` the
# smallest priority in the heap, discovered but unsettled nodes have distance ≥ dmin
# and undiscovered nodes have distance ≥ dmin + min_tt (the minimum transition time).
# (TGLib uses the minimum time difference between consecutive edges instead, which
# is not always a valid bound.)

@inline _upper_bound(exact, num_t, rest, dmin, min_tt) = exact + num_t / dmin + rest / (dmin + min_tt)

function _topk_candidate!(topk::TopkResult, ws::IncidentListsWorkspace{T,C}, dist::Vector{C},
                          discovered::Vector{Bool}, g::IncidentLists{V,T}, s::Int, ti::TimeInterval{T},
                          dt::DistanceType, min_tt::T) where {V,T,C}
    a0, b0 = ti
    n = num_nodes(g)
    fill!(dist, typemax(C))
    dist[s] = zero(C)
    fill!(discovered, false)
    _reset!(ws)
    _push_label!(ws, s, typemin(T), zero(C), 0, 0)
    push!(ws.heap, (zero(C), 1))
    deleted = ws.deleted
    mark_deleted = y -> (@inbounds deleted[y[3]] = true; nothing)
    num_f = 0
    num_t = 0
    exact = 0.0
    @inbounds while !isempty(ws.heap) && num_f < n
        _, l = pop!(ws.heap)
        deleted[l] && continue
        u = ws.node[l]
        if !ws.visited[u]
            ws.visited[u] = true
            num_f += 1
            discovered[u] && (num_t -= 1)
            dist[u] > 0 && (exact += 1.0 / dist[u])
        end
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
                val < dist[v] && (dist[v] = val)
                push!(ws.heap, (_priority(dt, arr, c), id))
                if !discovered[v] && !ws.visited[v]
                    discovered[v] = true
                    num_t += 1
                end
            end
        end
        if _is_full(topk) && !isempty(ws.heap)
            dmin = first(ws.heap)[1]
            bound = _upper_bound(exact, num_t, n - num_f - num_t, dmin, min_tt)
            bound < _min_topk(topk) && return topk
        end
    end
    insert!(topk, s, exact)
    return topk
end

"""
    compute_topk_closeness(g, k, dt::DistanceType, ti = time_interval(g))

The `k` nodes with the highest temporal closeness as a vector of `(node, closeness)`
pairs by decreasing closeness; nodes with equal closeness are all returned.

For [`Fastest`](@ref) and [`MinimumTransitionTimes`](@ref) the pruning algorithm of
Oettershagen and Mutzel (ICDM 2020) is used, for other distance types the closeness
of all nodes is computed.

# References

- L. Oettershagen and P. Mutzel. *Efficient top-k temporal closeness calculation in temporal networks.* ICDM, 2020. [DOI](https://doi.org/10.1109/ICDM50108.2020.00049)
"""
function compute_topk_closeness(g::IncidentLists{V,T}, k::Integer, dt::DistanceType, ti=time_interval(g)) where {V,T}
    ti = _interval(T, ti)
    topk = TopkResult(k)
    if dt isa Fastest || dt isa MinimumTransitionTimes
        n = num_nodes(g)
        min_tt = isempty(g.edges) ? zero(T) : minimum(e.tt for e in g.edges)
        ws = distance_workspace(g, dt)
        dist = Vector{distance_eltype(g, dt)}(undef, n)
        discovered = zeros(Bool, n)
        for s in sortperm([out_degree(g, u) for u in 1:n]; rev=true)
            _topk_candidate!(topk, ws, dist, discovered, g, s, ti, dt, min_tt)
        end
    else
        c = temporal_closeness(g, dt, ti)
        for s in sortperm(c; rev=true)
            insert!(topk, s, c[s])
        end
    end
    return results(topk)
end

function compute_topk_closeness(g::OrderedEdgeList, k::Integer, dt::DistanceType, ti=time_interval(g))
    (dt isa Fastest || dt isa MinimumTransitionTimes) && return compute_topk_closeness(to_incident_lists(g), k, dt, ti)
    topk = TopkResult(k)
    c = temporal_closeness(g, dt, ti)
    for s in sortperm(c; rev=true)
        insert!(topk, s, c[s])
    end
    return results(topk)
end
