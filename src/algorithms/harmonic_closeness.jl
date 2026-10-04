# Temporal closeness integrated over the departure time, following Crescenzi,
# Magnien and Marino, "Finding top-k nodes for temporal closeness in large temporal
# graphs" (Algorithms, 2020):
#
#     C(s) = 1 / ((n - 1)(B - A)) Σ_{v ≠ s} ∫_A^B dτ / ℓ_τ(s, v),
#
# where ℓ_τ(s, v) is the latency from s to v when leaving s at time τ (earliest
# arrival at v over the paths whose first edge departs at or after τ, minus τ), and
# unreachable nodes contribute 0.
#
# The latency is a step function of τ: if (a_1, r_1), (a_2, r_2), ... are the Pareto
# optimal pairs (earliest arrival, latest departure) of the paths from s to v, sorted
# by increasing arrival (and departure), then ℓ_τ = a_k - τ for τ in (r_{k-1}, r_k], so
#
#     ∫_A^B dτ / ℓ_τ = Σ_k log((a_k - max(A, r_{k-1})) / (a_k - r_k)).
#
# The Pareto fronts are the labels of the fastest path scan on the edge stream, kept
# without pruning. For the sampling approximation the fronts towards a target are
# computed on the reversed temporal graph, where a path from u to d with departure r
# and arrival a becomes a path from d to u with departure A + B - a and arrival
# A + B - r.

# Scan from s keeping at every node the complete front of labels (arrival, -departure).
function _latency_fronts!(ws::StreamWorkspace{T,T}, g::OrderedEdgeList{V,T}, s::Int, ti::TimeInterval{T}) where {V,T}
    _reset!(ws)
    fronts = ws.fronts
    _scan_stream(g, ti) do e
        @inline
        @inbounds begin
            t = e.t
            u = Int(e.u)
            if u == s
                fu = _list!(fronts, u)
                changed = _front_insert!(fu, 1, (t, zero(T) - t))
            else
                _has_list(fronts, u) || return false
                fu = fronts.lists[u]
                changed = false
            end
            i = _last_leq(fu, 1, t)
            i < 1 && return changed
            return _front_insert!(_list!(fronts, Int(e.v)), 1, (t + e.tt, fu[i][2])) | changed
        end
    end
    return ws
end

# ∫_A^B dτ / ℓ_τ for the front f of labels (arrival, -departure); with `reversed`
# the labels are those of the reversed graph.
function _latency_integral(f::Vector{Tuple{T,T}}, A::T, B::T, reversed::Bool) where {T}
    total = 0.0
    lo = A
    if reversed  # (A + B - r, -(A + B - a)): increasing r is decreasing position
        @inbounds for k in length(f):-1:1
            r = A + B - f[k][1]
            a = A + B + f[k][2]
            r > lo && (total += log((a - lo) / (a - r)))
            lo = max(lo, r)
        end
    else
        @inbounds for k in eachindex(f)
            a, r = f[k][1], -f[k][2]
            r > lo && (total += log((a - lo) / (a - r)))
            lo = max(lo, r)
        end
    end
    return total
end

function _check_positive_tt(g::OrderedEdgeList, ti)
    T = time_type(g)
    a, b = _interval(T, ti)
    any(e -> a <= e.t && e.t + e.tt <= b && !(e.tt > 0), edges(g)) && throw(ArgumentError(
        "the integrated temporal closeness requires positive transition times"))
    b > a || throw(ArgumentError("the time interval must have positive length"))
    return (a, b)
end

"""
    temporal_harmonic_closeness(g::OrderedEdgeList, ti = time_interval(g))

Temporal closeness of all nodes averaged over the departure time (Crescenzi, Magnien
and Marino, 2020):

    C(s) = 1 / ((n - 1)(b - a)) Σ_{v ≠ s} ∫_a^b dτ / ℓ_τ(s, v),

where `(a, b) = ti` and the *latency* `ℓ_τ(s, v)` is the earliest arrival time at `v`
when leaving `s` at time `τ` or later, minus `τ`. Unreachable nodes contribute 0.
Only edges inside `ti` are used and transition times must be positive. The integral
is computed exactly from the Pareto fronts of (arrival, departure) pairs in one scan
of the edges per node; nodes are processed in parallel.

# References

- P. Crescenzi, C. Magnien, and A. Marino. *Finding top-k nodes for temporal closeness in large temporal graphs.* Algorithms 13(9), 2020. [DOI](https://doi.org/10.3390/a13090211)
"""
function temporal_harmonic_closeness(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    A, B = _check_positive_tt(g, ti)
    n = num_nodes(g)
    n < 2 && return zeros(Float64, n)
    result = Vector{Float64}(undef, n)
    norm = (n - 1) * Float64(B - A)
    @tasks for s in 1:n
        @local ws = StreamWorkspace{T,T}(n)
        _latency_fronts!(ws, g, s, (A, B))
        c = 0.0
        for v in ws.fronts.touched
            v != s && (c += _latency_integral(ws.fronts.lists[v], A, B, false))
        end
        result[s] = c / norm
    end
    return result
end

"""
    temporal_harmonic_closeness(g::OrderedEdgeList, s, ti = time_interval(g))

Temporal closeness of node `s` averaged over the departure time, see above.
"""
function temporal_harmonic_closeness(g::OrderedEdgeList{V,T}, s::Integer, ti=time_interval(g)) where {V,T}
    A, B = _check_positive_tt(g, ti)
    s = _check_node(g, s)
    n = num_nodes(g)
    n < 2 && return 0.0
    ws = StreamWorkspace{T,T}(n)
    _latency_fronts!(ws, g, s, (A, B))
    c = sum((_latency_integral(ws.fronts.lists[v], A, B, false) for v in ws.fronts.touched if v != s); init=0.0)
    return c / ((n - 1) * Float64(B - A))
end

"""
    temporal_harmonic_closeness_approximation(g::OrderedEdgeList, samples, ti = time_interval(g);
                                              rng = Random.default_rng())

Unbiased estimate of [`temporal_harmonic_closeness`](@ref) of all nodes (Crescenzi,
Magnien and Marino, 2020): for `samples` uniformly random target nodes `d`, the
contributions of `d` to the closeness of all nodes are computed with one backward
scan (a forward scan of the reversed temporal graph), and the sums are scaled by
`n / samples`. The estimates are most accurate for the nodes with the highest
closeness, which makes them suited to find the top-k nodes.

# References

- P. Crescenzi, C. Magnien, and A. Marino. *Finding top-k nodes for temporal closeness in large temporal graphs.* Algorithms 13(9), 2020. [DOI](https://doi.org/10.3390/a13090211)
"""
function temporal_harmonic_closeness_approximation(g::OrderedEdgeList{V,T}, samples::Integer, ti=time_interval(g);
                                                   rng::AbstractRNG=Random.default_rng()) where {V,T}
    A, B = _check_positive_tt(g, ti)
    samples >= 1 || throw(ArgumentError("at least one sample is needed"))
    n = num_nodes(g)
    n < 2 && return zeros(Float64, n)
    rg = reverse(g, (A, B))
    targets = [rand(rng, 1:n) for _ in 1:samples]
    total = tmapreduce(.+, chunks(targets; n=min(samples, 4 * Threads.nthreads()))) do ds
        ws = StreamWorkspace{T,T}(n)
        acc = zeros(Float64, n)
        for d in ds
            _latency_fronts!(ws, rg, d, (A, B))
            for u in ws.fronts.touched
                u != d && (acc[u] += _latency_integral(ws.fronts.lists[u], A, B, true))
            end
        end
        acc
    end
    return total .* (n / (samples * (n - 1) * Float64(B - A)))
end
