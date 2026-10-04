# Progressive sampling approximation of the temporal betweenness and sampled
# temporal distance statistics, following Cruciani, "MANTRA: temporal betweenness
# centrality approximation through sampling" (ECML PKDD 2024, arXiv:2304.08356).
#
# Every sample gives an unbiased estimate x(v) ∈ [0, 1] of the normalized temporal
# betweenness b(v) = B(v) / (n (n - 1)) of every node:
# - :pairs (the estimator of ONBRA): a uniform pair s ≠ z, x(v) = σ_sz(v) / σ_sz;
# - :sources: a uniform node s, x(v) = Σ_{z ≠ s} σ_sz(v) / σ_sz / (n - 1).
# The sample grows geometrically until a bound on the supremum deviation from
# Monte Carlo empirical Rademacher averages (Theorem 1 of the paper) is at most ε,
# or until a sufficient sample size computed from the variance of the estimators
# and the average path length (Theorem 3) is reached.

# State of the per-sample computations.
struct _Sampler{W,G,P}
    criterion::DistanceType
    wg::G                 # _WalkGraph (walk criteria)
    work::W               # _WalkWork
    pfm::P                # prefix foremost data
    n::Int
    x::Vector{Float64}    # estimate of the current sample
    touched::Vector{Int}
end

# prefix foremost paths: edges sorted by arrival and per-source buffers
struct _PrefixForemostData{V,T}
    es::Vector{TemporalEdge{V,T}}
    byarr::Vector{Int}
    ea::Vector{T}
    σ::Vector{Float64}
    τ::Vector{Float64}    # number of paths from each node to the target (pairs)
    δ::Vector{Float64}
end

function _sampler(g::OrderedEdgeList{V,T}, criterion) where {V,T}
    n = num_nodes(g)
    if criterion isa PrefixForemost
        es = _strict_edges(g, time_interval(g))
        pfm = _PrefixForemostData{V,T}(es, sortperm([e.t + e.tt for e in es]), Vector{T}(undef, n),
                                       zeros(n), zeros(n), zeros(n))
        return _Sampler(criterion, nothing, nothing, pfm, n, zeros(n), Int[])
    end
    wg = _walk_graph(g, time_interval(g), Inf)
    return _Sampler(criterion, wg, _WalkWork(wg, criterion, Float64), nothing, n, zeros(n), Int[])
end

@inline function _add_x!(sm::_Sampler, v, val)
    @inbounds if val != 0
        sm.x[v] == 0 && push!(sm.touched, v)
        sm.x[v] += val
    end
    return nothing
end

function _clear_x!(sm::_Sampler)
    for v in sm.touched
        sm.x[v] = 0.0
    end
    empty!(sm.touched)
    return sm
end

# earliest arrival times and numbers of prefix foremost paths from s
function _pfm_forward!(p::_PrefixForemostData{V,T}, s::Int) where {V,T}
    ea, σ, es = p.ea, p.σ, p.es
    fill!(ea, typemax(T))
    ea[s] = typemin(T)
    @inbounds for e in es
        (ea[e.u] <= e.t && e.t + e.tt < ea[e.v]) && (ea[e.v] = e.t + e.tt)
    end
    fill!(σ, 0.0)
    σ[s] = 1.0
    @inbounds for i in p.byarr
        e = es[i]
        (e.v != s && ea[e.u] <= e.t && e.t + e.tt == ea[e.v]) && (σ[e.v] += σ[e.u])
    end
    return p
end

@inline _pfm_edge(p, e, s) = e.v != s && p.ea[e.u] <= e.t && e.t + e.tt == p.ea[e.v]

# Estimate of the sample (s, z) (z = 0: all targets, divided by n - 1) in sm.x.
function _sample!(sm::_Sampler, s::Int, z::Int)
    _clear_x!(sm)
    n = sm.n
    if sm.criterion isa PrefixForemost
        p = sm.pfm
        _pfm_forward!(p, s)
        if z != 0  # paths from s to z: σ_sv · τ_v / σ_sz
            p.ea[z] == typemax(eltype(p.ea)) && return sm
            fill!(p.τ, 0.0)
            p.τ[z] = 1.0
            @inbounds for i in Iterators.reverse(p.byarr)
                e = p.es[i]
                _pfm_edge(p, e, s) && (p.τ[e.u] += p.τ[e.v])
            end
            @inbounds for v in 1:n
                (v != s && v != z && p.τ[v] > 0) && _add_x!(sm, v, p.σ[v] * p.τ[v] / p.σ[z])
            end
        else
            fill!(p.δ, 0.0)
            @inbounds for i in Iterators.reverse(p.byarr)
                e = p.es[i]
                _pfm_edge(p, e, s) && (p.δ[e.u] += p.σ[e.u] / p.σ[e.v] * (1 + p.δ[e.v]))
            end
            @inbounds for v in 1:n
                v != s && _add_x!(sm, v, p.δ[v] / (n - 1))
            end
        end
        return sm
    end
    w, wg, dt = sm.work, sm.wg, sm.criterion
    _walk_forward!(w, wg, s, dt)
    _walk_targets!(w, wg, s, dt, z)
    scale = z == 0 ? 1.0 / (n - 1) : 1.0
    _walk_backward!(w, wg, s, dt) do v, val
        v != z && _add_x!(sm, v, Float64(val) * scale)
    end
    if z == 0
        @inbounds for v in 1:n
            (v != s && w.σstar[v] > 0) && _add_x!(sm, v, -scale)
        end
    end
    return sm
end

_draw(rng, n, estimator) = estimator === :pairs ?
    (s = rand(rng, 1:n); z = rand(rng, 1:n-1); (s, z >= s ? z + 1 : z)) : (rand(rng, 1:n), 0)

# Upper bound on the maximum variance of the estimators from the empirical wimpy
# variance W (Proposition 1 of the paper).
function _variance_bound(W, r, δ)
    L = log(1 / δ) / r
    return min(0.25, W + L + sqrt(L^2 + 2W * L))
end

# Sufficient sample size of Theorem 3 (numerical approximation given in the paper).
_sufficient_samples(v̂, ρ, ε, δ) = ceil(Int, (2v̂ + 2ε / 3) / ε^2 * (log(2max(ρ, v̂) / v̂) + log(1 / δ)))

# Bound on the supremum deviation of Theorem 1 from the c-MCERA `mcera` and the
# empirical wimpy variance W of a sample of size r.
function _sd_bound(mcera, W, v̂, r, c, δ)
    L = log(4 / δ)
    R̃ = mcera + sqrt(4W * L / (c * r))
    R = R̃ + L / r + sqrt((L / r)^2 + 2L * R̃ / r)
    return 2R + sqrt(2L * (v̂ + 4R) / r) + L / (3r)
end

"""
    temporal_betweenness_mantra(g::OrderedEdgeList, criterion = MinimumHops();
                                ε = 0.01, δ = 0.1, estimator = :pairs, mc_trials = 25,
                                max_samples = typemax(Int), rng = Random.default_rng())

Estimate the *normalized* temporal betweenness `b(v) = B(v) / (n (n - 1))` of all
nodes up to an absolute error `ε` with probability at least `1 - δ`, with the MANTRA
progressive sampling algorithm (Cruciani, 2024). `criterion` is `MinimumHops()`
(shortest paths), `ShortestForemost()` or `PrefixForemost()`; for these criteria the
optimal walks are paths and the estimates of each sample lie in `[0, 1]`.

`estimator = :pairs` samples pairs of nodes (the estimator of ONBRA), `:sources`
samples nodes and uses all their optimal paths. After a bootstrap phase that bounds
the variance of the estimators and the average number of internal nodes of the
optimal paths, the sample grows by a factor 1.2 until the bound on the maximum
error obtained from `mc_trials` Monte Carlo empirical Rademacher averages is at most
`ε`, or a sufficient sample size is reached. Samples are evaluated in parallel.

The exact [`temporal_betweenness`](@ref) costs about as much as `n` samples, so
sampling pays off when the required number of samples (which grows as `1/ε²`) is
much smaller than the number of nodes.

Returns a named tuple `(estimates, error_bound, samples)`.

# References

- A. Cruciani. *MANTRA: temporal betweenness centrality approximation through sampling.* ECML PKDD, 2024. [DOI](https://doi.org/10.1007/978-3-031-70341-6_8), [arXiv](https://arxiv.org/abs/2304.08356)
"""
function temporal_betweenness_mantra(g::OrderedEdgeList, criterion::DistanceType=MinimumHops();
                                     ε::Real=0.01, δ::Real=0.1, estimator::Symbol=:pairs, mc_trials::Integer=25,
                                     max_samples::Integer=typemax(Int), rng::AbstractRNG=Random.default_rng())
    criterion isa Union{MinimumHops,ShortestForemost,PrefixForemost} ||
        throw(ArgumentError("MANTRA supports MinimumHops(), ShortestForemost() and PrefixForemost()"))
    estimator in (:pairs, :sources) || throw(ArgumentError("estimator must be :pairs or :sources"))
    (0 < ε < 1 && 0 < δ < 1) || throw(ArgumentError("ε and δ must be in (0, 1)"))
    n = num_nodes(g)
    n < 3 && return (estimates=zeros(Float64, n), error_bound=0.0, samples=0)
    sm = _sampler(g, criterion)
    c = Int(mc_trials)

    # bootstrap: bounds on the variance of the estimators and on ρ = Σ_v b(v)
    s0 = min(max_samples, max(2, ceil(Int, log(1 / δ) / ε)))
    sumsq = zeros(Float64, n)
    ρ = 0.0
    for _ in 1:s0
        _sample!(sm, _draw(rng, n, estimator)...)
        for v in sm.touched
            sumsq[v] += sm.x[v]^2
            ρ += sm.x[v]
        end
    end
    v̂ = max(_variance_bound(maximum(sumsq) / s0, s0, δ / 4), eps())
    ρ = max(ρ / s0, v̂)
    ω = min(max_samples, _sufficient_samples(v̂, ρ, ε, δ / 4))

    # first sample size: smallest s with the bound of Theorem 1 for R = 0 at most ε
    first_bound(s) = sqrt(2log(8 / δ) * v̂ / s) + log(8 / δ) / (3s)
    lo, hi = s0, max(s0, ω)
    while lo < hi
        mid = (lo + hi) ÷ 2
        first_bound(mid) <= ε ? (hi = mid) : (lo = mid + 1)
    end
    target = lo

    B = zeros(Float64, n)
    W = zeros(Float64, n)
    M = zeros(Float64, c, n)       # Σ_i λ_{j,i} x_i(v) for each Monte Carlo trial j
    r = 0
    ξ = 1.0
    iteration = 0
    while true
        iteration += 1
        # draw the new samples and Rademacher signs (serially, for reproducibility),
        # then evaluate them in parallel with one accumulator per chunk
        k = target - r
        draws = [_draw(rng, n, estimator) for _ in 1:k]
        signs = rand(rng, Bool, c, k)
        nchunks = min(k, 4 * Threads.nthreads())
        parts = _unwrap_task_errors() do
          tmap(chunks(1:k; n=nchunks)) do idx
            local_sm = _sampler(g, criterion)
            b, w, m = zeros(n), zeros(n), zeros(c, n)
            for i in idx
                _sample!(local_sm, draws[i]...)
                @inbounds for v in local_sm.touched
                    xv = local_sm.x[v]
                    b[v] += xv
                    w[v] += xv^2
                    for j in 1:c
                        m[j, v] += signs[j, i] ? xv : -xv
                    end
                end
            end
            (b, w, m)
          end
        end
        for (b, w, m) in parts
            B .+= b
            W .+= w
            M .+= m
        end
        r = target
        δi = δ / 2^(iteration + 1)
        mcera = sum(max(0.0, maximum(view(M, j, :))) for j in 1:c) / (c * r)
        wimpy = maximum(W) / r
        ξ = _sd_bound(mcera, wimpy, _variance_bound(wimpy, r, δi), r, c, δi)
        (ξ <= ε || r >= ω) && break
        target = max(r + 1, ceil(Int, 1.2r))
    end
    return (estimates=B ./ r, error_bound=min(ξ, ε), samples=r)
end

# hop distances of the optimal paths from s (typemax for unreachable nodes)
function _hop_distances!(d::Vector{Int}, sm::_Sampler, s::Int)
    n = sm.n
    fill!(d, typemax(Int))
    d[s] = 0
    if sm.criterion isa PrefixForemost
        p = sm.pfm
        _pfm_forward!(p, s)
        @inbounds for i in p.byarr  # fewest edges among the prefix foremost paths
            e = p.es[i]
            if _pfm_edge(p, e, s) && d[e.u] != typemax(Int)
                d[e.v] = min(d[e.v], d[e.u] + 1)
            end
        end
    else
        w, wg, dt = sm.work, sm.wg, sm.criterion
        _walk_forward!(w, wg, s, dt)
        _walk_targets!(w, wg, s, dt)
        @inbounds for v in 1:n
            v == s && continue
            c = w.cstar[v]
            c == _bw_inf(typeof(c)) && continue
            d[v] = dt isa MinimumHops ? c : c[2]
        end
    end
    return d
end

"""
    temporal_distance_statistics(g::OrderedEdgeList, criterion = MinimumHops();
                                 samples = nothing, τ = 0.9, rng = Random.default_rng())

Statistics of the lengths (numbers of edges) of the optimal temporal paths of type
`criterion` (`MinimumHops()`, `ShortestForemost()` or `PrefixForemost()`, for which the
fewest edges among the prefix foremost paths is used): a named tuple with

- `diameter`: the largest length,
- `effective_diameter`: the smallest `h` such that a fraction `τ` of the connected
  pairs are at distance at most `h`,
- `connectivity_rate`: the fraction of ordered pairs of distinct nodes that are
  temporally connected,
- `average_distance`: the average length over the connected pairs.

By default all nodes are used as sources (exact values, in parallel). With
`samples = r` the statistics are estimated from `r` random sources as in Algorithm 2
of Cruciani (2024); `r = Θ(log(n) / ε²)` sources give the connectivity rate within
`ε` with high probability (the diameter estimate is a lower bound).

# References

- A. Cruciani. *MANTRA: temporal betweenness centrality approximation through sampling.* ECML PKDD, 2024. [DOI](https://doi.org/10.1007/978-3-031-70341-6_8), [arXiv](https://arxiv.org/abs/2304.08356)
"""
function temporal_distance_statistics(g::OrderedEdgeList, criterion::DistanceType=MinimumHops();
                                      samples::Union{Nothing,Integer}=nothing, τ::Real=0.9,
                                      rng::AbstractRNG=Random.default_rng())
    criterion isa Union{MinimumHops,ShortestForemost,PrefixForemost} ||
        throw(ArgumentError("supported criteria: MinimumHops(), ShortestForemost() and PrefixForemost()"))
    0 < τ <= 1 || throw(ArgumentError("τ must be in (0, 1]"))
    n = num_nodes(g)
    sources = samples === nothing ? collect(1:n) : rand(rng, 1:n, samples)
    isempty(sources) && return (diameter=0, effective_diameter=0, connectivity_rate=0.0, average_distance=0.0)
    counts = tmapreduce((a, b) -> (length(a) < length(b) && ((a, b) = (b, a)); a[1:length(b)] .+= b; a),
                        chunks(sources; n=min(length(sources), 4 * Threads.nthreads()))) do ss
        sm = _sampler(g, criterion)
        d = Vector{Int}(undef, n)
        hist = zeros(Float64, 1)   # hist[h] = pairs at distance h
        for s in ss
            _hop_distances!(d, sm, s)
            for v in 1:n
                (v == s || d[v] == typemax(Int)) && continue
                length(hist) < d[v] && append!(hist, zeros(d[v] - length(hist)))
                hist[d[v]] += 1
            end
        end
        hist
    end
    scale = n / length(sources)
    reach = cumsum(counts) .* scale
    total = isempty(reach) ? 0.0 : reach[end]
    total == 0 && return (diameter=0, effective_diameter=0, connectivity_rate=0.0, average_distance=0.0)
    D = findlast(>(0), counts)
    return (diameter=D, effective_diameter=findfirst(r -> r / total >= τ, reach),
            connectivity_rate=total / (n * (n - 1)),
            average_distance=sum(h * counts[h] for h in eachindex(counts)) * scale / total)
end

# ω* of Eq. (4) of Zhang et al. (the bound of Riondato and Upfal on the empirical
# Rademacher average): min over c > 0 of (1/c) log Σ_v exp(c² ‖v‖² / (2 r²)), for the
# squared norms `sq` of the vectors of the sample values of the nodes. Zero vectors
# count once; the other vectors are counted as distinct, which can only increase the
# bound.
function _abra_omega(sq::Vector{Float64}, r::Int)
    a = [x / (2r^2) for x in sq if x > 0]
    any(<=(0), sq) && push!(a, 0.0)
    isempty(a) && return 0.0
    f(logc) = (c = exp(logc); m = maximum(a) * c^2; (m + log(sum(exp(x * c^2 - m) for x in a))) / c)
    lo, hi = -30.0, 30.0  # golden section search on log c (the function is unimodal)
    φ = (sqrt(5) - 1) / 2
    x1, x2 = hi - φ * (hi - lo), lo + φ * (hi - lo)
    f1, f2 = f(x1), f(x2)
    for _ in 1:200
        if f1 < f2
            hi, x2, f2 = x2, x1, f1
            x1 = hi - φ * (hi - lo)
            f1 = f(x1)
        else
            lo, x1, f1 = x1, x2, f2
            x2 = lo + φ * (hi - lo)
            f2 = f(x2)
        end
    end
    return min(f1, f2)
end

# upper bound on the maximum deviation, Eqs. (2) and (3) of Zhang et al.
function _abra_bound(ω, r, η)
    L = log(2 / η)
    α = L / (L + sqrt((2r * ω + L) * L))
    return ω / (1 - α) + L / (2r * α * (1 - α)) + sqrt(L / (2r))
end

"""
    temporal_betweenness_atbc(g::OrderedEdgeList, criterion = MinimumHops();
                              ε = 0.01, δ = 0.1, max_samples = typemax(Int),
                              rng = Random.default_rng())

Estimate the *normalized* temporal betweenness `B(v) / (n (n - 1))` of all nodes up to
an absolute error `ε` with probability at least `1 - δ`, with the ATBC progressive
sampling algorithm of Zhang et al. (WWW 2024): pairs of nodes are sampled until the
bound of Riondato and Upfal (ABRA) on the maximum deviation, computed from the
empirical Rademacher average of the sample, is at most `ε`. The paper considers
shortest (`MinimumHops()`), earliest (`EarliestArrival()`), fastest (`Fastest()`) and
shortest earliest (`ShortestForemost()`) paths; strict and non-strict paths correspond
to positive and zero transition times. The sample grows by a factor 1.2 between
iterations (the paper uses the sample size schedule of ABRA). Samples are evaluated in
parallel.

The guarantee assumes that every sample value `σ_sz(v) / σ_sz` is at most 1, which
holds for shortest walks; for earliest and fastest walks that visit a node several
times (counted with multiplicity) it can be slightly larger, in which case a warning
is issued.

Returns a named tuple `(estimates, error_bound, samples)`.

# References

- T. Zhang, Y. Gao, J. Zhao, L. Chen, L. Jin, Z. Yang, B. Cao, and J. Fan. *Efficient exact and approximate betweenness centrality computation for temporal graphs.* The Web Conference (WWW), 2024. [DOI](https://doi.org/10.1145/3589334.3645438)
- M. Riondato and E. Upfal. *ABRA: approximating betweenness centrality in static and dynamic graphs with Rademacher averages.* ACM Transactions on Knowledge Discovery from Data 12(5), 2018. [DOI](https://doi.org/10.1145/3208351)
"""
function temporal_betweenness_atbc(g::OrderedEdgeList, criterion::DistanceType=MinimumHops();
                                   ε::Real=0.01, δ::Real=0.1, max_samples::Integer=typemax(Int),
                                   rng::AbstractRNG=Random.default_rng())
    criterion isa Union{MinimumHops,EarliestArrival,Fastest,ShortestForemost} ||
        throw(ArgumentError("ATBC supports MinimumHops(), EarliestArrival(), Fastest() and ShortestForemost()"))
    (0 < ε < 1 && 0 < δ < 1) || throw(ArgumentError("ε and δ must be in (0, 1)"))
    n = num_nodes(g)
    n < 3 && return (estimates=zeros(Float64, n), error_bound=0.0, samples=0)
    _sampler(g, criterion)  # validate the graph (e.g., cycles of instantaneous edges)
    target = min(max_samples, ceil(Int, (1 + 8ε + sqrt(1 + 16ε)) * log(2 / δ) / (4ε^2)))
    B = zeros(Float64, n)
    W = zeros(Float64, n)
    r = 0
    ub = Inf
    iteration = 0
    largest = 0.0
    while true
        iteration += 1
        k = target - r
        draws = [_draw(rng, n, :pairs) for _ in 1:k]
        parts = _unwrap_task_errors() do
            tmap(chunks(1:k; n=min(k, 4 * Threads.nthreads()))) do idx
                sm = _sampler(g, criterion)
                b, w, mx = zeros(n), zeros(n), 0.0
                for i in idx
                    _sample!(sm, draws[i]...)
                    @inbounds for v in sm.touched
                        xv = sm.x[v]
                        b[v] += xv
                        w[v] += xv^2
                        mx = max(mx, xv)
                    end
                end
                (b, w, mx)
            end
        end
        for (b, w, mx) in parts
            B .+= b
            W .+= w
            largest = max(largest, mx)
        end
        r = target
        ub = _abra_bound(_abra_omega(W, r), r, δ / 2^iteration)
        (ub <= ε || r >= max_samples) && break
        target = min(max_samples, max(r + 1, ceil(Int, 1.2r)))
    end
    largest > 1 + 1e-9 && @warn "Some sample values exceed 1 (walks visiting a node several times), " *
                                "the error bound is not guaranteed."
    return (estimates=B ./ r, error_bound=ub, samples=r)
end
