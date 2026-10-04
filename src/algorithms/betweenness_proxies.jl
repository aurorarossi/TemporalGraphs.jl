# Prefix foremost betweenness, local proxies of the temporal betweenness and its
# sampling approximation. See Becker, Crescenzi, Cruciani and Kodric, "Proxying
# betweenness centrality rankings in temporal networks" (SEA 2023), Buß, Molter,
# Niedermeier and Rymar, "Algorithmic aspects of temporal betweenness" (KDD 2020), and
# Santoro and Sarpe, "ONBRA: rigorous estimation of the temporal betweenness
# centrality in temporal networks" (WWW 2022).

function _strict_edges(g::OrderedEdgeList{V,T}, ti) where {V,T}
    a, b = _interval(T, ti)
    es = [e for e in edges(g) if a <= e.t && e.t + e.tt <= b]
    any(e -> !(e.tt > 0), es) && throw(ArgumentError(
        "temporal betweenness requires positive transition times (strict temporal walks)"))
    return es
end

# Prefix foremost paths reach every node at its earliest arrival time, so the edges
# that can be used (from a node reached at time ea[u] to a node first reached when
# the edge arrives) form a DAG ordered by arrival times. Brandes' accumulation on this
# DAG takes O(M) time per source after one earliest arrival scan.
function _prefix_foremost_betweenness(g::OrderedEdgeList{V,T}, ti, ::Type{R}; parallel::Bool=true) where {V,T,R}
    es = _strict_edges(g, ti)
    n = num_nodes(g)
    n == 0 && return Float64[], false
    byarr = sortperm([e.t + e.tt for e in es])
    function accumulate(sources)
        ea = Vector{T}(undef, n)
        σ = Vector{R}(undef, n)
        δ = Vector{R}(undef, n)
        bc = zeros(Float64, n)
        overflow = false
        for s in sources
            fill!(ea, typemax(T))
            ea[s] = typemin(T)
            @inbounds for e in es  # sorted by departure, transition times are positive
                (ea[e.u] <= e.t && e.t + e.tt < ea[e.v]) && (ea[e.v] = e.t + e.tt)
            end
            fill!(σ, zero(R))
            fill!(δ, zero(R))
            σ[s] = one(R)
            @inbounds for i in byarr
                e = es[i]
                (e.v != s && ea[e.u] <= e.t && e.t + e.tt == ea[e.v]) && (σ[e.v] += σ[e.u])
            end
            @inbounds for i in Iterators.reverse(byarr)
                e = es[i]
                if e.v != s && ea[e.u] <= e.t && e.t + e.tt == ea[e.v]
                    δ[e.u] += σ[e.u] / σ[e.v] * (1 + δ[e.v])
                end
            end
            @inbounds for v in 1:n
                v != s && (bc[v] += Float64(δ[v]))
                overflow |= !isfinite(σ[v])
            end
        end
        return bc, overflow
    end
    parallel || return accumulate(1:n)
    return tmapreduce(accumulate, (x, y) -> (x[1] .+ y[1], x[2] | y[2]), chunks(1:n; n=min(n, 4 * Threads.nthreads())))
end

"""
    temporal_ego_betweenness(g::OrderedEdgeList, criterion = MinimumHops(), ti = time_interval(g);
                             β = Inf, count_type = nothing)

Temporal ego-betweenness (Becker et al., 2023): for every node `v`, the temporal
betweenness of `v` (see [`temporal_betweenness`](@ref), same criteria) in its
*ego network*, the temporal subgraph induced by `v` and its in- and out-neighbors.
A local proxy of the temporal betweenness. Nodes are processed in parallel.

# References

- R. Becker, P. Crescenzi, A. Cruciani, and B. Kodric. *Proxying betweenness centrality rankings in temporal networks.* SEA, 2023. [DOI](https://doi.org/10.4230/LIPIcs.SEA.2023.6)
"""
function temporal_ego_betweenness(g::OrderedEdgeList{V,T}, criterion::DistanceType=MinimumHops(), ti=time_interval(g);
                                  β::Real=Inf, count_type::Union{Nothing,Type{<:Real}}=nothing) where {V,T}
    R = count_type === nothing ? Float64 : count_type
    return _with_count_type(R -> _ego_betweenness(g, criterion, ti, β, R), R, count_type !== nothing)
end

function _ego_betweenness(g::OrderedEdgeList{V,T}, criterion, ti, β, ::Type{R}) where {V,T,R}
    (criterion isa _WalkCriterion || criterion isa PrefixForemost) ||
        throw(ArgumentError("temporal betweenness is not defined for $(typeof(criterion))"))
    criterion isa PrefixForemost && isfinite(β) &&
        throw(ArgumentError("prefix foremost betweenness has no waiting constraint"))
    es = _strict_edges(g, ti)
    n = num_nodes(g)
    il = IncidentLists(n, es, _interval(T, ti))
    nbrs = [Int[] for _ in 1:n]
    for e in es
        e.u == e.v && continue
        push!(nbrs[e.u], e.v)
        push!(nbrs[e.v], e.u)
    end
    foreach(x -> unique!(sort!(x)), nbrs)
    result = zeros(Float64, n)
    overflowed = zeros(Bool, n)
    @tasks for v in 1:n
        @local localid = zeros(Int, n)
        nodes = sort!(push!(copy(nbrs[v]), v))
        for (k, x) in enumerate(nodes)
            localid[x] = k
        end
        sub = TemporalEdge{V,T}[]
        for x in nodes, e in out_edges(il, x)
            localid[e.v] != 0 && push!(sub, TemporalEdge{V,T}(localid[x], localid[e.v], e.t, e.tt))
        end
        sg = OrderedEdgeList(length(nodes), sub, _interval(T, ti))
        bc, over = if criterion isa PrefixForemost
            _prefix_foremost_betweenness(sg, time_interval(sg), R; parallel=false)
        else
            _walk_betweenness(_walk_graph(sg, time_interval(sg), β), criterion, R; parallel=false)
        end
        overflowed[v] = over
        result[v] = bc[localid[v]]
        for x in nodes
            localid[x] = 0
        end
    end
    return result, any(overflowed)
end

"""
    temporal_pass_through_degree(g::OrderedEdgeList, ti = time_interval(g))

Temporal pass-through degree (Becker et al., 2023): for every node `u`, the square
root of the number of ordered pairs `(v, w)` of nodes different from `u` (possibly
`v = w`) such that some edge from `v` to `u` is followed in time by some edge from
`u` to `w`, i.e., the earliest arrival of an edge `v → u` is not later than the
latest departure of an edge `u → w`. A local proxy of the temporal betweenness
computable in `O(M log M)` time.

# References

- R. Becker, P. Crescenzi, A. Cruciani, and B. Kodric. *Proxying betweenness centrality rankings in temporal networks.* SEA, 2023. [DOI](https://doi.org/10.4230/LIPIcs.SEA.2023.6)
"""
function temporal_pass_through_degree(g::OrderedEdgeList{V,T}, ti=time_interval(g)) where {V,T}
    a, b = _interval(T, ti)
    n = num_nodes(g)
    first_arrival = Dict{Tuple{V,V},T}()
    last_departure = Dict{Tuple{V,V},T}()
    for e in edges(g)
        (a <= e.t && e.t + e.tt <= b && e.u != e.v) || continue
        key = (e.u, e.v)
        first_arrival[key] = min(get(first_arrival, key, typemax(T)), e.t + e.tt)
        last_departure[key] = max(get(last_departure, key, typemin(T)), e.t)
    end
    arrivals = [T[] for _ in 1:n]  # earliest arrivals of the static arcs entering each node
    for ((_, v), t) in first_arrival
        push!(arrivals[v], t)
    end
    foreach(sort!, arrivals)
    count = zeros(Int, n)
    for ((u, _), t) in last_departure
        count[u] += searchsortedlast(arrivals[u], t)
    end
    return sqrt.(Float64.(count))
end

"""
    temporal_betweenness_approximation(g::OrderedEdgeList, samples, criterion = MinimumHops(),
                                       ti = time_interval(g); β = Inf, η = 0.1,
                                       rng = Random.default_rng(), count_type = Float64)

Estimate the *normalized* temporal betweenness `b(v) / (n (n - 1))` of all nodes
with the ONBRA sampling scheme (Santoro and Sarpe, 2022): `samples` ordered pairs
`(s, z)` of distinct nodes are drawn uniformly at random and, for each of them, every
node other than `s` and `z` receives the fraction of optimal `s`-`z` walks through it.
Returns `(estimates, ε)` where, with probability at least `1 - η`, every estimate is
within `ε` of the exact value (empirical Bernstein bound). The criteria and `β` are
those of [`temporal_betweenness`](@ref); multiplicities are counted as there.

# References

- D. Santoro and I. Sarpe. *ONBRA: rigorous estimation of the temporal betweenness centrality in temporal networks.* The Web Conference (WWW), 2022. [DOI](https://doi.org/10.1145/3485447.3512204), [arXiv](https://arxiv.org/abs/2203.00653)
"""
function temporal_betweenness_approximation(g::OrderedEdgeList, samples::Integer, criterion::DistanceType=MinimumHops(),
                                            ti=time_interval(g); β::Real=Inf, η::Real=0.1,
                                            rng::AbstractRNG=Random.default_rng(),
                                            count_type::Type{R}=Float64) where {R<:Real}
    criterion isa _WalkCriterion ||
        throw(ArgumentError("the approximation supports the walk criteria of temporal_betweenness"))
    samples >= 2 || throw(ArgumentError("at least two samples are needed"))
    0 < η < 1 || throw(ArgumentError("η must be in (0, 1)"))
    wg = _walk_graph(g, ti, β)
    n = wg.n
    n < 2 && return zeros(Float64, n), 0.0
    pairs = Vector{Tuple{Int,Int}}(undef, samples)
    for i in 1:samples
        s = rand(rng, 1:n)
        z = rand(rng, 1:n-1)
        pairs[i] = (s, z >= s ? z + 1 : z)
    end
    sort!(pairs)  # samples with the same source share the forward pass
    groups = Vector{UnitRange{Int}}()
    i = 1
    while i <= samples
        j = i
        while j < samples && pairs[j+1][1] == pairs[i][1]
            j += 1
        end
        push!(groups, i:j)
        i = j + 1
    end
    total, total2, overflow = _unwrap_task_errors() do
        tmapreduce((x, y) -> (x[1] .+ y[1], x[2] .+ y[2], x[3] | y[3]),
                                         chunks(groups; n=min(length(groups), 4 * Threads.nthreads()))) do gs
        w = _WalkWork(wg, criterion, R)
        x = zeros(Float64, n)
        sum1 = zeros(Float64, n)
        sum2 = zeros(Float64, n)
        over = false
        for grp in gs
            s = pairs[first(grp)][1]
            over |= _walk_forward!(w, wg, s, criterion)
            for k in grp
                z = pairs[k][2]
                _walk_targets!(w, wg, s, criterion, z)
                fill!(x, 0.0)
                _walk_backward!(w, wg, s, criterion) do v, val
                    @inbounds v != z && (x[v] += Float64(val))
                end
                sum1 .+= x
                sum2 .+= x .^ 2
            end
        end
        (sum1, sum2, over)
        end
    end
    overflow && throw(OverflowError(
        "the numbers of optimal walks overflow $R; use count_type = BigFloat or Rational{BigInt}"))
    ℓ = samples
    estimates = total ./ ℓ
    variance = max.((total2 .- ℓ .* estimates .^ 2) ./ (ℓ - 1), 0.0)
    L = log(4n / η)
    ε = maximum(sqrt.(2 .* variance .* L ./ ℓ)) + 7L / (3(ℓ - 1))
    return estimates, ε
end
