"""
    temporal_katz_centrality(g::OrderedEdgeList, beta)

Temporal Katz centrality (Béres et al., 2018): for every node the sum of `beta^k`
over all temporal walks of length `k` ending at the node. Consecutive edges of a walk
have strictly increasing time stamps; transition times are ignored.

# References

- F. Béres, R. Pálovics, A. Oláh, and A. A. Benczúr. *Temporal walk based centrality metric for graph streams.* Applied Network Science 3, 2018. [DOI](https://doi.org/10.1007/s41109-018-0080-5)
"""
function temporal_katz_centrality(g::OrderedEdgeList{V,T}, beta::Real) where {V,T}
    n = num_nodes(g)
    r = zeros(Float64, n)
    base = zeros(Float64, n)       # value of r before time `last[v]`
    last = fill(typemin(T), n)
    @inbounds for e in edges(g)
        u, v, t = Int(e.u), Int(e.v), e.t
        before = last[u] == t ? base[u] : r[u]
        if last[v] != t
            base[v] = r[v]
            last[v] = t
        end
        r[v] += beta * (1 + before)
    end
    return r
end

"""
    temporal_pagerank(g::OrderedEdgeList, alpha, beta, gamma)

Temporal PageRank (Rozenshtein and Gionis, 2016) as implemented in TGLib, with
restart probability `1 - alpha`, transition probability `beta` and decay `gamma`.

# References

- P. Rozenshtein and A. Gionis. *Temporal PageRank.* ECML PKDD, 2016. [DOI](https://doi.org/10.1007/978-3-319-46227-1_42)
"""
function temporal_pagerank(g::OrderedEdgeList, alpha::Real, beta::Real, gamma::Real)
    n = num_nodes(g)
    m = num_edges(g)
    degrees = zeros(Float64, n)
    for e in edges(g)
        degrees[e.u] += 1
    end
    degrees ./= m
    beta == 1 && (beta = 0.0)
    r = zeros(Float64, n)
    s = zeros(Float64, n)
    c1 = (1.0 - alpha) * alpha
    c2 = (2.0 - alpha) * alpha * (1.0 - beta)
    @inbounds for e in edges(g)
        u, v = Int(e.u), Int(e.v)
        x = c1 * degrees[u]
        su = s[u]
        r[u] = r[u] * gamma + x
        r[v] = r[v] * gamma + su + x
        s[v] += su * c2
        s[u] *= beta
    end
    return r
end

# Weighted walk counts per node as (time, weight) pairs, appended in scan order.
# Only the last entry of a node can still change, and queries exclude exactly the
# entry of the current time stamp, so every edge is processed in constant time
# (TGLib iterates over all entries of the node for every edge).
struct WalkCounts{T}
    keys::Vector{Vector{T}}
    vals::Vector{Vector{Float64}}
    total::Vector{Float64}
    base::Vector{Float64}   # total without the last entry
end
WalkCounts{T}(n) where {T} = WalkCounts{T}([T[] for _ in 1:n], [Float64[] for _ in 1:n], zeros(n), zeros(n))

# sum of the weights of u, excluding the entry with time `key`
@inline function _sum_excluding(w::WalkCounts, u, key)
    ks = @inbounds w.keys[u]
    return (!isempty(ks) && @inbounds(ks[end]) == key) ? @inbounds(w.base[u]) : @inbounds(w.total[u])
end

@inline function _add!(w::WalkCounts, v, key, x)
    @inbounds begin
        ks = w.keys[v]
        if !isempty(ks) && ks[end] == key
            w.vals[v][end] += x
        else
            push!(ks, key)
            push!(w.vals[v], x)
            w.base[v] = w.total[v]
        end
        w.total[v] += x
    end
    return nothing
end

"""
    temporal_walk_centrality(g::OrderedEdgeList, alpha, beta)

Temporal walk centrality (Oettershagen, Mutzel and Kriege, 2022): for every node the
weighted number of temporal walks passing through it, where incoming walks are
weighted by `alpha` and outgoing walks by `beta` per edge. Consecutive edges of a
walk have strictly increasing time stamps. Runs in O(m + n) time.

# References

- L. Oettershagen, P. Mutzel, and N. M. Kriege. *Temporal walk centrality: ranking nodes in evolving networks.* The Web Conference (WWW), 2022. [DOI](https://doi.org/10.1145/3485447.3512210), [arXiv](https://arxiv.org/abs/2202.03706)
"""
function temporal_walk_centrality(g::OrderedEdgeList{V,T}, alpha::Real, beta::Real) where {V,T}
    n = num_nodes(g)
    es = edges(g)
    inw = WalkCounts{T}(n)     # walks ending at v with last edge at time t
    @inbounds for e in es
        _add!(inw, Int(e.v), e.t, 1 + alpha * _sum_excluding(inw, Int(e.u), e.t))
    end
    outw = WalkCounts{T}(n)    # walks starting at u with first edge at time t
    @inbounds for i in length(es):-1:1
        e = es[i]
        _add!(outw, Int(e.u), e.t, 1 + beta * _sum_excluding(outw, Int(e.v), e.t))
    end
    centrality = zeros(Float64, n)
    @inbounds for v in 1:n
        ik, iv = inw.keys[v], inw.vals[v]      # increasing time
        ok, ov = outw.keys[v], outw.vals[v]    # decreasing time
        insum = 0.0
        i = 1
        c = 0.0
        for j in length(ok):-1:1
            while i <= length(ik) && ik[i] < ok[j]
                insum += iv[i]
                i += 1
            end
            c += insum * ov[j]
        end
        centrality[v] = c
    end
    return centrality
end
