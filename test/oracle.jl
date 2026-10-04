# Brute-force reference implementations used to validate the algorithms.

const TG = TemporalGraphs

usable(e, ti) = ti[1] <= e.t && e.t + e.tt <= ti[2]

# All temporal distances from s by exhaustive search over the graph of edges
# (e -> f iff f leaves the head of e not before e arrives).
function oracle_distances(g::OrderedEdgeList, s, ti=time_interval(g))
    n = num_nodes(g)
    es = [e for e in edges(g) if usable(e, ti)]
    m = length(es)
    succ = [[j for j in 1:m if es[j].u == es[i].v && es[j].t >= es[i].t + es[i].tt] for i in 1:m]
    ea = fill(INF, n); ea[s] = 0
    fastest = fill(INF, n); fastest[s] = 0
    ld = fill(typemin(Int), n); ld[s] = ti[2]
    for i in 1:m
        es[i].u == s || continue
        seen = falses(m)
        stack = [i]
        seen[i] = true
        while !isempty(stack)
            j = pop!(stack)
            f = es[j]
            v = f.v
            v != s && (ea[v] = min(ea[v], f.t + f.tt))
            fastest[v] = min(fastest[v], f.t + f.tt - es[i].t)
            ld[v] = max(ld[v], es[i].t)
            for k in succ[j]
                seen[k] || (seen[k] = true; push!(stack, k))
            end
        end
    end
    # hops and transition times: Bellman-Ford on the edge graph
    function cost(w)
        c = fill(INF, m)
        for i in 1:m
            es[i].u == s && (c[i] = w(es[i]))
        end
        for _ in 1:m, j in 1:m
            c[j] == INF && continue
            for k in succ[j]
                c[k] = min(c[k], c[j] + w(es[k]))
            end
        end
        d = fill(INF, n); d[s] = 0
        for j in 1:m
            d[es[j].v] = min(d[es[j].v], c[j])
        end
        return d
    end
    return (ea=ea, fastest=fastest, ld=ld, hops=cost(e -> 1), mtt=cost(e -> e.tt))
end

# All temporal walks with strictly increasing time stamps.
function all_walks(g::OrderedEdgeList)
    es = edges(g)
    walks = Vector{Vector{Int}}()
    function extend(w)
        push!(walks, copy(w))
        last = es[w[end]]
        for j in eachindex(es)
            if es[j].u == last.v && es[j].t > last.t
                push!(w, j)
                extend(w)
                pop!(w)
            end
        end
    end
    for i in eachindex(es)
        extend([i])
    end
    return walks
end

function oracle_katz(g, beta)
    r = zeros(num_nodes(g))
    for w in all_walks(g)
        r[edges(g)[w[end]].v] += beta^length(w)
    end
    return r
end

function oracle_walk_centrality(g, alpha, beta)
    es = edges(g)
    walks = all_walks(g)
    c = zeros(num_nodes(g))
    for w1 in walks, w2 in walks
        last, first = es[w1[end]], es[w2[1]]
        if last.v == first.u && last.t + 1 <= first.t
            c[last.v] += alpha^(length(w1) - 1) * beta^(length(w2) - 1)
        end
    end
    return c
end

# Core numbers by definition: the largest k such that v is in the k-core.
function oracle_kcores(pairs, n)
    adj = [Set{Int}() for _ in 1:n]
    for (u, v) in pairs
        u == v && continue
        push!(adj[u], v); push!(adj[v], u)
    end
    core = zeros(Int, n)
    for k in 1:n
        alive = trues(n)
        changed = true
        while changed
            changed = false
            for v in 1:n
                if alive[v] && count(u -> alive[u], adj[v]) < k
                    alive[v] = false
                    changed = true
                end
            end
        end
        for v in 1:n
            alive[v] && (core[v] = k)
        end
    end
    return core
end

# Brandes-free betweenness of the directed line graph: counts of shortest paths.
function oracle_betweenness(dlg)
    N = num_nodes(dlg)
    d = fill(-1, N, N)
    σ = zeros(N, N)
    for s in 1:N
        d[s, s] = 0; σ[s, s] = 1
        frontier = [s]
        while !isempty(frontier)
            next = Int[]
            for u in frontier, v in out_edges(dlg, u)
                if d[s, v] < 0
                    d[s, v] = d[s, u] + 1
                    push!(next, v)
                end
                d[s, v] == d[s, u] + 1 && (σ[s, v] += σ[s, u])
            end
            frontier = next
        end
    end
    bc = zeros(N)
    for s in 1:N, t in 1:N, w in 1:N
        (s == t || w == s || w == t || d[s, t] < 0 || d[s, w] < 0 || d[w, t] < 0) && continue
        d[s, w] + d[w, t] == d[s, t] && (bc[w] += σ[s, w] * σ[w, t] / σ[s, t])
    end
    return bc
end

# All restless walks (β = maximum waiting time) with at most `maxlen` edges, as
# vectors of edge indices (with transition times 0 cycles at one time stamp would
# give infinitely many walks).
function all_restless_walks(es, β; maxlen=typemax(Int))
    walks = Vector{Vector{Int}}()
    function extend(w)
        push!(walks, copy(w))
        length(w) >= maxlen && return
        last = es[w[end]]
        arr = last.t + last.tt
        for j in eachindex(es)
            if es[j].u == last.v && arr <= es[j].t <= arr + β
                push!(w, j)
                extend(w)
                pop!(w)
            end
        end
    end
    for i in eachindex(es)
        extend([i])
    end
    return walks
end

# Temporal betweenness by definition (Brunelli et al.): for every pair s ≠ t the
# optimal s-t walks, and for every node u ≠ s, t the fraction of them through u,
# counting every edge entering u. Latest departure walks can revisit their target;
# as in Fact 3 of Brunelli et al., these earlier visits of t also count for t.
function oracle_temporal_betweenness(g, criterion, β; maxlen=typemax(Int))
    es = edges(g)
    n = num_nodes(g)
    walks = all_restless_walks(es, β; maxlen=maxlen)
    dep(w) = es[w[1]].t
    arr(w) = es[w[end]].t + es[w[end]].tt
    key(w) = criterion isa MinimumHops ? length(w) :
             criterion isa EarliestArrival ? arr(w) :
             criterion isa LatestDeparture ? -dep(w) :
             criterion isa Fastest ? arr(w) - dep(w) :
             criterion isa ShortestForemost ? (arr(w), length(w)) :
             criterion isa ShortestLatest ? (-dep(w), length(w)) :
             (arr(w) - dep(w), length(w))
    bc = zeros(n)
    for s in 1:n, t in 1:n
        s == t && continue
        W = [w for w in walks if es[w[1]].u == s && es[w[end]].v == t]
        isempty(W) && continue
        best = minimum(key, W)
        opt = [w for w in W if key(w) == best]
        for w in opt, j in 1:length(w)-1
            u = es[w[j]].v
            (u != s && (u != t || criterion isa LatestDeparture)) && (bc[u] += 1 / length(opt))
        end
    end
    return bc
end

# Prefix foremost betweenness by definition (Buß et al.): strict paths from s to t
# whose every prefix is a foremost path, i.e., every node is reached at its earliest
# arrival time from s.
function oracle_prefix_foremost_betweenness(g)
    es = edges(g)
    n = num_nodes(g)
    bc = zeros(n)
    for s in 1:n
        ea = fill(typemax(Int), n)
        for w in all_restless_walks(es, 10^9)
            es[w[1]].u == s || continue
            for e in w
                v = es[e].v
                v != s && (ea[v] = min(ea[v], es[e].t + es[e].tt))
            end
        end
        paths = Vector{Vector{Int}}()
        function extend(p, visited)
            push!(paths, copy(p))
            last = es[p[end]]
            for j in eachindex(es)
                f = es[j]
                if f.u == last.v && f.t >= last.t + last.tt && !(f.v in visited) && f.t + f.tt == ea[f.v]
                    push!(p, j); push!(visited, f.v)
                    extend(p, visited)
                    pop!(p); delete!(visited, f.v)
                end
            end
        end
        for i in eachindex(es)
            e = es[i]
            (e.u == s && e.v != s && e.t + e.tt == ea[e.v]) && extend([i], Set([s, e.v]))
        end
        for t in 1:n
            t == s && continue
            P = [p for p in paths if es[p[end]].v == t]
            for p in P, j in 1:length(p)-1
                bc[es[p[j]].v] += 1 / length(P)
            end
        end
    end
    return bc
end

# Pass-through degree by definition: pairs (v, w) of nodes ≠ u with an edge v → u
# followed in time by an edge u → w.
function oracle_pass_through_degree(g)
    es = edges(g)
    n = num_nodes(g)
    return [sqrt(count(((v, w),) -> any(e1 -> e1.u == v && e1.v == u &&
                 any(e2 -> e2.u == u && e2.v == w && e2.t >= e1.t + e1.tt, es), es),
                 [(v, w) for v in 1:n, w in 1:n if v != u && w != u])) for u in 1:n]
end

# Ego betweenness by definition: the betweenness of v in the subgraph induced by v
# and its neighbors.
function oracle_ego_betweenness(g, criterion)
    es = edges(g)
    n = num_nodes(g)
    map(1:n) do v
        nodes = sort!(unique!([v; [e.u for e in es if e.v == v && e.u != v]; [e.v for e in es if e.u == v && e.v != v]]))
        id = Dict(x => i for (i, x) in enumerate(nodes))
        sub = [TemporalEdge(id[e.u], id[e.v], e.t, e.tt) for e in es if haskey(id, e.u) && haskey(id, e.v)]
        sg = OrderedEdgeList(length(nodes), sub)
        b = criterion isa PrefixForemost ? oracle_prefix_foremost_betweenness(sg) :
            oracle_temporal_betweenness(sg, criterion, 10^9)
        b[id[v]]
    end
end

# Closeness integrated over the departure time by definition: for every departure
# time τ the earliest arrival over all walks leaving s at or after τ, integrated
# exactly between consecutive departure times.
function oracle_harmonic_closeness(g, ti=time_interval(g))
    A, B = ti
    es = [e for e in edges(g) if A <= e.t && e.t + e.tt <= B]
    n = num_nodes(g)
    walks = all_restless_walks(es, 10^9)
    c = zeros(n)
    for s in 1:n, v in 1:n
        s == v && continue
        J = [(es[w[1]].t, es[w[end]].t + es[w[end]].tt) for w in walks if es[w[1]].u == s && es[w[end]].v == v]
        isempty(J) && continue
        R = sort!(unique!(first.(J)))
        prev = A
        for r in R
            a = minimum(last(j) for j in J if first(j) >= r)
            r > prev && (c[s] += log((a - prev) / (a - r)))
            prev = max(prev, r)
        end
    end
    return c ./ ((n - 1) * (B - A))
end

# Hop lengths of the optimal paths by enumeration: shortest walks, shortest among the
# foremost walks, shortest among the prefix foremost walks (every prefix ends at the
# earliest arrival time of its last node).
function oracle_hop_distances(g, criterion)
    es = edges(g)
    n = num_nodes(g)
    walks = all_restless_walks(es, 10^9)
    d = fill(typemax(Int), n, n)
    for s in 1:n
        d[s, s] = 0
        W = [w for w in walks if es[w[1]].u == s]
        arr(w) = es[w[end]].t + es[w[end]].tt
        ea = fill(typemax(Int), n)
        for w in W, k in 1:length(w)
            v = es[w[k]].v
            v != s && (ea[v] = min(ea[v], es[w[k]].t + es[w[k]].tt))
        end
        for t in 1:n
            t == s && continue
            Wt = [w for w in W if es[w[end]].v == t && !any(es[e].v == s for e in w)]
            if criterion isa ShortestForemost
                Wt = [w for w in Wt if arr(w) == ea[t]]
            elseif criterion isa PrefixForemost
                Wt = [w for w in Wt if all(es[w[k]].t + es[w[k]].tt == ea[es[w[k]].v] for k in eachindex(w))]
            end
            isempty(Wt) || (d[s, t] = minimum(length, Wt))
        end
    end
    return d
end

# δ-temporal motifs with 3 edges by enumeration of all triples of edges, in the 6×6
# layout of Paranjape et al.: ties are ordered as in the edge list unless `strict`.
function oracle_motif_counts(g, δ; strict=true)
    es = edges(g)
    m = length(es)
    M = zeros(Int, 6, 6)
    order = ((1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2))
    for i in 1:m, j in i+1:m, k in j+1:m
        a, b, c = es[i], es[j], es[k]
        any(e -> e.u == e.v, (a, b, c)) && continue
        c.t - a.t <= δ || continue
        strict && !(a.t < b.t < c.t) && continue
        nodes = unique([a.u, a.v, b.u, b.v, c.u, c.v])
        length(nodes) <= 3 || continue
        lab(x) = x == a.u ? 1 : x == a.v ? 2 : 3
        pos(e) = findfirst(==((lab(e.u), lab(e.v))), order)
        M[7-pos(b), pos(c)] += 1
    end
    return M
end

# Instances of an arbitrary motif by enumeration of all l-subsets of edges.
function oracle_motif_count(g, motif, δ; strict=true)
    es = [e for e in edges(g)]
    m = length(es)
    l = length(motif)
    count = 0
    function rec(chosen, from)
        if length(chosen) == l
            ts = [es[i].t for i in chosen]
            ts[end] - ts[1] <= δ || return
            strict && any(ts[r] >= ts[r+1] for r in 1:l-1) && return
            f = Dict{Int,Int}()
            for (r, i) in enumerate(chosen)
                x, y = motif[r]
                for (p, q) in ((x, es[i].u), (y, es[i].v))
                    get!(f, p, q) == q || return
                end
            end
            length(unique(values(f))) == length(f) && (count += 1)
            return
        end
        for i in from:m
            push!(chosen, i)
            rec(chosen, i + 1)
            pop!(chosen)
        end
    end
    rec(Int[], 1)
    return count
end

# Restless walk distances by enumeration of the walks (β = maximum waiting time).
function oracle_restless_distances(g, s, dt, β)
    es = edges(g)
    n = num_nodes(g)
    W = [w for w in all_restless_walks(es, β; maxlen=length(es)) if es[w[1]].u == s]
    D = distance_eltype(g, dt)
    b = time_interval(g)[2]
    val(w) = dt isa EarliestArrival ? es[w[end]].t + es[w[end]].tt :
             dt isa LatestDeparture ? es[w[1]].t :
             dt isa Fastest ? es[w[end]].t + es[w[end]].tt - es[w[1]].t :
             dt isa MinimumHops ? length(w) : sum(es[i].tt for i in w)
    dist = fill(dt isa LatestDeparture ? typemin(D) : typemax(D), n)
    for w in W
        v = es[w[end]].v
        dist[v] = dt isa LatestDeparture ? max(dist[v], val(w)) : min(dist[v], val(w))
    end
    dist[s] = dt isa LatestDeparture ? b : zero(D)
    return dist
end

# Temporal walks (or paths) from s to z as vectors of edge indices.
function oracle_sz_walks(es, s, z, β, n; paths=true)
    W = all_restless_walks(es, β; maxlen=paths ? n : length(es) + 1)
    return [w for w in W if es[w[1]].u == s && es[w[end]].v == z &&
                            (!paths || allunique([es[w[1]].u; [es[i].v for i in w]]))]
end

# Largest number of pairwise compatible walks.
function oracle_packing(P, conflict)
    best = 0
    function rec(chosen, from)
        best = max(best, length(chosen))
        for i in from:length(P)
            all(!conflict(P[i], P[j]) for j in chosen) || continue
            push!(chosen, i)
            rec(chosen, i + 1)
            pop!(chosen)
        end
    end
    rec(Int[], 1)
    return best
end

# Minimum capacity of a set of edges meeting all walks of P.
function oracle_min_cut(P, m, cap)
    best = Inf
    for mask in 0:(2^m-1)
        S = [i for i in 1:m if (mask >> (i - 1)) & 1 == 1]
        all(w -> any(in(S), w), P) && (best = min(best, sum(cap[S]; init=0)))
    end
    return best
end

# Arrival times along a temporal spanning tree rooted at r, or false if it is not one.
function oracle_tree_arrivals(g, r, tree)
    n = num_nodes(g)
    a = time_interval(g)[1]
    ea = earliest_arrival_times(g, r)
    reach = [v != r && ea[v] < typemax(eltype(ea)) for v in 1:n]
    length(tree) == count(reach) || return false
    inc = Dict{Int,eltype(tree)}()
    for e in tree
        (haskey(inc, e.v) || e.v == r || !reach[e.v]) && return false
        inc[e.v] = e
    end
    arr = Dict{Int,Int}()
    function arrival(v)
        v == r && return a
        haskey(arr, v) && return arr[v]
        haskey(inc, v) || return -1
        arr[v] = -1
        e = inc[v]
        x = arrival(Int(e.u))
        arr[v] = (x >= 0 && x <= e.t) ? e.t + e.tt : -1
    end
    all(v -> arrival(v) >= 0, keys(inc)) || return false
    return arr
end

# Temporal graph isomorphisms by their definitions (Heeg et al.): all node bijections
# and, for :event, all bijections of the parallel edges of every pair of nodes.
function _oracle_perms(v)
    length(v) <= 1 && return [copy(v)]
    out = Vector{Vector{eltype(v)}}()
    for i in eachindex(v), p in _oracle_perms(deleteat!(copy(v), i))
        push!(out, [v[i]; p])
    end
    return out
end
_oracle_follows(e, f, β) = e.v == f.u && e.t + e.tt <= f.t <= e.t + e.tt + β
function oracle_temporal_isomorphic(g1, g2, kind, β)
    n = num_nodes(g1)
    n == num_nodes(g2) || return false
    e1, e2 = edges(g1), edges(g2)
    length(e1) == length(e2) || return false
    t1 = isempty(e1) ? 0 : minimum(e.t for e in e1)
    t2 = isempty(e2) ? 0 : minimum(e.t for e in e2)
    lab(g, t0, u, v) = (ts = sort([e.t - t0 for e in edges(g) if (e.u, e.v) == (u, v)]); kind == :aggregated ? length(ts) : ts)
    for π in _oracle_perms(collect(1:n))
        if kind == :event
            groups = Dict{Tuple{Int,Int},Vector{Int}}()
            for (i, e) in enumerate(e1)
                push!(get!(groups, (e.u, e.v), Int[]), i)
            end
            ks = collect(keys(groups))
            targets = [[j for (j, f) in enumerate(e2) if (f.u, f.v) == (π[k[1]], π[k[2]])] for k in ks]
            all(length(t) == length(groups[k]) for (k, t) in zip(ks, targets)) || continue
            for combo in Iterators.product((_oracle_perms(t) for t in targets)...)
                πE = zeros(Int, length(e1))
                for (k, c) in zip(ks, combo), (i, j) in zip(groups[k], c)
                    πE[i] = j
                end
                all(_oracle_follows(e1[i], e1[j], β) == _oracle_follows(e2[πE[i]], e2[πE[j]], β)
                    for i in eachindex(e1), j in eachindex(e1)) && return true
            end
        else
            all(lab(g1, t1, u, v) == lab(g2, t2, π[u], π[v]) for u in 1:n, v in 1:n) && return true
        end
    end
    return false
end
