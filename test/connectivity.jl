# Reachability, temporal components, restless walks, dissemination times and static expansions.

@testset "connectivity: reachability, components, restless walks and paths" begin
    TG = TemporalGraphs
    rng = MersenneTwister(3)
    for _ in 1:100
        n = rand(rng, 1:150)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:15), rand(rng, 0:3)) for _ in 1:rand(rng, 0:400)])
        ti = rand(rng, Bool) ? time_interval(g) : (rand(rng, 0:5), rand(rng, 8:20))
        R = temporal_reachability(g, ti)
        @test all(R[s, :] == (earliest_arrival_times(g, s, ti) .< typemax(Int)) .| ((1:n) .== s) for s in 1:n)
        @test is_temporally_connected(g, ti) == all(R)
    end
    # components against all subsets of nodes
    function connected(g, X, closed, uni)
        h = closed ? TG._induced(g, X, time_interval(g)) : g
        R = temporal_reachability(h)
        idx = closed ? collect(eachindex(X)) : X
        return all(uni ? (R[a, b] || R[b, a]) : (R[a, b] && R[b, a]) for a in idx, b in idx)
    end
    for _ in 1:60
        n = rand(rng, 1:7)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:8), rand(rng, 0:2)) for _ in 1:rand(rng, 0:16)])
        subsets = [[i for i in 1:n if (mask >> (i - 1)) & 1 == 1] for mask in 1:(2^n-1)]
        for uni in (false, true)
            conn = [X for X in subsets if connected(g, X, false, uni)]
            maximal = sort([X for X in conn if !any(Y -> length(Y) > length(X) && issubset(X, Y), conn)]; by=c -> (-length(c), c))
            @test temporal_connected_components(g; unilateral=uni) == maximal
            for closed in (false, true)
                L = largest_temporal_connected_component(g; closed=closed, unilateral=uni)
                @test connected(g, L, closed, uni)
                @test length(L) == maximum(length(X) for X in subsets if connected(g, X, closed, uni))
            end
        end
    end
    # {1, 2} is connected through 3 and 4 but not inside: open but not closed
    g = OrderedEdgeList(4, [(1, 3, 1, 0), (3, 2, 2, 0), (2, 4, 3, 0), (4, 1, 4, 0)])
    @test temporal_connected_components(g)[1] == [1, 2]
    @test length(largest_temporal_connected_component(g; closed=true)) == 1
    # restless walks and paths against enumeration
    for _ in 1:80
        n = rand(rng, 2:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:10), rand(rng, 0:2)) for _ in 1:rand(rng, 1:10)])
        es = edges(g)
        for dt in (EarliestArrival(), LatestDeparture(), Fastest(), MinimumHops(), MinimumTransitionTimes()), s in 1:n
            for β in (0, 1, 3)
                @test temporal_distances(g, s, dt; β=β) == oracle_restless_distances(g, s, dt, β)
            end
            @test temporal_distances(g, s, dt; β=10^6) == temporal_distances(g, s, dt)
        end
        for β in (0, 2), s in 1:n, z in 1:n
            s == z && continue
            p = restless_path(g, s, z, β)
            @test (p !== nothing) == !isempty(oracle_sz_walks(es, s, z, β, n))
            if p !== nothing
                @test p[1].u == s && p[end].v == z && allunique([p[1].u; [e.v for e in p]])
                @test all(p[k].v == p[k+1].u && p[k].t + p[k].tt <= p[k+1].t <= p[k].t + p[k].tt + β for k in 1:length(p)-1)
            end
        end
    end
    # a restless walk that needs to revisit a node: s=1, b=2, c=3, d=4, z=5
    g = OrderedEdgeList(5, [(1, 2, 1, 0), (2, 3, 2, 0), (3, 4, 3, 0), (4, 2, 4, 0), (2, 5, 5, 0)])
    @test earliest_arrival_times(g, 1)[5] == 5
    @test temporal_distances(g, 1, EarliestArrival(); β=1)[5] == 5
    @test restless_path(g, 1, 5, 1) === nothing
    @test restless_path(g, 1, 5, 4) == [TemporalEdge(1, 2, 1, 0), TemporalEdge(2, 5, 5, 0)]
    @test_throws ArgumentError temporal_distances(g, 1, EarliestArrival(); β=-1)
    @test_throws ArgumentError temporal_distances(to_incident_lists(g), 1, EarliestArrival(); β=1)
end

@testset "dissemination times, static expansions and further components" begin
    TG = TemporalGraphs
    bfs(h, srcs) = (seen = falses(Graphs.nv(h)); q = collect(srcs); seen[q] .= true;
                    while !isempty(q); u = popfirst!(q); for w in Graphs.outneighbors(h, u); seen[w] || (seen[w] = true; push!(q, w)); end; end; seen)
    rng = MersenneTwister(17)
    for _ in 1:80
        n = rand(rng, 1:7)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:8), rand(rng, 0:2)) for _ in 1:rand(rng, 0:14)])
        ti = rand(rng, Bool) ? time_interval(g) : (rand(rng, 0:3), rand(rng, 5:10))
        es = edges(g)
        a, b = ti
        A = earliest_arrival_matrix(g, ti)
        @test all(A[s, :] == earliest_arrival_times(g, s, ti) for s in 1:n)
        @test temporal_flooding_times(g, ti) == [something(temporal_flooding_time(g, s, ti), typemax(Int)) for s in 1:n]
        gt = temporal_gossip_time(g, ti)
        if is_temporally_connected(g, ti)
            cands = sort(unique([a; [e.t + e.tt for e in es if a <= e.t && e.t + e.tt <= b]]))
            @test gt == cands[findfirst(t -> is_temporally_connected(g, (a, t)), cands)] - a
        else
            @test gt === nothing
        end
        R = temporal_reachability(g, ti)
        for compressed in (true, false)
            X = vertex_expansion(g, ti; compressed=compressed)
            Rx = [i == j for i in 1:n, j in 1:n]
            for s in 1:n
                first_copy = findfirst(x -> x[1] == s, X.nodes)
                first_copy === nothing && continue
                reach = bfs(X.graph, [first_copy])
                for (i, (v, _)) in enumerate(X.nodes)
                    reach[i] && (Rx[s, v] = true)
                end
            end
            @test Rx == R
        end
        for (β, hubs) in ((Inf, false), (Inf, true), (0, false), (1, false))
            E = edge_expansion(g, ti; β=β, hubs=hubs)
            Re = [i == j for i in 1:n, j in 1:n]
            for s in 1:n
                starts = [x for (x, i) in enumerate(E.edges) if es[i].u == s]
                isempty(starts) && continue
                reach = bfs(E.graph, starts)
                for (x, i) in enumerate(E.edges)
                    reach[x] && (Re[s, es[i].v] = true)
                end
            end
            Rβ = isfinite(β) ? permutedims(reduce(hcat, [(temporal_distances(g, s, EarliestArrival(), ti; β=β) .< typemax(Int)) .| ((1:n) .== s) for s in 1:n])) : R
            @test Re == Rβ
        end
        subsets = [[i for i in 1:n if (mask >> (i - 1)) & 1 == 1] for mask in 0:(2^n-1)]
        for from in (true, false), closed in (false, true)
            S = unique(rand(rng, 1:n, rand(rng, 1:min(n, 2))))
            function valid(X)
                closed && !isempty(X) && !issubset(S, X) && return false
                h = closed ? TG._induced(g, X, ti) : g
                Rh = temporal_reachability(h, ti)
                id(v) = closed ? findfirst(==(v), X) : v
                return all(from ? Rh[id(s), id(x)] : Rh[id(x), id(s)] for s in S, x in X)
            end
            C = from ? source_component(g, S, ti; closed=closed) : sink_component(g, S, ti; closed=closed)
            @test valid(C) && all(X -> !valid(X) || issubset(X, C), subsets)
        end
        h = Graphs.SimpleGraph(n)
        for e in es
            (a <= e.t <= b && e.u != e.v) && Graphs.add_edge!(h, e.u, e.v)
        end
        @test sort(window_components(g, ti)) == sort(sort.(Graphs.connected_components(h)))
        # snapshots on the integer grid of ti
        L = b - a + 1
        snap = [Graphs.SimpleGraph(n) for _ in 1:L]
        for e in es
            (a <= e.t <= b && e.u != e.v) && Graphs.add_edge!(snap[e.t-a+1], e.u, e.v)
        end
        iscomp(X, k) = X in sort.(Graphs.connected_components(snap[k]))
        expected = Set{Tuple{Vector{Int},Tuple{Int,Int}}}()
        for i in 1:L, c in sort.(Graphs.connected_components(snap[i]))
            (length(c) < 2 || (i > 1 && iscomp(c, i - 1))) && continue
            j = i
            while j < L && iscomp(c, j + 1)
                j += 1
            end
            push!(expected, (c, (a + i - 1, a + j - 1)))
        end
        @test Set((p.nodes, p.interval) for p in persistent_components(g, ti; resolution=1)) == expected
        for T in 1:3
            H = [Graphs.SimpleGraph(n) for _ in 1:max(L-T+1, 0)]
            for (k, hw) in enumerate(H), e in Graphs.edges(snap[k])
                all(j -> Graphs.has_edge(snap[j], Graphs.src(e), Graphs.dst(e)), k:k+T-1) && Graphs.add_edge!(hw, Graphs.src(e), Graphs.dst(e))
            end
            okset(X) = T > L ? length(X) == 1 : all(hw -> Graphs.is_connected(Graphs.induced_subgraph(hw, X)[1]), H)
            good = [X for X in subsets if !isempty(X) && okset(X)]
            maxim = sort([X for X in good if !any(Y -> length(Y) > length(X) && issubset(X, Y), good)])
            @test sort(interval_connected_components(g, T, ti; resolution=1)) == maxim
        end
    end
    # a path 1 → 2 → 3 floods in 2 time units from 1; gossip is impossible (3 reaches nobody)
    g = OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1)])
    @test temporal_flooding_time(g, 1) == 2 && temporal_gossip_time(g) === nothing
    @test source_component(g, 1) == [1, 2, 3] && sink_component(g, 3) == [1, 2, 3] && sink_component(g, [2, 3]) == [1, 2]
    @test_throws ArgumentError edge_expansion(g; β=1, hubs=true)
    @test_throws ArgumentError persistent_components(OrderedEdgeList(2, [(1, 2, 0.5, 0.0)]))
    @test_throws ArgumentError interval_connected_components(g, 0)
end

@testset "Δ-temporal components, stream components and seed selection" begin
    rng = MersenneTwister(23)
    for _ in 1:60
        n = rand(rng, 1:6)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:8), rand(rng, 0:2)) for _ in 1:rand(rng, 0:14)])
        a, b = time_interval(g)
        subsets = [[i for i in 1:n if (mask >> (i - 1)) & 1 == 1] for mask in 1:(2^n-1)]
        # Δ-temporal components against all subsets and all windows on the integer grid
        for Δ in (0, 1, 3, 20), uni in (false, true)
            Rs = [temporal_reachability(g, (τ, min(τ + Δ, b))) for τ in a:(b - Δ < a ? a : b - Δ)]
            ok(X) = all(R -> all(uni ? (R[x, y] || R[y, x]) : (R[x, y] && R[y, x]) for x in X, y in X), Rs)
            good = filter(ok, subsets)
            maximal = sort([X for X in good if !any(Y -> length(Y) > length(X) && issubset(X, Y), good)]; by=c -> (-length(c), c))
            @test delta_temporal_connected_components(g, Δ; resolution=1, unilateral=uni) == maximal
        end
        @test delta_temporal_connected_components(g, b - a; resolution=1) == temporal_connected_components(g)
        # stream components against the components of the covered points of the
        # half-integer grid (doubled to integers), linked within a node if at most
        # max(gap, 1/2) apart and across an edge at its time stamp
        ti = rand(rng, Bool) ? (a, b) : (rand(rng, 0:3), rand(rng, 5:10))
        for gap in (0, 1, 2)
            es = [e for e in edges(g) if ti[1] <= e.t && e.t + e.tt <= ti[2]]
            points = unique!(sort!([(v, x) for e in es for v in (e.u, e.v) for x in 2e.t:2(e.t+e.tt)]))
            id = Dict(p => i for (i, p) in enumerate(points))
            h = Graphs.SimpleGraph(length(points))
            for (i, (v, x)) in enumerate(points), (j, (w, y)) in enumerate(points)
                v == w && 0 < y - x <= max(2gap, 1) && Graphs.add_edge!(h, i, j)
            end
            for e in es
                Graphs.add_edge!(h, id[(e.u, 2e.t)], id[(e.v, 2e.t)])
            end
            oracle = zeros(Int, length(points))
            for (c, comp) in enumerate(Graphs.connected_components(h))
                oracle[comp] .= c
            end
            C = stream_components(g, ti; gap=gap)
            mine = map(points) do (v, x)
                cs = [c for (c, comp) in enumerate(C) for (w, (s, f)) in comp if w == v && 2s <= x <= 2f]
                length(cs) == 1 ? cs[1] : 0
            end
            @test all(>(0), mine)
            @test length(unique(zip(mine, oracle))) == length(C) == length(unique(oracle))
            @test all(comp -> all(((v, (s, f)),) -> haskey(id, (v, 2s)) && haskey(id, (v, 2f)), comp), C)
        end
    end
    # seed selection against all subsets
    for _ in 1:60
        n = rand(rng, 1:8)
        g = OrderedEdgeList(n, [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:8), rand(rng, 0:2)) for _ in 1:rand(rng, 0:16)])
        R = temporal_reachability(g)
        reach(S) = count(v -> any(s -> R[s, v], S), 1:n)
        subsets = [[i for i in 1:n if (mask >> (i - 1)) & 1 == 1] for mask in 0:(2^n-1)]
        opt = minimum(length(S) for S in subsets if reach(S) == n)
        D = reachability_dominating_set(g)
        @test reach(D) == n && issorted(D) && allunique(D) && length(D) <= opt * (log(n) + 1)
        E = reachability_dominating_set(g; exact=true)
        @test reach(E) == n && length(E) == opt
        for k in 0:3
            best = maximum(reach(S) for S in subsets if length(S) <= k)
            r = max_reach_seeds(g, k)
            @test r.reached == reach(r.seeds) && length(r.seeds) <= k && r.reached >= (1 - 1 / ℯ) * best
            x = max_reach_seeds(g, k; exact=true)
            @test x.reached == best == reach(x.seeds) && length(x.seeds) <= k && issorted(x.seeds)
        end
    end
    g = random_simple_temporal_graph(60, 0.08; rng=rng)
    D, E = reachability_dominating_set(g), reachability_dominating_set(g; exact=true)
    R = temporal_reachability(g)
    @test length(E) <= length(D) && all(v -> any(s -> R[s, v], E), 1:60)
    @test max_reach_seeds(g, 3; exact=true).reached >= max_reach_seeds(g, 3).reached
    # a path 1 → 2 → 3 is dominated by 1; node 4 has no edges and must be a seed
    @test reachability_dominating_set(OrderedEdgeList(3, [(1, 2, 1, 1), (2, 3, 2, 1)])) == [1]
    @test max_reach_seeds(OrderedEdgeList(4, [(1, 2, 1, 1), (2, 3, 2, 1)]), 2) == (seeds=[1, 4], reached=4)
    @test_throws ArgumentError max_reach_seeds(OrderedEdgeList(2, [(1, 2, 1, 1)]), -1)
    @test_throws ArgumentError delta_temporal_connected_components(OrderedEdgeList(2, [(1, 2, 1, 1)]), -1)
    @test_throws ArgumentError stream_components(OrderedEdgeList(2, [(1, 2, 1, 1)]); gap=-1)
    # two meetings of 1 and 2 separated by more than the gap are two temporal nodes each
    s = OrderedEdgeList(3, [(1, 2, 0, 1), (2, 3, 1, 1), (1, 2, 5, 1)])
    @test stream_components(s) == [[(1, (0, 1)), (2, (0, 2)), (3, (1, 2))], [(1, (5, 6)), (2, (5, 6))]]
    @test length(stream_components(s; gap=3)) == 1
end
