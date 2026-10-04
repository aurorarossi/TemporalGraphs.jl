# Randomized reference models.

@testset "randomized reference models (Gauvin et al.)" begin
    key(e, dir) = dir ? (Int(e.u), Int(e.v)) : minmax(Int(e.u), Int(e.v))
    links(g; dir=true) = Set(key(e, dir) for e in edges(g))
    function lines(g; dir=true)
        d = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
        for e in edges(g)
            push!(get!(d, key(e, dir), Tuple{Int,Int}[]), (e.t, e.tt))
        end
        foreach(sort!, values(d))
        return d
    end
    pool(g; dir=true) = sort!(collect(values(lines(g; dir=dir))))
    weights(g) = Dict(k => length(v) for (k, v) in lines(g))
    times(g) = sort!([(e.t, e.tt) for e in edges(g)])
    function degrees(es, n)
        o, i = zeros(Int, n), zeros(Int, n)
        for (u, v) in Set((Int(e.u), Int(e.v)) for e in es)
            o[u] += 1
            i[v] += 1
        end
        return o, i
    end
    function snaps(g)
        d = Dict{Int,Vector{Tuple{Int,Int}}}()
        for e in edges(g)
            push!(get!(d, e.t, Tuple{Int,Int}[]), (e.u, e.v))
        end
        foreach(sort!, values(d))
        return d
    end
    function instdegrees(g)
        d = Dict{Tuple{Int,Int,Int},Int}()
        for e in edges(g)
            d[(e.t, e.u, 1)] = get(d, (e.t, e.u, 1), 0) + 1
            d[(e.t, e.v, 2)] = get(d, (e.t, e.v, 2), 0) + 1
        end
        return d
    end
    gaps(g) = Dict(k => (v[1][1], sort(diff(first.(v)))) for (k, v) in lines(g))
    rng = MersenneTwister(11)
    for _ in 1:60
        n = rand(rng, 3:7)
        es = unique([TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:20), rand(rng, 0:2)) for _ in 1:rand(rng, 1:15)])
        es = [e for e in es if e.u != e.v]
        isempty(es) && continue
        g = OrderedEdgeList(n, es)
        a, b = time_interval(g)
        h = randomize(g, LinkShuffling(); rng=rng)
        @test length(links(h)) == length(links(g)) && pool(h) == pool(g)
        h = randomize(g, DegreeLinkShuffling(); rng=rng)
        @test degrees(edges(h), n) == degrees(edges(g), n) && pool(h) == pool(g)
        h = randomize(g, TopologyLinkShuffling(); rng=rng)
        @test links(h) == links(g) && pool(h) == pool(g)
        h = randomize(g, WeightLinkShuffling(); rng=rng)
        @test weights(h) == weights(g) && pool(h) == pool(g)
        h = randomize(g, RandomTimes(); rng=rng)
        @test weights(h) == weights(g) && all(a <= e.t && e.t + e.tt <= b for e in edges(h))
        @test all(allunique, values(lines(h)))
        h = randomize(g, InterEventShuffling(); rng=rng)
        @test gaps(h) == gaps(g)
        h = randomize(g, TimelineShifting(); rng=rng)
        @test weights(h) == weights(g) && all(a <= e.t <= b for e in edges(h))
        h = randomize(g, TimestampShuffling(); rng=rng)
        @test weights(h) == weights(g) && times(h) == times(g)
        h = randomize(g, TimestampShuffling(simple=true); rng=rng)
        @test weights(h) == weights(g) && times(h) == times(g) && allunique(edges(h))
        for active in (true, false)
            h = randomize(g, SequenceShuffling(active=active); rng=rng)
            @test sort!(collect(values(snaps(h)))) == sort!(collect(values(snaps(g))))
            active && @test Set(keys(snaps(h))) == Set(keys(snaps(g)))
        end
        h = randomize(g, IsomorphicSnapshotShuffling(); rng=rng)
        @test Dict(t => length(s) for (t, s) in snaps(h)) == Dict(t => length(s) for (t, s) in snaps(g))
        h = randomize(g, DegreeSnapshotShuffling(); rng=rng)
        @test instdegrees(h) == instdegrees(g)
        h = randomize(g, EventShuffling(); rng=rng)
        @test allunique(edges(h)) && sort([e.tt for e in edges(h)]) == sort([e.tt for e in edges(g)])
        h = randomize(g, TopologyLinkShuffling(directed=false), DegreeLinkShuffling(directed=false); rng=rng)
        @test pool(h; dir=false) == pool(g; dir=false)
        for m in (LinkShuffling(), SnapshotShuffling(), EventShuffling())  # need enough pairs of nodes
            h = try
                randomize(g, m; rng=rng)
            catch err
                err isa ArgumentError || rethrow()
                nothing
            end
            h === nothing && continue
            @test num_edges(h) == num_edges(g) && time_interval(h) == time_interval(g) && all(e.u != e.v for e in edges(h))
            m isa SnapshotShuffling && @test times(h) == times(g) && all(allunique, values(snaps(h)))
        end
    end
    # uniformity on small outcome spaces: χ² within 4 standard deviations of its mean
    function χ²(g, m, N)
        count = Dict{Vector{Tuple{Int,Int,Int,Int}},Int}()
        r = MersenneTwister(7)
        for _ in 1:N
            k = sort([(Int(e.u), Int(e.v), Int(e.t), Int(e.tt)) for e in edges(randomize(g, m; rng=r))])
            count[k] = get(count, k, 0) + 1
        end
        e = N / length(count)
        return length(count), sum((x - e)^2 / e for x in values(count))
    end
    g1 = OrderedEdgeList(3, [(1, 2, 0, 0), (2, 3, 0, 0), (3, 1, 1, 0)])
    g2 = OrderedEdgeList(3, [(1, 2, 0, 0), (1, 2, 2, 0), (2, 3, 1, 0)])
    g3 = OrderedEdgeList(4, [(1, 2, 0, 1), (2, 3, 1, 1), (1, 2, 3, 1), (3, 4, 3, 1), (4, 1, 4, 1)])
    for (g, m, outcomes) in ((g1, SnapshotShuffling(), 90), (g1, EventShuffling(), 220), (g1, IsomorphicSnapshotShuffling(), 36),
                             (g1, LinkShuffling(), 60), (g2, RandomTimes(), 9), (g3, TimestampShuffling(simple=true), 27),
                             (g3, DegreeLinkShuffling(), 216), (g3, TimelineShifting(), 648), (g3, SequenceShuffling(active=false), 360))
        k, x = χ²(g, m, 40000)
        @test k == outcomes
        @test x < k - 1 + 4 * sqrt(2 * (k - 1))
    end
    # parallel samples are reproducible
    g = random_graph(MersenneTwister(2), 20, 300; ti=(0, 300))
    s1 = reference_samples(length ∘ unique ∘ edges, g, TimestampShuffling(); samples=8, rng=MersenneTwister(1))
    s2 = reference_samples(length ∘ unique ∘ edges, g, TimestampShuffling(); samples=8, rng=MersenneTwister(1))
    @test s1 == s2 && length(s1) == 8
    @test reference_samples(h -> sum(temporal_motif_counts(h, 10)), g, SnapshotShuffling(); samples=3) isa Vector{Int}
    @test randomize(g) == g
    @test_throws ArgumentError randomize(OrderedEdgeList(2, [(1, 2, 0, 0), (2, 1, 0, 0), (1, 2, 1, 0)]), SnapshotShuffling(directed=false))
    @test_throws ArgumentError randomize(OrderedEdgeList(2, [(1, 2, 0.5, 0.0)]), SequenceShuffling(active=false))
    h = randomize(OrderedEdgeList(3, [(1, 2, 0, 0), (1, 2, 0, 0), (2, 3, 1, 0), (1, 3, 2, 0)]), TimestampShuffling(simple=true))
    @test count(e -> (e.u, e.v, e.t) == (1, 2, 0), edges(h)) <= 2 && sort([e.t for e in edges(h)]) == [0, 0, 1, 2]
end
