# Test instances of TGLib (test/TemporalGraphsInstances.cpp) with 1-based node ids.

tg(n, es, ti=nothing) = ti === nothing ?
    OrderedEdgeList(n, [TemporalEdge(u + 1, v + 1, t, tt) for (u, v, t, tt) in es]) :
    OrderedEdgeList(n, [TemporalEdge(u + 1, v + 1, t, tt) for (u, v, t, tt) in es], ti)

example_from_paper() = tg(4, [(0, 3, 1, 5), (0, 1, 2, 1), (0, 1, 5, 2), (2, 1, 6, 1),
                              (3, 2, 6, 2), (1, 3, 7, 2), (3, 2, 8, 4)], (1, 12))

example_from_paper_with_loops_and_multi_edges() =
    tg(4, [(0, 3, 1, 5), (0, 3, 1, 5), (0, 1, 2, 1), (0, 1, 5, 2), (0, 1, 5, 2), (0, 0, 5, 2),
           (2, 1, 6, 1), (3, 2, 6, 2), (1, 3, 7, 2), (1, 3, 7, 2), (3, 2, 8, 4), (3, 3, 8, 4)], (1, 12))

example3() = tg(4, [(0, 1, 1, 1), (1, 2, 2, 1), (2, 3, 3, 1), (1, 3, 4, 1)], (0, 5))

chain() = tg(8, [(i, i + 1, i + 1, 1) for i in 0:6])

betweenness_test() = tg(7, [(0, 1, 1, 1), (1, 2, 2, 1), (2, 3, 3, 1), (3, 4, 4, 1), (3, 5, 4, 1), (5, 6, 5, 1)])

complete_graph() = tg(5, [(u, v, 1, 1) for u in 0:4 for v in 0:4 if u != v])

lkcore_tg() = tg(5, [(0, 1, 1, 1), (0, 2, 1, 1), (1, 2, 1, 1), (0, 1, 2, 1), (0, 2, 2, 1),
                     (1, 2, 2, 1), (0, 3, 2, 1), (1, 4, 3, 1), (2, 4, 4, 1)])

non_bursty(m) = tg(2, vcat([[(0, 1, i, 1), (1, 0, i, 1)] for i in 0:m-1]...))

function bursty(rng, m)
    es = NTuple{4,Int}[]
    for _ in 1:m
        t = rand(rng, 0:999) + (rand(rng, Bool) ? 1_000_000 : 0)
        push!(es, (0, 1, t, 1), (1, 0, t, 1))
    end
    return tg(2, es)
end

no_topological_overlap(m) = tg(3, vcat([[(0, 1, i, 1), (0, 2, i + 1, 1)] for i in 0:2:m-1]...))
full_topological_overlap(m) = tg(2, [(0, 1, i, 1) for i in 0:m-1])

# Random graph without loops and with unit transition times (getRandomGraphExample).
function random_graph(rng, n, m; ti=(0, 1_000_000_000))
    es = Set{TemporalEdge{Int32,Int64}}()
    for _ in 1:m
        u, v = rand(rng, 1:n), rand(rng, 1:n)
        u == v && continue
        push!(es, TemporalEdge(u, v, rand(rng, 1:m), 1))
    end
    return OrderedEdgeList(n, sort!(collect(es)), ti)
end

# Random graph with loops, multi-edges and arbitrary transition times (≥ min_tt).
function random_general_graph(rng, n, m; tmax=20, ttmax=4, min_tt=0)
    es = [TemporalEdge(rand(rng, 1:n), rand(rng, 1:n), rand(rng, 0:tmax), rand(rng, min_tt:ttmax)) for _ in 1:m]
    return OrderedEdgeList(n, es)
end
