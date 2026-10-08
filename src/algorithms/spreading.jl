# Seed selection by temporal reachability: the smallest sets of nodes that reach all
# nodes (temporal reachability dominating sets) and the k nodes that reach the most
# nodes (influence maximization when every contact transmits). Both are covering
# problems on the rows of the reachability matrix, stored as bitsets.

# Rows of the reachability matrix as bitsets: cover[s] holds the nodes reached from s.
function _reach_sets(g::OrderedEdgeList, ti)
    R = temporal_reachability(g, ti)
    n = size(R, 1)
    cover = [zeros(UInt64, cld(n, 64)) for _ in 1:n]
    for v in 1:n, s in 1:n
        R[s, v] && _setbit!(cover[s], v)
    end
    return cover
end

# Number of nodes of c that are not in covered.
function _gain(c::Vector{UInt64}, covered::Vector{UInt64})
    k = 0
    @inbounds for w in eachindex(c)
        k += count_ones(c[w] & ~covered[w])
    end
    return k
end

# Greedy covering: repeatedly the node that reaches the most nodes not reached yet,
# until k seeds are chosen or all nodes are reached. Gains only decrease, so they are
# evaluated lazily: a node is chosen when its gain, computed in the current round, is
# the largest stored one. Returns the seeds and the number of nodes reached.
function _greedy_cover(cover::Vector{Vector{UInt64}}, k::Int)
    n = length(cover)
    covered = zeros(UInt64, cld(n, 64))
    heap = BinaryMinHeap{Tuple{Int,Int,Int}}()  # (-gain, node, round of the gain)
    for s in 1:n
        push!(heap, (-_gain(cover[s], covered), s, 0))
    end
    seeds = Int[]
    reached = 0
    while length(seeds) < k && reached < n
        ng, s, round = pop!(heap)
        if round == length(seeds)
            push!(seeds, s)
            reached -= ng
            covered .|= cover[s]
        else
            push!(heap, (-_gain(cover[s], covered), s, length(seeds)))
        end
    end
    return seeds, reached
end

# Smallest set of nodes reaching all nodes, by branch and bound starting from the
# solution `best`: some source of the unreached node with the fewest allowed sources is
# a seed, and after its branch a source is excluded from the next branches. Pruned with
# |seeds| + ⌈unreached / largest gain⌉.
function _exact_dominating(cover::Vector{Vector{UInt64}}, best::Vector{Int})
    n = length(cover)
    sources = [Int[] for _ in 1:n]  # the nodes reaching every node
    for s in 1:n, v in 1:n
        _bit(cover[s], v) && push!(sources[v], s)
    end
    allowed = trues(n)
    chosen = Int[]
    function search(covered::Vector{UInt64}, reached::Int)
        if reached == n
            length(chosen) < length(best) && (best = copy(chosen))
            return
        end
        gains = [allowed[s] ? _gain(cover[s], covered) : 0 for s in 1:n]
        gmax = maximum(gains)
        (gmax == 0 || length(chosen) + cld(n - reached, gmax) >= length(best)) && return
        v, few = 0, typemax(Int)
        for x in 1:n
            _bit(covered, x) && continue
            c = count(s -> allowed[s], sources[x])
            c < few && ((v, few) = (x, c))
        end
        few == 0 && return
        excluded = Int[]
        for s in sort!(filter(s -> allowed[s], sources[v]); by=s -> -gains[s])
            push!(chosen, s)
            search(covered .| cover[s], reached + gains[s])
            pop!(chosen)
            allowed[s] = false
            push!(excluded, s)
        end
        allowed[excluded] .= true
    end
    search(zeros(UInt64, cld(n, 64)), 0)
    return best
end

# At most k nodes reaching the most nodes, by branch and bound from the solution
# (best, bestval): every node, by decreasing number of reached nodes, is a seed or not;
# pruned with the sum of the largest gains of the remaining nodes.
function _exact_max_reach(cover::Vector{Vector{UInt64}}, k::Int, best::Vector{Int}, bestval::Int)
    n = length(cover)
    empty = zeros(UInt64, cld(n, 64))
    order = sortperm([_gain(c, empty) for c in cover]; rev=true)
    chosen = Int[]
    function search(i::Int, covered::Vector{UInt64}, reached::Int)
        if reached > bestval
            best, bestval = copy(chosen), reached
        end
        (length(chosen) == k || i > n || reached == n) && return
        gains = sort!([_gain(cover[order[j]], covered) for j in i:n]; rev=true)
        reached + sum(view(gains, 1:min(k - length(chosen), length(gains)))) <= bestval && return
        s = order[i]
        push!(chosen, s)
        search(i + 1, covered .| cover[s], reached + _gain(cover[s], covered))
        pop!(chosen)
        search(i + 1, covered, reached)
    end
    search(1, empty, 0)
    return best, bestval
end

"""
    reachability_dominating_set(g::OrderedEdgeList, ti = time_interval(g); exact = false)

A small *temporal reachability dominating set* (TaRDiS) of `g` (Kutner and
Larios-Jones, 2026): a set `S` of nodes such that every node is in `S` or is reached
from a node of `S` by a temporal path inside `ti`. If the nodes of `S` are infected at
the start of `ti` and every contact transmits, every node gets infected.

Finding a smallest one is NP-hard (it is a set cover problem on the rows of the
[`temporal_reachability`](@ref) matrix). By default the greedy algorithm repeatedly
chooses the node that reaches the most nodes not reached yet, which is at most
`ln n + 1` times larger than optimal (Chvátal, 1979). With `exact = true` a smallest
set is found by branch and bound from the greedy solution (exponential in the worst
case). Returns the nodes of the set, sorted.

# References

- D. C. Kutner and L. Larios-Jones. *Temporal reachability dominating sets: contagion in temporal graphs.* Journal of Computer and System Sciences 155, 2026. [DOI](https://doi.org/10.1016/j.jcss.2025.103701), [arXiv](https://arxiv.org/abs/2306.06999)
- V. Chvátal. *A greedy heuristic for the set-covering problem.* Mathematics of Operations Research 4(3), 1979. [DOI](https://doi.org/10.1287/moor.4.3.233)
"""
function reachability_dominating_set(g::OrderedEdgeList, ti=time_interval(g); exact::Bool=false)
    cover = _reach_sets(g, ti)
    seeds, _ = _greedy_cover(cover, num_nodes(g))
    exact && (seeds = _exact_dominating(cover, seeds))
    return sort!(seeds)
end

"""
    max_reach_seeds(g::OrderedEdgeList, k, ti = time_interval(g); exact = false)

At most `k` nodes that together reach the most nodes of `g` by temporal paths inside
`ti`: influence maximization (Kempe, Kleinberg and Tardos, 2003) in the deterministic
model where the seeds are informed at the start of `ti` and every contact transmits.
Returns `(seeds, reached)`: the seeds, sorted, and the number of nodes they reach
(including themselves). Fewer than `k` seeds are returned only if they reach all
nodes.

The problem is NP-hard (maximum coverage on the rows of the
[`temporal_reachability`](@ref) matrix). By default the greedy algorithm, with lazy
evaluation of the gains, reaches at least `1 - 1/e ≈ 63%` of the optimum (Nemhauser,
Wolsey and Fisher, 1978); with `exact = true` an optimal set is found by branch and
bound (exponential in the worst case).

# References

- D. Kempe, J. Kleinberg, and É. Tardos. *Maximizing the spread of influence through a social network.* KDD, 2003. [DOI](https://doi.org/10.1145/956750.956769)
- G. L. Nemhauser, L. A. Wolsey, and M. L. Fisher. *An analysis of approximations for maximizing submodular set functions—I.* Mathematical Programming 14, 1978. [DOI](https://doi.org/10.1007/BF01588971)
"""
function max_reach_seeds(g::OrderedEdgeList, k::Integer, ti=time_interval(g); exact::Bool=false)
    k >= 0 || throw(ArgumentError("k must be non-negative"))
    cover = _reach_sets(g, ti)
    seeds, reached = _greedy_cover(cover, Int(k))
    exact && ((seeds, reached) = _exact_max_reach(cover, Int(k), seeds, reached))
    return (seeds=sort!(seeds), reached=reached)
end
