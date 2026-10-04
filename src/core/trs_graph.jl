# Edge of the TRS graph: head time-node, its temporal graph node and the transition
# time, packed together so that a traversal reads one contiguous record per edge.
struct TRSArc{V<:Integer,T<:Real}
    tt::T
    head::V
    node::V
end

"""
    TRSGraph

Time-respecting static graph: a static DAG equivalent to a temporal graph with
respect to reachability, used for temporal distance computations.

Every TRS node is a time-node `(v, τ)`, i.e., a node `v` of the temporal graph at a
time `τ` at which `v` has an outgoing edge (plus the latest arrival time at `v`).
Time-nodes of the same node are consecutive and sorted by time, and they are
chained by waiting edges with transition time 0. A temporal edge `(u, v, t, tt)`
becomes an edge from `(u, t)` to the earliest time-node `(v, τ)` with `τ ≥ t + tt`.
Construct it with [`to_trs_graph`](@ref).
"""
struct TRSGraph{V<:Integer,T<:Real}
    num_tg_nodes::Int
    tg_node::Vector{V}        # temporal graph node of each TRS node
    time::Vector{T}           # time of each TRS node
    offsets::Vector{Int}      # CSR adjacency
    arcs::Vector{TRSArc{V,T}}
    first_pos::Vector{Int}    # first TRS node of each temporal graph node (0 if none)
    last_pos::Vector{Int}     # last TRS node of each temporal graph node (0 if none)
    ti::TimeInterval{T}
end

num_nodes(g::TRSGraph) = g.num_tg_nodes
num_edges(g::TRSGraph) = length(g.arcs)
time_interval(g::TRSGraph) = g.ti
time_type(::TRSGraph{V,T}) where {V,T} = T
node_type(::TRSGraph{V,T}) where {V,T} = V

"""
    num_trs_nodes(g::TRSGraph)

Number of time-nodes of the TRS graph.
"""
num_trs_nodes(g::TRSGraph) = length(g.time)

"""
    trs_node(g::TRSGraph, i)

The temporal graph node of time-node `i`.
"""
trs_node(g::TRSGraph, i::Integer) = Int(g.tg_node[i])

"""
    trs_time(g::TRSGraph, i)

The time of time-node `i`.
"""
trs_time(g::TRSGraph, i::Integer) = g.time[i]

"""
    trs_neighbors(g::TRSGraph, i)

The out-neighbors of time-node `i` and the corresponding transition times, as a
tuple of two vectors. The waiting edge, if any, comes first.
"""
function trs_neighbors(g::TRSGraph, i::Integer)
    arcs = view(g.arcs, g.offsets[i]:g.offsets[i+1]-1)
    return [Int(a.head) for a in arcs], [a.tt for a in arcs]
end

# index of the first time-node of v with time ≥ τ
function _trs_position(times, first_pos, last_pos, v, τ)
    lo, hi = first_pos[v], last_pos[v]
    return lo - 1 + searchsortedfirst(view(times, lo:hi), τ)
end

"""
    to_trs_graph(g::OrderedEdgeList)

Transform `g` into its time-respecting static graph representation [`TRSGraph`](@ref).
"""
function to_trs_graph(g::OrderedEdgeList{V,T}) where {V,T}
    n = num_nodes(g)
    es = edges(g)
    max_arrival = fill(typemin(T), n)
    time_nodes = Vector{Tuple{V,T}}()
    sizehint!(time_nodes, length(es) + n)
    for e in es
        push!(time_nodes, (e.u, e.t))
        max_arrival[e.v] = max(max_arrival[e.v], e.t + e.tt)
    end
    for v in 1:n
        max_arrival[v] > typemin(T) && push!(time_nodes, (V(v), max_arrival[v]))
    end
    sort!(time_nodes)
    unique!(time_nodes)
    N = length(time_nodes)
    N <= typemax(V) || throw(OverflowError("the TRS graph has $N time-nodes, use a larger node id type than $V"))
    tg_node = first.(time_nodes)
    times = last.(time_nodes)

    first_pos = zeros(Int, n)
    last_pos = zeros(Int, n)
    for i in N:-1:1
        first_pos[tg_node[i]] = i
    end
    for i in 1:N
        last_pos[tg_node[i]] = i
    end

    tails = Vector{Int}(undef, length(es))
    arrivals = Vector{Int}(undef, length(es))
    degree = zeros(Int, N)
    for i in 1:N
        (i < N && tg_node[i+1] == tg_node[i]) && (degree[i] += 1)  # waiting edge
    end
    for (k, e) in enumerate(es)
        tails[k] = _trs_position(times, first_pos, last_pos, e.u, e.t)
        arrivals[k] = _trs_position(times, first_pos, last_pos, e.v, e.t + e.tt)
        degree[tails[k]] += 1
    end
    offsets = Vector{Int}(undef, N + 1)
    offsets[1] = 1
    for i in 1:N
        offsets[i+1] = offsets[i] + degree[i]
    end
    arcs = Vector{TRSArc{V,T}}(undef, offsets[end] - 1)
    pos = offsets[1:N]
    for i in 1:N
        if i < N && tg_node[i+1] == tg_node[i]
            arcs[pos[i]] = TRSArc{V,T}(zero(T), i + 1, tg_node[i])
            pos[i] += 1
        end
    end
    for (k, e) in enumerate(es)
        p = pos[tails[k]]
        arcs[p] = TRSArc{V,T}(e.tt, arrivals[k], e.v)
        pos[tails[k]] += 1
    end
    return TRSGraph{V,T}(n, tg_node, times, offsets, arcs, first_pos, last_pos, time_interval(g))
end

function Base.show(io::IO, g::TRSGraph{V,T}) where {V,T}
    print(io, "TRSGraph{", V, ", ", T, "} with ", num_trs_nodes(g), " time-nodes and ", num_edges(g),
          " edges for a temporal graph with ", g.num_tg_nodes, " nodes")
end
