# Maximum flow on a static directed network with Dinic's algorithm, used on the
# time-expanded networks of temporal graphs. Arcs are stored in pairs (arc, reverse
# arc) in CSR order after `_build!`.

mutable struct _FlowNetwork{C<:Real}
    n::Int
    tail::Vector{Int}
    head::Vector{Int}
    cap::Vector{C}       # residual capacities
    off::Vector{Int}     # CSR of the arcs leaving each node
    arcs::Vector{Int}
    level::Vector{Int}
    iter::Vector{Int}
end

_FlowNetwork{C}(n::Integer) where {C} = _FlowNetwork{C}(Int(n), Int[], Int[], C[], Int[], Int[], Int[], Int[])

# Add the arc u → v with capacity c (and its reverse arc); returns the arc index.
function _add_arc!(net::_FlowNetwork{C}, u::Integer, v::Integer, c) where {C}
    push!(net.tail, u, v)
    push!(net.head, v, u)
    push!(net.cap, convert(C, c), zero(C))
    return length(net.head) - 1
end

_reverse_arc(a::Int) = isodd(a) ? a + 1 : a - 1

function _build!(net::_FlowNetwork)
    n = net.n
    off = zeros(Int, n + 1)
    for u in net.tail
        off[u+1] += 1
    end
    off[1] = 1
    for u in 1:n
        off[u+1] += off[u]
    end
    pos = off[1:n]
    arcs = Vector{Int}(undef, length(net.tail))
    for (a, u) in enumerate(net.tail)
        arcs[pos[u]] = a
        pos[u] += 1
    end
    net.off, net.arcs = off, arcs
    net.level = zeros(Int, n)
    net.iter = zeros(Int, n)
    return net
end

function _bfs_levels!(net::_FlowNetwork, s::Int, z::Int, queue::Vector{Int})
    fill!(net.level, -1)
    net.level[s] = 0
    empty!(queue)
    push!(queue, s)
    h = 1
    @inbounds while h <= length(queue)
        u = queue[h]
        h += 1
        for x in net.off[u]:net.off[u+1]-1
            a = net.arcs[x]
            v = net.head[a]
            if net.cap[a] > 0 && net.level[v] < 0
                net.level[v] = net.level[u] + 1
                push!(queue, v)
            end
        end
    end
    return net.level[z] >= 0
end

# Blocking flow by iterative depth-first search with current-arc pointers.
function _blocking_flow!(net::_FlowNetwork{C}, s::Int, z::Int, stack::Vector{Int}) where {C}
    total = zero(C)
    @inbounds for u in 1:net.n
        net.iter[u] = net.off[u]
    end
    @inbounds while true
        empty!(stack)  # arcs of the current path
        u = s
        while u != z
            advanced = false
            while net.iter[u] < net.off[u+1]
                a = net.arcs[net.iter[u]]
                v = net.head[a]
                if net.cap[a] > 0 && net.level[v] == net.level[u] + 1
                    push!(stack, a)
                    u = v
                    advanced = true
                    break
                end
                net.iter[u] += 1
            end
            if !advanced
                u == s && return total
                net.level[u] = -1  # dead end
                a = pop!(stack)
                u = net.tail[a]
                net.iter[u] += 1
            end
        end
        δ = minimum(net.cap[a] for a in stack)
        for a in stack
            net.cap[a] -= δ
            net.cap[_reverse_arc(a)] += δ
        end
        total += δ
    end
end

# Maximum flow value from s to z (the residual capacities are left in net).
function _max_flow!(net::_FlowNetwork{C}, s::Int, z::Int) where {C}
    s == z && return zero(C)
    queue = Int[]
    stack = Int[]
    flow = zero(C)
    while _bfs_levels!(net, s, z, queue)
        flow += _blocking_flow!(net, s, z, stack)
    end
    return flow
end

# Nodes reachable from s in the residual network (the source side of a minimum cut).
function _residual_reachable(net::_FlowNetwork, s::Int)
    seen = falses(net.n)
    seen[s] = true
    queue = [s]
    h = 1
    while h <= length(queue)
        u = queue[h]
        h += 1
        for x in net.off[u]:net.off[u+1]-1
            a = net.arcs[x]
            v = net.head[a]
            if net.cap[a] > 0 && !seen[v]
                seen[v] = true
                push!(queue, v)
            end
        end
    end
    return seen
end
