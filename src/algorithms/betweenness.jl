# Brandes' accumulation for the sources in `range`, added into `bc`.
function _brandes!(bc::Vector{Float64}, dlg::DirectedLineGraph, range)
    N = num_nodes(dlg)
    d = fill(-1, N)
    sigma = zeros(Float64, N)
    delta = zeros(Float64, N)
    order = Int[]
    @inbounds for s in range
        empty!(order)
        push!(order, s)
        d[s] = 0
        sigma[s] = 1.0
        i = 1
        while i <= length(order)
            u = order[i]
            i += 1
            for v in out_edges(dlg, u)
                if d[v] < 0
                    d[v] = d[u] + 1
                    push!(order, v)
                end
                d[v] == d[u] + 1 && (sigma[v] += sigma[u])
            end
        end
        for j in length(order):-1:1
            w = order[j]
            for v in out_edges(dlg, w)
                d[v] == d[w] + 1 && (delta[w] += sigma[w] / sigma[v] * (1 + delta[v]))
            end
            w != s && (bc[w] += delta[w])
        end
        for v in order  # reset only what this source touched
            d[v] = -1
            sigma[v] = 0.0
            delta[v] = 0.0
        end
    end
    return bc
end

"""
    temporal_edge_betweenness(dlg::DirectedLineGraph)
    temporal_edge_betweenness(g::OrderedEdgeList)

Temporal edge betweenness: Brandes' betweenness of the nodes of the directed line
graph, i.e., for every temporal edge the sum over all pairs of edges `(e, f)` of the
fraction of shortest (minimum hop) temporal paths from `e` to `f` that use it.
Entry `i` belongs to `edges(dlg)[i]`, respectively `edges(g)[i]`. The sources are
processed in parallel chunks with one accumulator per chunk.
"""
function temporal_edge_betweenness(dlg::DirectedLineGraph)
    N = num_nodes(dlg)
    N == 0 && return Float64[]
    nchunks = min(N, 4 * Threads.nthreads())
    return tmapreduce(+, chunks(1:N; n=nchunks)) do range
        _brandes!(zeros(Float64, N), dlg, range)
    end
end

temporal_edge_betweenness(g::OrderedEdgeList) = temporal_edge_betweenness(to_directed_line_graph(g))
