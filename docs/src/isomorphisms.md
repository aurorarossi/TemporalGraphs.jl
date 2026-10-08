# Temporal graph isomorphisms

When are two temporal graphs "the same"? Heeg, Sauer, Mutzel and Scholtes (2025)
compare three notions, from the weakest to the strongest:

| `kind` | The graphs are isomorphic if… |
|:---|:---|
| `:aggregated` | the static graphs with the number of temporal edges of every pair are isomorphic (*time-aggregated isomorphism*) |
| `:event` | a bijection of the nodes and of the temporal edges preserves all time-respecting paths (*time-respecting path isomorphism*) |
| `:concatenated` | the static graphs with the time stamps of every pair (relative to the first one) are isomorphic (*time-concatenated* or *timewise isomorphism*, Wałęga and Rawson, 2025) |

## Example

The `:event` notion only cares about the causal structure: time stamps may change, as
long as the same time-respecting paths exist. In the example of Figure 2 of the paper,
`G3` swaps the times of two edges that lie on no common time-respecting path. The
paper uses a maximum time difference ``δ = 2`` between consecutive edges with unit
transition times, which is `β = δ - 1 = 1` here.

The graph `G1`:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(5, [(1, 2, 1, 1), (2, 4, 2, 1), (3, 4, 3, 1), (4, 5, 4, 1)])) # hide
```

`G3` swaps the times of `(2, 4)` and `(3, 4)`. No time-respecting path uses both, so
`G3` has the same paths as `G1`:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(5, [(1, 2, 1, 1), (2, 4, 3, 1), (3, 4, 2, 1), (4, 5, 4, 1)])) # hide
```

`G4` moves `(4, 5)` to time 1, which destroys the paths `1 → 2 → 4 → 5` and
`3 → 4 → 5` but keeps the aggregated graph:

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(5, [(1, 2, 1, 1), (2, 4, 2, 1), (3, 4, 3, 1), (4, 5, 1, 1)])) # hide
```

```jldoctest iso
julia> using TemporalGraphs

julia> G(es) = OrderedEdgeList(5, [(u, v, t, 1) for (u, v, t) in es]);

julia> G1 = G([(1, 2, 1), (2, 4, 2), (3, 4, 3), (4, 5, 4)]);

julia> G3 = G([(1, 2, 1), (2, 4, 3), (3, 4, 2), (4, 5, 4)]);

julia> G4 = G([(1, 2, 1), (2, 4, 2), (3, 4, 3), (4, 5, 1)]);

julia> kinds = (:aggregated, :event, :concatenated);

julia> [is_temporally_isomorphic(G1, H; kind = k, β = 1) for k in kinds, H in (G3, G4)]
3×2 Matrix{Bool}:
 1  1
 1  0
 0  0

julia> temporal_isomorphism(G1, G3; β = 1)   # nodes and edges of G1 → G3
(nodes = [1, 2, 3, 4, 5], edges = [1, 3, 2, 4])
```

## How it is computed

Time-respecting path isomorphism is equivalent to static
isomorphism of the [`augmented_event_graph`](@ref)s (Theorem 1 of the paper). This
graph has a node for every temporal edge, an arc between consecutive edges of
time-respecting paths, and the original nodes connected to their edges.

[`temporal_isomorphism`](@ref) decides isomorphism exactly, by individualization and
refinement:

- the two graphs are colored together with the directed Weisfeiler–Leman algorithm
  (D-WL);
- when the coloring stops changing, one node per graph is given a fresh color, and
  the coloring is refined incrementally;
- the search backtracks through an undo trail if the color classes of the two graphs
  differ.

Graph isomorphism has no known polynomial algorithm, but this is fast on temporal
networks. On the whole `sp-workplace-2013` dataset of [`load_dataset`](@ref), with an
hour as maximum time difference, an isomorphism to a copy with renamed nodes and
event graph components rearranged in time is found in less than a second.

## Weisfeiler–Leman test and kernel

- [`temporal_wl_equivalent`](@ref) is the one-sided test of the paper. If D-WL
  distinguishes the augmented event graphs, the temporal graphs are not isomorphic;
  otherwise neither D-WL nor message passing neural networks on the augmented event
  graph can tell them apart.
- [`temporal_wl_kernel`](@ref) computes the Weisfeiler–Leman subtree kernel of a
  collection of temporal graphs, for temporal graph classification. With
  `augmented = false` it runs on the event graph alone, as the temporal graph kernel
  of Oettershagen, Kriege, Morris and Mutzel (2020). That variant does not see which
  nodes the temporal edges join.

## References

- F. Heeg, J. Sauer, P. Mutzel, and I. Scholtes. *Weisfeiler and Leman follow the arrow of time: expressive power of message passing in temporal event graphs.* 2025. [arXiv](https://arxiv.org/abs/2505.24438)
- P. A. Wałęga and M. Rawson. *Expressive power of temporal message passing.* AAAI, 2025. [DOI](https://doi.org/10.1609/aaai.v39i20.35396)
- L. Oettershagen, N. M. Kriege, C. Morris, and P. Mutzel. *Temporal graph kernels for classifying dissemination processes.* SDM, 2020. [DOI](https://doi.org/10.1137/1.9781611976236.56), [arXiv](https://arxiv.org/abs/1911.05496)
