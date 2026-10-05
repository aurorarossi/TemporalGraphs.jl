# Motifs and reference models

## Temporal motifs

A *δ-temporal motif* (Paranjape, Benson and Leskovec, 2017) is a sequence of ``l``
directed edges on ``k`` nodes, such as "``a`` writes to ``b``, then ``b`` answers,
then ``a`` writes to ``c``". An *instance* is a sequence of ``l`` temporal edges of the
graph that

- maps onto the motif by a bijection of the nodes,
- is ordered in time, and
- spans at most ``δ`` time units, ``t_l - t_1 ≤ δ``.

The edges of an instance need not be consecutive, and other edges may occur in
between. Counting instances at several time scales ``δ`` shows how communication
patterns form: for example, the difference between the counts for ``δ = 1`` hour
and ``δ = 30`` minutes counts the patterns that take between 30 minutes and an hour
to form.

[`temporal_motif_counts`](@ref) counts the 36 motifs with 3 edges on 2 or 3 nodes
and returns a ``6 × 6`` matrix. The layout follows Figure 3 of the paper and the SNAP
implementation. Call ``a → b`` the first edge and ``c`` the third node:

- column `j` is given by the third edge: ``a → b, b → a, a → c, c → a, b → c, c → b``;
- row `i` is given by the second edge, in the reverse order.

In this layout:

- `M[5:6, 1:2]` are the four 2-node motifs;
- `M[1:2, 3:4]` and `M[3:4, 5:6]` are the eight triangles;
- the other 24 entries are stars.

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(3, [(1, 2, 10, 0), (2, 1, 20, 0), (1, 3, 25, 0), (3, 2, 40, 0), (1, 2, 70, 0)])) # hide
```

```jldoctest motifs
julia> using TemporalGraphs

julia> g = OrderedEdgeList(3, [(1, 2, 10, 0), (2, 1, 20, 0), (1, 3, 25, 0), (3, 2, 40, 0), (1, 2, 70, 0)]);

julia> temporal_motif_counts(g, 30)
6×6 Matrix{Int64}:
 0  0  0  0  0  0
 0  0  0  1  0  0
 0  0  0  0  0  0
 0  0  0  0  0  1
 0  0  1  0  0  1
 0  0  0  0  0  0
```

The four instances are:

| Instance | Cell | Motif |
|---|---|---|
| ``1→2, 2→1, 1→3`` | `M[5, 3]` | a star |
| ``1→2, 2→1, 3→2`` | `M[5, 6]` | a star |
| ``1→2, 1→3, 3→2`` | `M[4, 6]` | a triangle |
| ``2→1, 1→3, 3→2`` | `M[2, 4]` | a cyclic triangle |

[`temporal_motif_count`](@ref) counts any motif, given as a list of edges on the
nodes `1:k` in temporal order:

```jldoctest motifs
julia> temporal_motif_count(g, [(1, 2), (2, 1)], 30)          # message and reply
1

julia> temporal_motif_count(g, [(1, 2), (1, 3), (3, 2)], 30)  # triangle
1
```

**Equal time stamps.** The paper requires strictly increasing times, which is the
default `strict = true`. With `strict = false`, edges with the same time stamp are
ordered as in `edges(g)`. SNAP also counts tied edges, but orders them differently in
its star, triangle and 2-node subroutines. On data with ties it therefore gives the
same total as `strict = false`, but slightly different counts in single cells. Only
departure times are used, and self loops are ignored.

**Algorithms.**

- **Stars and 2-node motifs:** one sweep over the edges of each node, keeping counts
  of the edges in the time windows before and after the current edge, and of the
  pairs of edges with the same neighbor (Algorithms 2 and 3 of the paper). The pairs
  are also kept per neighbor, which gives the 2-node motifs in the same pass. This
  takes ``O(m)`` time.
- **Triangles:** every triangle of the static graph is assigned to its pair of nodes
  with the most edges. One sweep per pair counts all the triangles assigned to it
  (Algorithms 4 and 5), in ``O(m\sqrt{τ})`` time, where ``τ`` is the number of static
  triangles. The edges of every pair are collected in time order by a single pass
  over the edge list, without sorting.
- **General motifs:** [`temporal_motif_count`](@ref) enumerates the embeddings of the
  static graph of the motif. For each embedding, it counts the subsequences of its
  edges that spell the motif with a sliding-window dynamic program. This works for
  any motif, but can be slow for motifs with many embeddings, such as 4-cycles in
  dense graphs.

Nodes, pairs and embeddings are processed in parallel. On the CollegeMsg data, the
counts agree exactly with SNAP for edges with distinct time stamps.

## Randomized reference models

Is a pattern seen in a temporal network, such as a motif count, a centrality or a
spreading time, explained by the network's static topology, by the number of events
on each link, or by the burstiness of the links? Randomized reference models answer
this question. Each model generates random networks that keep some features of the
input and are otherwise maximally random. Gauvin et al. (2022) classify the models
used in the literature:

- a *microcanonical* model ``P[x]`` returns a network chosen uniformly at random
  among those with the same features ``x`` as the input;
- the models are ordered by how much they randomize.

[`randomize`](@ref) draws a network from a model:

| Model | Canonical name | Common name |
|---|---|---|
| [`LinkShuffling`](@ref) | ``P[p_\mathcal{L}(Θ)]`` | link shuffling |
| [`DegreeLinkShuffling`](@ref) | ``P[\mathbf{k}, p_\mathcal{L}(Θ)]`` | degree-constrained link shuffling |
| [`TopologyLinkShuffling`](@ref) | ``P[\mathcal{L}, p_\mathcal{L}(Θ)]`` | topology-constrained link shuffling |
| [`WeightLinkShuffling`](@ref) | ``P[\mathbf{w}, p_\mathcal{L}(Θ)]`` | weight-constrained link shuffling |
| [`RandomTimes`](@ref) | ``P[\mathbf{w}]`` | random times |
| [`InterEventShuffling`](@ref) | ``P[π_\mathcal{L}(Δτ), \mathbf{t}_1]`` | inter-event shuffling |
| [`TimelineShifting`](@ref) | ``P[\mathrm{per}(Θ)]`` | timeline shifting |
| [`TimestampShuffling`](@ref) | ``P[\mathbf{w}, \mathbf{t}]`` | timestamp shuffling |
| [`SequenceShuffling`](@ref) | ``P[p_\mathcal{T}(Γ), χ_{\mathbb{N}^+}(\mathbf{A})]`` | sequence shuffling |
| [`SnapshotShuffling`](@ref) | ``P[\mathbf{t}]`` | snapshot shuffling |
| [`DegreeSnapshotShuffling`](@ref) | ``P[\mathbf{d}]`` | degree-constrained snapshot shuffling |
| [`IsomorphicSnapshotShuffling`](@ref) | ``P[\mathrm{iso}(Γ)]`` | isomorphic snapshot shuffling |
| [`EventShuffling`](@ref) | ``P[p(τ)]`` | event shuffling |

The names use these terms:

- an event is an edge ``(u, v, t, τ)``, where ``τ`` is its duration;
- a *link* is a pair of nodes with at least one event;
- the *timeline* of a link is the sequence of its events;
- a *snapshot* is the set of events at one time;
- ``\mathbf{w}`` are the numbers of events per link;
- ``\mathbf{k}`` are the degrees in the static graph;
- ``\mathbf{d}`` are the degrees in each snapshot.

The four groups of models randomize different things:

- **Link shufflings** change the topology but keep every timeline intact.
- **Timeline shufflings** keep the static graph and move the events in time.
- **Sequence and snapshot shufflings** permute or rewire the snapshots.

Durations always move with their events. New start times are drawn from the time
grid of the observation interval [`time_interval`](@ref)`(g)`. For integer times, the
grid step is the greatest common divisor of the time stamps; for floating point
times, the times are continuous. Pass `resolution` to set the step.

Applying a link shuffling and then a timeline shuffling gives another microcanonical
model; so does a sequence shuffling followed by a snapshot shuffling (Propositions
V.5 and V.6 of the paper). For example, `randomize(g, TopologyLinkShuffling(),
TimestampShuffling())` is ``P[\mathcal{L}, p(\mathbf{w}), \mathbf{t}]``. It keeps the
static graph, the distribution of link weights and the times, but destroys the
correlations between weights and topology.

```jldoctest motifs
julia> h = randomize(g, TimestampShuffling());

julia> sort([e.t for e in edges(h)]) == sort([e.t for e in edges(g)])   # same times
true

julia> sort([(e.u, e.v) for e in edges(h)]) == sort([(e.u, e.v) for e in edges(g)])   # same links and weights
true
```

[`reference_samples`](@ref) evaluates a statistic on many random networks in
parallel. With a seeded generator, the result is reproducible. A typical use is the
z-scores of the motif counts against a null model:

```julia
using Statistics, Random
M = temporal_motif_counts(g, 3600)
null = reference_samples(h -> temporal_motif_counts(h, 3600), g, TimestampShuffling();
                         samples = 200, rng = Xoshiro(1))
z = (M .- mean(null)) ./ std(null)
```

**Exactness.** All models except two sample uniformly at random from the networks
with the features they keep, with distinct events. The test suite checks this by a
χ² test over all outcomes on small graphs. The two exceptions run a Markov chain
over the networks and are close to uniform after enough swaps:

- the degree-constrained models (`swaps` double edge swaps per link or event);
- `TimestampShuffling(simple = true)`.

The usual implementation of [`TimestampShuffling`](@ref), a random permutation of the
time stamps, can put two events with the same time on the same link. Such networks
are less likely than the others. `TimestampShuffling(simple = true)` avoids this:
it exchanges times between events and rejects the exchanges that would repeat an
event.

**Undirected networks.** Links are ordered pairs of nodes. For undirected data,
store every contact once and pass `directed = false` to the models, so that links
are unordered pairs and degrees are undirected. Then apply
[`make_undirected`](@ref) to the result if needed.

## References

- A. Paranjape, A. R. Benson, and J. Leskovec. *Motifs in temporal networks.* WSDM, 2017. [DOI](https://doi.org/10.1145/3018661.3018731), [arXiv](https://arxiv.org/abs/1612.09259)
- L. Gauvin, M. Génois, M. Karsai, M. Kivelä, T. Takaguchi, E. Valdano, and C. L. Vestergaard. *Randomized reference models for temporal networks.* SIAM Review 64(4), 2022. [DOI](https://doi.org/10.1137/19M1242252), [arXiv](https://arxiv.org/abs/1806.04032)
- S. Maslov and K. Sneppen. *Specificity and stability in topology of protein networks.* Science 296(5569), 2002. [DOI](https://doi.org/10.1126/science.1065103), [arXiv](https://arxiv.org/abs/cond-mat/0205380)
