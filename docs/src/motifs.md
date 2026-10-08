# Motifs, reference models and generators

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

julia> g = OrderedEdgeList(3, [(1, 2, 10, 0), (2, 1, 20, 0), (1, 3, 25, 0),
                               (3, 2, 40, 0), (1, 2, 70, 0)]);

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
|:---|:---|:---|
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
  (Algorithms 4 and 5), in ``O(m\sqrt{T_△})`` time, where ``T_△`` is the number of
  triangles of the static graph. The edges of every pair are collected in time order by a single pass
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
on each link, or by the burstiness of the links? Randomized reference models (*null
models*) answer this question. Each model generates random networks that keep some features of the
input and are otherwise maximally random. Gauvin et al. (2022) classify the models
used in the literature: a *microcanonical* model ``P[x]`` returns a network chosen
uniformly at random among those with the same features ``x`` as the input, and the
models are ordered by how much they randomize.

[`randomize`](@ref) draws a network from a model:

| Model | Canonical name | Common name |
|:---|:---|:---|
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

- an *event* is an edge `(u, v, t, tt)`, whose transition time `tt` is its duration;
- a *link* is a pair of nodes with at least one event;
- the *timeline* of a link is the sequence of its events;
- a *snapshot* is the set of events at one time;
- ``\mathbf{w}`` are the numbers of events per link;
- ``\mathbf{k}`` are the degrees in the static graph;
- ``\mathbf{d}`` are the degrees in each snapshot.

The groups of models randomize different things:

- **Link shufflings** change the topology but keep every timeline intact.
- **Timeline shufflings** keep the static graph and move the events in time.
- **Sequence and snapshot shufflings** permute or rewire the snapshots.
- **Event shuffling**, the most random model, gives every event random nodes and a
  random time, keeping only the durations.

Durations always move with their events. New start times are drawn from the time
grid of the observation interval [`time_interval`](@ref)`(g)`: for integer times the
grid step is the [resolution](@ref "Snapshot graphs") of the time stamps, for floating
point times the times are continuous. The models that draw new times
([`RandomTimes`](@ref), [`TimelineShifting`](@ref), [`SequenceShuffling`](@ref) and
[`EventShuffling`](@ref)) take a `resolution` keyword to set the step.

Applying a link shuffling and then a timeline shuffling gives another microcanonical
model; so does a sequence shuffling followed by a snapshot shuffling (Propositions
V.5 and V.6 of Gauvin et al.). For example, `randomize(g, TopologyLinkShuffling(),
TimestampShuffling())` is ``P[\mathcal{L}, p(\mathbf{w}), \mathbf{t}]``. It keeps the
static graph, the distribution of link weights and the times, but destroys the
correlations between weights and topology.

```jldoctest motifs
julia> h = randomize(g, TimestampShuffling());

julia> sort([e.t for e in edges(h)]) == sort([e.t for e in edges(g)])   # same times
true

julia> links(x) = sort([(e.u, e.v) for e in edges(x)]);

julia> links(h) == links(g)   # same links and weights
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

**Exactness.** The models sample uniformly at random among the networks without
repeated events that have the features they keep; the test suite checks this by a χ²
test over all outcomes on small graphs. The exceptions run a Markov chain over the
networks and are close to uniform after enough swaps:

- the degree-constrained models, with `swaps` double edge swaps (Maslov and Sneppen,
  2002) per link or event;
- `TimestampShuffling(simple = true)`.

The usual implementation of [`TimestampShuffling`](@ref), a random permutation of the
time stamps, can put two events with the same time on the same link. Such networks
are less likely than the others. `TimestampShuffling(simple = true)` avoids this:
it exchanges times between events and rejects the exchanges that would repeat an
event.

**Undirected networks.** Links are ordered pairs of nodes. For undirected data,
store every contact once and pass `directed = false` to the models that have this
keyword, so that links are unordered pairs and degrees are undirected. Then apply
[`make_undirected`](@ref) to the result if needed.

## Generators

Generators build temporal graphs from scratch, for tests, benchmarks and the study of
graph classes (see the *Classes* and *Constructions* pages of the [Temporal Graph
Wiki](https://temporalgraph.notion.site/)). Undirected graphs are stored with both directions of every edge.

- [`random_temporal_graph`](@ref)`(n, p, T)`: `T` independent Erdős–Rényi snapshots
  `G(n, p)` at the times `1, ..., T`.
- [`random_simple_temporal_graph`](@ref)`(n, p)`: a *random simple temporal graph*
  (Casteigts, Raskin, Renken and Zamaraev, 2024), an Erdős–Rényi graph whose edges get
  distinct times in random order; with `p = 1` a random temporal clique.
- [`random_temporal_labeling`](@ref)`(G, T)`: random times for the edges of any
  Graphs.jl graph, e.g. temporal paths, stars, trees and grids.
- [`round_robin_temporal_clique`](@ref)`(n)`: the clique scheduled as a round-robin
  tournament, one matching per time stamp.
- [`temporal_hypercube`](@ref)`(d)` and [`temporal_knodel_graph`](@ref)`(d)`: the
  hypercube (Kempe, Kleinberg and Kumar, 2002) and the Knödel graph on `2^d` nodes with
  the edges of dimension `k` at time `k`. Every pair of nodes is joined by exactly one
  temporal path, so these temporally connected graphs with `d 2^(d-1)` undirected edges
  (`d 2^d` temporal edges, as both directions are stored) have no smaller temporal
  spanner.

```jldoctest
julia> using TemporalGraphs

julia> h = temporal_hypercube(3)
OrderedEdgeList{Int32, Int64} with 8 nodes, 24 temporal edges, time interval (1, 4)

julia> is_temporally_connected(h), length(temporal_spanner(h))   # every edge is needed
(true, 24)
```

The random generators take an `rng` keyword argument for reproducible graphs:

```julia
using Random, Graphs
g = random_temporal_graph(1000, 0.001, 100; rng = MersenneTwister(1))
s = random_simple_temporal_graph(200, 3 * log(200) / 200)   # connectivity threshold
t = random_temporal_labeling(uniform_tree(50), 20)          # random tree with random times
```

## References

- A. Paranjape, A. R. Benson, and J. Leskovec. *Motifs in temporal networks.* WSDM, 2017. [DOI](https://doi.org/10.1145/3018661.3018731), [arXiv](https://arxiv.org/abs/1612.09259)
- L. Gauvin, M. Génois, M. Karsai, M. Kivelä, T. Takaguchi, E. Valdano, and C. L. Vestergaard. *Randomized reference models for temporal networks.* SIAM Review 64(4), 2022. [DOI](https://doi.org/10.1137/19M1242252), [arXiv](https://arxiv.org/abs/1806.04032)
- S. Maslov and K. Sneppen. *Specificity and stability in topology of protein networks.* Science 296(5569), 2002. [DOI](https://doi.org/10.1126/science.1065103), [arXiv](https://arxiv.org/abs/cond-mat/0205380)
- A. Casteigts, M. Raskin, M. Renken, and V. Zamaraev. *Sharp thresholds in random simple temporal graphs.* SIAM Journal on Computing, 2024. [DOI](https://doi.org/10.1137/22M1511916), [arXiv](https://arxiv.org/abs/2011.03738)
- D. Kempe, J. Kleinberg, and A. Kumar. *Connectivity and inference problems for temporal networks.* Journal of Computer and System Sciences 64(4), 2002. [DOI](https://doi.org/10.1006/jcss.2002.1829)
- W. Knödel. *New gossips and telephones.* Discrete Mathematics 13(1), 1975. [DOI](https://doi.org/10.1016/0012-365X(75)90090-4)
