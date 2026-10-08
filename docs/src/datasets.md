# Datasets and machine learning

## Public temporal networks

[`load_dataset`](@ref) loads public temporal networks from SNAP and SocioPatterns.

- **Download and cache.** On first use the file is downloaded from the original
  source and cached in [`dataset_dir`](@ref). Set the environment variable
  `TEMPORALGRAPHS_DATA` to change the cache directory.
- **No redistribution.** The data are not redistributed. Every entry of
  [`temporal_datasets`](@ref) records the license of the data and the paper to cite,
  and both are printed when the file is downloaded.

```julia
using TemporalGraphs
temporal_datasets()                       # the catalog
g = load_dataset("CollegeMsg")            # 1,899 nodes, 59,835 edges
h = load_dataset("sp-workplace-2015")     # face-to-face contacts, both directions
snapshots(h)                              # 20 s snapshots
```

Node ids are renumbered `1:n`, and [`original_ids`](@ref) gives the ids of the
source. Every edge gets the transition time `transition_time`, which defaults to 1.
SocioPatterns contacts are undirected, so they are loaded in both directions.

| Name | Source | Directed | Description | License |
|:---|:---|:---|:---|:---|
| `CollegeMsg` | SNAP | yes | Private messages on an online social network at the University of California, Irvine (1,899 nodes, 59,835 edges, 193 days). | SNAP |
| `email-Eu-core-temporal` | SNAP | yes | Emails between members of a European research institution (986 nodes, 332,334 edges, 803 days). | SNAP |
| `email-Eu-core-temporal-Dept1` | SNAP | yes | Emails inside department 1 of a European research institution. | SNAP |
| `email-Eu-core-temporal-Dept2` | SNAP | yes | Emails inside department 2 of a European research institution. | SNAP |
| `email-Eu-core-temporal-Dept3` | SNAP | yes | Emails inside department 3 of a European research institution. | SNAP |
| `email-Eu-core-temporal-Dept4` | SNAP | yes | Emails inside department 4 of a European research institution. | SNAP |
| `sx-mathoverflow` | SNAP | yes | Interactions on the Stack Exchange site mathoverflow, all interactions. | SNAP |
| `sx-mathoverflow-a2q` | SNAP | yes | Interactions on the Stack Exchange site mathoverflow, answers to questions. | SNAP |
| `sx-mathoverflow-c2q` | SNAP | yes | Interactions on the Stack Exchange site mathoverflow, comments to questions. | SNAP |
| `sx-mathoverflow-c2a` | SNAP | yes | Interactions on the Stack Exchange site mathoverflow, comments to answers. | SNAP |
| `sx-askubuntu` | SNAP | yes | Interactions on the Stack Exchange site askubuntu, all interactions. | SNAP |
| `sx-askubuntu-a2q` | SNAP | yes | Interactions on the Stack Exchange site askubuntu, answers to questions. | SNAP |
| `sx-askubuntu-c2q` | SNAP | yes | Interactions on the Stack Exchange site askubuntu, comments to questions. | SNAP |
| `sx-askubuntu-c2a` | SNAP | yes | Interactions on the Stack Exchange site askubuntu, comments to answers. | SNAP |
| `sx-superuser` | SNAP | yes | Interactions on the Stack Exchange site superuser, all interactions. | SNAP |
| `sx-superuser-a2q` | SNAP | yes | Interactions on the Stack Exchange site superuser, answers to questions. | SNAP |
| `sx-superuser-c2q` | SNAP | yes | Interactions on the Stack Exchange site superuser, comments to questions. | SNAP |
| `sx-superuser-c2a` | SNAP | yes | Interactions on the Stack Exchange site superuser, comments to answers. | SNAP |
| `sx-stackoverflow` | SNAP | yes | Interactions on the Stack Exchange site stackoverflow, all interactions. | SNAP |
| `sx-stackoverflow-a2q` | SNAP | yes | Interactions on the Stack Exchange site stackoverflow, answers to questions. | SNAP |
| `sx-stackoverflow-c2q` | SNAP | yes | Interactions on the Stack Exchange site stackoverflow, comments to questions. | SNAP |
| `sx-stackoverflow-c2a` | SNAP | yes | Interactions on the Stack Exchange site stackoverflow, comments to answers. | SNAP |
| `wiki-talk-temporal` | SNAP | yes | Edits of user talk pages of the English Wikipedia (1,140,149 nodes, 7,833,140 edges). | SNAP |
| `soc-sign-bitcoin-otc` | SNAP | yes | Who-trusts-whom network of the Bitcoin OTC platform (5,881 nodes, 35,592 rated edges; ratings are ignored). | SNAP |
| `soc-sign-bitcoin-alpha` | SNAP | yes | Who-trusts-whom network of the Bitcoin Alpha platform (3,783 nodes, 24,186 rated edges; ratings are ignored). | SNAP |
| `sp-workplace-2013` | SocioPatterns | no | Face-to-face contacts in an office building in France, 2013 (92 nodes, 20 s resolution). | CC0 |
| `sp-workplace-2015` | SocioPatterns | no | Face-to-face contacts in an office building in France, 2015 (217 nodes, 20 s resolution). | CC0 |
| `sp-sfhh` | SocioPatterns | no | Face-to-face contacts at the SFHH conference in Nice, 2009 (403 nodes, 20 s resolution). | CC0 |
| `sp-primary-school` | SocioPatterns | no | Face-to-face contacts in a primary school in Lyon, 2009 (242 nodes, 20 s resolution). | CC BY-NC-SA |
| `sp-high-school-2011` | SocioPatterns | no | Face-to-face contacts in a high school in Marseille, 2011 (126 nodes). | CC BY-NC-SA |
| `sp-high-school-2012` | SocioPatterns | no | Face-to-face contacts in a high school in Marseille, 2012 (180 nodes). | CC BY-NC-SA |
| `sp-hospital` | SocioPatterns | no | Face-to-face contacts in a hospital ward in Lyon, 2010 (75 nodes, 20 s resolution). | CC BY-NC-SA |
| `sp-hypertext-2009` | SocioPatterns | no | Face-to-face contacts at the ACM Hypertext 2009 conference (113 nodes, 20 s resolution). | CC BY-NC-SA |

The SocioPatterns data are licensed CC0 or CC BY-NC-SA (non-commercial use); see
the `license` and `citation` fields. The node counts and edge counts of all
datasets, except the largest Stack Overflow and Wikipedia ones, were checked against
the numbers published by the sources.

## GraphNeuralNetworks.jl

After `using GNNGraphs` (or `GraphNeuralNetworks`), the conversions to and from
GNNGraphs.jl are available.

- **Continuous time.** `GNNGraph(g)` is the multigraph of all temporal edges, with
  their times as edge features `t` and `tt`. `OrderedEdgeList(h::GNNGraph)` converts
  back.
- **Discrete time.** `TemporalSnapshotsGNNGraph(g; resolution)` returns one
  `GNNGraph` per snapshot, with the times in `tgdata.times`.
  `OrderedEdgeList(tg::TemporalSnapshotsGNNGraph)` converts back.
- **Event graphs.** [`augmented_event_gnngraph`](@ref) is the augmented event graph
  of Heeg et al. (see [Temporal graph isomorphisms](@ref)), with one-hot node
  features. A directed GNN on this graph is as expressive as the temporal
  Weisfeiler–Leman test.

```@example
using TemporalGraphs # hide
draw_timelines(OrderedEdgeList(4, [(1, 2, 0, 1), (2, 3, 0, 1), (1, 3, 5, 1), (2, 1, 7, 1)])) # hide
```

```jldoctest gnn
julia> using TemporalGraphs, GNNGraphs

julia> g = OrderedEdgeList(4, [(1, 2, 0, 1), (2, 3, 0, 1), (1, 3, 5, 1), (2, 1, 7, 1)]);

julia> h = GNNGraph(g);

julia> h.num_edges, h.edata.t
(4, [0, 0, 5, 7])

julia> edges(OrderedEdgeList(h)) == edges(g)
true

julia> tg = TemporalSnapshotsGNNGraph(g; resolution = 5);

julia> tg.num_snapshots, tg.num_edges, tg.tgdata.times
(2, [2, 2], [0, 5])

julia> a = augmented_event_gnngraph(g; β = 2);

julia> a.num_nodes, size(a.ndata.x)
(8, (2, 8))
```

## MLDatasets.jl

After `using MLDatasets`, the graph containers of MLDatasets.jl can be converted.

- `OrderedEdgeList(tg::MLDatasets.TemporalSnapshotsGraph, times = 1:tg.num_snapshots)`
  converts a sequence of snapshots, for example the 1,000 brain networks of
  `TemporalBrains()` (a 1.6 GB download).
  `MLDatasets.TemporalSnapshotsGraph(g; resolution)` converts the other way.
- `OrderedEdgeList(g::MLDatasets.Graph; time = :timestamp)` reads a static graph
  whose edges carry their times as edge data.
- `OrderedEdgeList(h::MLDatasets.HeteroGraph, relation; time = :timestamp)` converts
  one relation of a heterogeneous graph. For example,
  `("user", "rating", "movie")` of `MovieLens` gives a bipartite temporal graph, with
  the source nodes numbered first.

```julia
using TemporalGraphs, MLDatasets, DataFrames, CSV   # MovieLens needs DataFrames and CSV
g = OrderedEdgeList(MovieLens("100k").graphs[1], ("user", "rating", "movie"))
# 2,625 nodes (943 users and 1,682 movies), 100,000 ratings
b = OrderedEdgeList(TemporalBrains()[1])             # 102 nodes, 27 snapshots
```
