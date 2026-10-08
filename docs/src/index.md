```@raw html
<div class="tg-hero">
<img class="tg-hero-logo" src="assets/logo.svg" alt="TemporalGraphs.jl logo">
<div>
<h1 class="tg-hero-title">TemporalGraphs.jl</h1>
<p class="tg-hero-tagline">Fast temporal graph analysis Julia</p>
```

Paths, centralities, components, flows, motifs, null models and generators for networks
whose edges appear over time: messages between users, contacts between people, flights
between airports, trades between accounts.

```@raw html
<div class="tg-hero-buttons">
```

[Get started](tutorial.md) [API reference](api.md) [GitHub](https://github.com/aurorarossi/TemporalGraphs.jl)

```@raw html
</div></div></div>
```

A **temporal graph** is a graph whose edges carry times. In TemporalGraphs.jl the nodes
`1:n` are fixed, and every temporal edge `(u, v, t, tt)` leaves `u` at time `t` and
reaches `v` at time `t + tt`. Temporal paths must follow time direction, and this
changes a lot with respect to static graphs: reachability is no longer transitive, connected components can overlap, and
several problems become NP-hard.

## Features

```@raw html
<div class="tg-cards">
<div class="tg-card">
```

**[Paths and distances](distances.md)**

Earliest arrival, latest departure, fastest, minimum transition time and minimum hop
paths, in time windows and with waiting constraints.

```@raw html
</div>
<div class="tg-card">
```

**[Centralities](centralities.md)**

Closeness (exact, top-k and sampled), betweenness of walks and paths, Katz, PageRank and
walk centrality.

```@raw html
</div>
<div class="tg-card">
```

**[Connectivity](theory.md)**

Reachability, temporal components and their variants, dissemination times, seed
selection and static expansions.

```@raw html
</div>
<div class="tg-card">
```

**[Flows and spanners](@ref "Flows, cuts and disjoint paths")**

Temporal flows and cuts, disjoint paths, separators, spanning trees and spanners.

```@raw html
</div>
<div class="tg-card">
```

**[Motifs and reference models](motifs.md)**

Counts of δ-temporal motifs and randomized reference models of temporal networks.

```@raw html
</div>
<div class="tg-card">
```

**[Generators and parameters](motifs.md#Generators)**

Random temporal graphs and classic constructions; [parameters](@ref "Temporal graph
parameters") such as the temporality and the interval-membership widths.

```@raw html
</div>
<div class="tg-card">
```

**[Datasets and machine learning](datasets.md)**

SNAP and SocioPatterns networks on demand, MLDatasets.jl and GraphNeuralNetworks.jl.

```@raw html
</div>
<div class="tg-card">
```

**[Isomorphisms](isomorphisms.md)**

Time-respecting path isomorphism, the Weisfeiler–Leman test and kernels for temporal
graph classification.

```@raw html
</div>
</div>
```

All of this runs on several [graph representations](representations.md), with
allocation-free kernels and multithreading: see the [Performance tips](performance.md).

!!! info "Made with Claude"
    The code of TemporalGraphs.jl was written by Claude Opus 5.5, with the help of a
    human.


## Installation

TemporalGraphs.jl requires Julia 1.10 or later. It is not registered yet, so install it
from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/aurorarossi/TemporalGraphs.jl")
```

## Quick example

The temporal graph below has 4 nodes and 7 edges, departing at the times 1 to 8. The
first picture labels every arrow with its `(t, tt)`; the second one draws a timeline for
every node, and every edge goes from the timeline of `u` at time `t` to the timeline of
`v` at time `t + tt`:

```@raw html
<div class="tg-figure">
<img src="assets/example_graph.svg" alt="The temporal graph of the quick example" width="300">
<img src="assets/example_timelines.svg" alt="The timelines of the temporal graph of the quick example" width="460">
</div>
```

```jldoctest quick
julia> using TemporalGraphs

julia> g = OrderedEdgeList(4, [(1, 4, 1, 5), (1, 2, 2, 1), (1, 2, 5, 2), (3, 2, 6, 1),
                               (4, 3, 6, 2), (2, 4, 7, 2), (4, 3, 8, 1)])
OrderedEdgeList{Int32, Int64} with 4 nodes, 7 temporal edges, time interval (1, 9)

julia> earliest_arrival_times(g, 1)      # when can node 1 reach the others?
4-element TemporalDistances{Int64}:
 0
 3
 8
 6

julia> minimum_durations(g, 1)           # how long do the fastest trips take?
4-element TemporalDistances{Int64}:
 0
 1
 7
 4

julia> minimum_duration_path(g, 1, 4)    # and which edges do they use?
2-element Vector{TemporalEdge{Int32, Int64}}:
 (1 2 5 2)
 (2 4 7 2)

julia> temporal_closeness(g, Fastest())   # sum of 1 / duration over reachable nodes
4-element Vector{Float64}:
 1.3928571428571428
 0.5
 1.3333333333333333
 1.0
```

Head to [Getting started](tutorial.md) for a tour of the package.

## Citing

If you use TemporalGraphs.jl in your work, please cite:

```bibtex
@software{aurorarossi2026temporalgraphs,
  author    = {Aurora Rossi},
  title     = {TemporalGraphs.jl: Fast Temporal Graph Analysis in Julia},
  year      = {2026},
  version   = {0.1.0},
  url       = {https://github.com/aurorarossi/TemporalGraphs.jl}
}
```

TemporalGraphs.jl was inspired by [TGLib](https://gitlab.com/tgpublic/tglib), the
C++/Python library of Lutz Oettershagen and Petra Mutzel.
