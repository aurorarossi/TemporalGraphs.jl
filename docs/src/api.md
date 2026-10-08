# API reference

Every exported function and type, grouped as in the manual. An [alphabetical index](@ref
"Index") is at the end of the page.

```@meta
CurrentModule = TemporalGraphs
```

```@docs
TemporalGraphs
```

## Temporal edges and basic types

```@docs
TemporalEdge
TimeInterval
INF
TemporalDistances
StaticWeightedEdge
DistanceType
EarliestArrival
Fastest
LatestDeparture
MinimumTransitionTimes
MinimumHops
ShortestForemost
ShortestFastest
ShortestLatest
PrefixForemost
```

## Graph representations

```@docs
OrderedEdgeList
IncidentLists
TRSGraph
DirectedLineGraph
num_nodes
num_edges
edges
time_interval
time_type
node_type
out_edges
out_degree
original_id
original_ids
node_map
num_trs_nodes
trs_node
trs_time
trs_neighbors
```

## Input and output

```@docs
load_ordered_edge_list
load_incident_lists
load_trs_graph
save_ordered_edge_list
OrderedEdgeList(::Any)
```

## Transformations

```@docs
to_ordered_edge_list
to_incident_lists
to_trs_graph
to_directed_line_graph
to_aggregated_edge_list
static_graph
snapshots
OrderedEdgeList(::AbstractVector{<:Graphs.AbstractGraph}, ::AbstractVector{<:Real})
aggregate_time
normalize_graph
scale_timestamps
unit_transition_times
make_undirected
Base.reverse(::OrderedEdgeList, ::Any)
```

## Drawing

```@docs
draw_graph
draw_timelines
TemporalGraphDrawing
```

## Distances

```@docs
temporal_distances
temporal_distances!
distance_workspace
distance_eltype
earliest_arrival_times
latest_departure_times
minimum_durations
minimum_transition_times
minimum_hops
number_of_reachable_nodes
StreamWorkspace
IncidentListsWorkspace
TRSWorkspace
```

## Paths

```@docs
earliest_arrival_path
minimum_duration_path
minimum_transition_time_path
minimum_hops_path
```

## Closeness, eccentricity, diameter and efficiency

```@docs
temporal_closeness
compute_topk_closeness
temporal_closeness_approximation
temporal_harmonic_closeness
temporal_harmonic_closeness_approximation
temporal_eccentricity
temporal_diameter
temporal_efficiency
temporal_distance_statistics
```

## Temporal betweenness

```@docs
temporal_betweenness
temporal_ego_betweenness
temporal_pass_through_degree
temporal_betweenness_approximation
temporal_betweenness_mantra
temporal_betweenness_atbc
temporal_edge_betweenness
```

## Walk-based centralities

```@docs
temporal_katz_centrality
temporal_pagerank
temporal_walk_centrality
```

## Statistics and parameters

```@docs
TemporalGraphStatistics
get_statistics
vertex_interval_membership_width
edge_interval_membership_width
is_simple
is_proper
edge_burstiness
node_burstiness
temporal_clustering_coefficient
topological_overlap
kcores
temporal_khcores
temporal_lkcores
```

## Reachability and dissemination

```@docs
temporal_reachability
is_temporally_connected
earliest_arrival_matrix
temporal_flooding_time
temporal_flooding_times
temporal_gossip_time
```

## Temporal components

```@docs
temporal_connected_components
largest_temporal_connected_component
source_component
sink_component
window_components
persistent_components
interval_connected_components
delta_temporal_connected_components
stream_components
```

## Seed selection

```@docs
reachability_dominating_set
max_reach_seeds
```

## Static expansions and restless paths

```@docs
vertex_expansion
edge_expansion
restless_path
```

## Flows, cuts and separators

```@docs
temporal_max_flow
temporal_min_cut
temporal_edge_disjoint_paths
temporal_out_disjoint_paths
temporal_vertex_separator
```

## Spanning trees and spanners

```@docs
earliest_arrival_tree
minimum_weight_spanning_tree
temporal_spanner
temporal_clique_spanner
```

## Isomorphisms

```@docs
augmented_event_graph
temporal_isomorphism
is_temporally_isomorphic
temporal_wl_equivalent
temporal_wl_kernel
augmented_event_gnngraph
```

## Motifs

```@docs
temporal_motif_counts
temporal_motif_count
```

## Randomized reference models

```@docs
randomize
reference_samples
ReferenceModel
LinkShuffling
DegreeLinkShuffling
TopologyLinkShuffling
WeightLinkShuffling
RandomTimes
InterEventShuffling
TimelineShifting
TimestampShuffling
SequenceShuffling
SnapshotShuffling
DegreeSnapshotShuffling
IsomorphicSnapshotShuffling
EventShuffling
```

## Generators

```@docs
random_temporal_graph
random_simple_temporal_graph
random_temporal_labeling
round_robin_temporal_clique
temporal_hypercube
temporal_knodel_graph
```

## Datasets

```@docs
TemporalDataset
temporal_datasets
load_dataset
dataset_dir
```

## Index

```@index
Pages = ["api.md"]
```
