"""
    TemporalGraphs

A Julia library for temporal graph analysis: temporal distances and paths,
temporal centralities (closeness, betweenness, Katz, PageRank, walk centrality),
connectivity, flows, spanners, motifs, isomorphisms, randomized reference models and
other local and global temporal graph statistics.

The package was inspired by TGLib (<https://gitlab.com/tgpublic/tglib>, Lutz
Oettershagen), whose functionality it covers with its own Julia design and
algorithms.

The graph types are parametric in the node id type `V` and the time type `T`
(integers or floating point numbers). Node ids are 1-based. Every temporal edge
`(u, v, t, tt)` goes from `u` to `v`,
departs at time `t` and arrives at time `t + tt`. Algorithms can be restricted to
a time interval `(a, b)`: an edge is usable iff `a ≤ t` and `t + tt ≤ b`.

# References

- L. Oettershagen and P. Mutzel. *TGLib: an open-source library for temporal graph analysis.* IEEE International Conference on Data Mining Workshops (ICDMW), 2022. [DOI](https://doi.org/10.1109/ICDMW58026.2022.00160), [arXiv](https://arxiv.org/abs/2209.12587)
"""
module TemporalGraphs

using DataStructures: BinaryMinHeap
using OhMyThreads: @tasks, @local, tmap, tmapreduce, chunks
using PrecompileTools: @compile_workload, @setup_workload
using Random: Random, AbstractRNG
using Statistics: mean, std
import Graphs
import Graphs: edges, nv, ne
import Tables
import Downloads
import p7zip_jll

export
    # basic types
    TemporalEdge, TimeInterval, StaticWeightedEdge, INF, time_type, node_type,
    DistanceType, EarliestArrival, Fastest, LatestDeparture, MinimumTransitionTimes, MinimumHops,
    ShortestForemost, ShortestFastest, ShortestLatest, PrefixForemost,
    # representations
    OrderedEdgeList, IncidentLists, TRSGraph, DirectedLineGraph,
    num_nodes, num_edges, edges, time_interval, out_edges, out_degree,
    original_id, original_ids, node_map,
    num_trs_nodes, trs_neighbors, trs_node, trs_time,
    # transformations
    to_ordered_edge_list, to_incident_lists, to_trs_graph, to_directed_line_graph,
    to_aggregated_edge_list, normalize_graph, scale_timestamps, unit_transition_times,
    make_undirected, static_graph, snapshots, aggregate_time,
    # statistics and IO
    TemporalGraphStatistics, get_statistics,
    load_ordered_edge_list, load_incident_lists, load_trs_graph, save_ordered_edge_list,
    # distances and paths
    earliest_arrival_times, latest_departure_times, minimum_durations, minimum_hops,
    minimum_transition_times, temporal_distances, temporal_distances!, distance_workspace,
    distance_eltype, number_of_reachable_nodes,
    earliest_arrival_path, minimum_duration_path, minimum_transition_time_path, minimum_hops_path,
    # centralities and statistics
    temporal_closeness, compute_topk_closeness, temporal_closeness_approximation,
    temporal_harmonic_closeness, temporal_harmonic_closeness_approximation,
    temporal_motif_counts, temporal_motif_count,
    ReferenceModel, LinkShuffling, DegreeLinkShuffling, TopologyLinkShuffling, WeightLinkShuffling,
    RandomTimes, InterEventShuffling, TimelineShifting, TimestampShuffling, SequenceShuffling,
    SnapshotShuffling, DegreeSnapshotShuffling, IsomorphicSnapshotShuffling, EventShuffling,
    randomize, reference_samples,
    temporal_reachability, restless_path, is_temporally_connected, temporal_connected_components,
    largest_temporal_connected_component,
    temporal_max_flow, temporal_min_cut, temporal_edge_disjoint_paths, temporal_out_disjoint_paths,
    temporal_vertex_separator,
    earliest_arrival_tree, minimum_weight_spanning_tree, temporal_spanner, temporal_clique_spanner,
    earliest_arrival_matrix, temporal_flooding_time, temporal_flooding_times, temporal_gossip_time,
    vertex_expansion, edge_expansion, source_component, sink_component, window_components,
    persistent_components, interval_connected_components,
    augmented_event_graph, temporal_isomorphism, is_temporally_isomorphic, temporal_wl_equivalent,
    temporal_wl_kernel, augmented_event_gnngraph,
    TemporalDataset, temporal_datasets, load_dataset, dataset_dir,
    temporal_eccentricity, temporal_diameter, temporal_efficiency,
    temporal_edge_betweenness, temporal_katz_centrality, temporal_pagerank,
    temporal_betweenness, temporal_ego_betweenness, temporal_pass_through_degree,
    temporal_betweenness_approximation, temporal_betweenness_mantra, temporal_betweenness_atbc, temporal_distance_statistics,
    temporal_walk_centrality, edge_burstiness, node_burstiness,
    temporal_clustering_coefficient, topological_overlap,
    kcores, temporal_khcores, temporal_lkcores

include("core/types.jl")
include("core/ordered_edge_list.jl")
include("core/incident_lists.jl")
include("core/trs_graph.jl")
include("core/directed_line_graph.jl")
include("core/transformations.jl")
include("core/statistics.jl")
include("core/snapshots.jl")
include("io.jl")

include("util/heap.jl")
include("util/topk.jl")
include("util/maxflow.jl")

include("algorithms/distances_stream.jl")
include("algorithms/distances_incident.jl")
include("algorithms/distances_trs.jl")
include("algorithms/distances.jl")
include("algorithms/paths.jl")
include("algorithms/closeness.jl")
include("algorithms/topk_closeness.jl")
include("algorithms/harmonic_closeness.jl")
include("algorithms/diameter_efficiency.jl")
include("algorithms/betweenness.jl")
include("algorithms/temporal_betweenness.jl")
include("algorithms/betweenness_proxies.jl")
include("algorithms/betweenness_sampling.jl")
include("algorithms/walks.jl")
include("algorithms/burstiness.jl")
include("algorithms/local_statistics.jl")
include("algorithms/cores.jl")
include("algorithms/motifs.jl")
include("algorithms/reference_models.jl")
include("algorithms/connectivity.jl")
include("algorithms/restless.jl")
include("algorithms/flows.jl")
include("algorithms/trees.jl")
include("algorithms/expansions.jl")
include("algorithms/components.jl")
include("algorithms/isomorphism.jl")
include("datasets.jl")

"""
    augmented_event_gnngraph(g::OrderedEdgeList, ti = time_interval(g); β = Inf)

The [`augmented_event_graph`](@ref) of `g` as a `GNNGraph` of GNNGraphs.jl, with node
data `label` (0 for the nodes of `g`, 1 for its temporal edges) and `x`, the one-hot
encoding of the labels: the input of the event graph neural network of Heeg, Sauer,
Mutzel and Scholtes (2025), e.g. a directed GNN of GraphNeuralNetworks.jl. Available
after `using GNNGraphs` (or `GraphNeuralNetworks`).
"""
function augmented_event_gnngraph end
include("precompile.jl")

end # module
