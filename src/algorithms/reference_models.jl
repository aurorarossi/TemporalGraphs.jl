# Microcanonical randomized reference models (MRRMs) of temporal networks, following
# Gauvin, Génois, Karsai, Kivelä, Takaguchi, Valdano and Vestergaard, "Randomized
# reference models for temporal networks" (SIAM Review, 2022). A MRRM returns a
# uniformly random network among those sharing some features with the input; it is
# named P[x] after the features x it keeps (see the docstrings). Every model is
# implemented by a shuffling method.
#
# A temporal edge (u, v, t, tt) is an event between u and v starting at t with
# duration tt; durations always move with their event. The links are the ordered
# pairs (u, v) with at least one event, or the unordered pairs with
# `directed = false`; the timeline of a link is the sequence of its events. The
# snapshots are the sets of events with the same start time. New start times are
# drawn on the time grid a, a + Δ, a + 2Δ, ... of the observation interval
# (a, b) = time_interval(g), so that t + tt ≤ b, where Δ is the resolution (by default
# the gcd of the times minus a for integer times, continuous for floating point
# times).

"""
    ReferenceModel

Supertype of the randomized reference models of temporal networks, see
[`randomize`](@ref).
"""
abstract type ReferenceModel end

"""
    LinkShuffling(; directed = true)

``P[p_\\mathcal{L}(Θ)]``, *link shuffling*: the static graph is replaced by a uniformly
random graph with the same number of links (Erdős–Rényi), and the timelines are
redistributed at random on the new links. Keeps the multiset of timelines (hence the
activity timeline and the inter-event times on links); randomizes the topology.
"""
Base.@kwdef struct LinkShuffling <: ReferenceModel
    directed::Bool = true
end

"""
    DegreeLinkShuffling(; directed = true, swaps = 10)

``P[\\mathbf{k}, p_\\mathcal{L}(Θ)]``, *degree-constrained link shuffling*: the static
graph is randomized keeping the degree of every node (in- and out-degrees if
`directed`), with `swaps` double edge swaps per link (Maslov and Sneppen), and the
timelines are redistributed at random on the links.

# References

- S. Maslov and K. Sneppen. *Specificity and stability in topology of protein networks.* Science 296(5569), 2002. [DOI](https://doi.org/10.1126/science.1065103), [arXiv](https://arxiv.org/abs/cond-mat/0205380)
"""
Base.@kwdef struct DegreeLinkShuffling <: ReferenceModel
    directed::Bool = true
    swaps::Int = 10
end

"""
    TopologyLinkShuffling(; directed = true)

``P[\\mathcal{L}, p_\\mathcal{L}(Θ)]``, *topology-constrained link shuffling*: the
timelines are permuted at random among the links of the static graph, which is
kept. Destroys the correlations between topology and dynamics.
"""
Base.@kwdef struct TopologyLinkShuffling <: ReferenceModel
    directed::Bool = true
end

"""
    WeightLinkShuffling(; directed = true)

``P[\\mathbf{w}, p_\\mathcal{L}(Θ)]``, *weight-constrained link shuffling*: the
timelines are permuted at random among the links with the same number of events.
"""
Base.@kwdef struct WeightLinkShuffling <: ReferenceModel
    directed::Bool = true
end

"""
    RandomTimes(; directed = true, resolution = nothing)

``P[\\mathbf{w}]``, *weight-constrained timeline shuffling* (also called *random
times* or *Poissonized inter-event intervals*): every event gets a uniformly random
start time in the observation interval, keeping its link; the events of a link have
distinct times. Keeps the static graph and the number of events of every link.
"""
Base.@kwdef struct RandomTimes{R} <: ReferenceModel
    directed::Bool = true
    resolution::R = nothing
end

"""
    InterEventShuffling(; directed = true)

``P[π_\\mathcal{L}(Δτ), \\mathbf{t}_1]``, *inter-event shuffling*: on every link the
intervals between the start times of consecutive events are permuted at random,
keeping the time of the first event. Keeps the static graph, the number of events
and the distribution of inter-event times of every link; destroys their order
(e.g. bursts). The durations keep their order on the link.
"""
Base.@kwdef struct InterEventShuffling <: ReferenceModel
    directed::Bool = true
end

"""
    TimelineShifting(; directed = true, resolution = nothing)

``P[\\mathrm{per}(Θ)]``, *timeline shifting*: every timeline is translated by a
uniformly random offset with periodic boundary conditions on the observation
interval. Keeps all intervals between the events of a link, randomizes the relative
timing of the links.
"""
Base.@kwdef struct TimelineShifting{R} <: ReferenceModel
    directed::Bool = true
    resolution::R = nothing
end

"""
    TimestampShuffling(; simple = false, swaps = 10)

``P[\\mathbf{w}, \\mathbf{t}]``, *timestamp shuffling* (also called *permuted times*
or *time reshuffle*): the start times (with the durations) are permuted at random
among all events, keeping their links. Keeps the static graph, the number of events
of every link and the activity timeline (the number of events at every time).

A random permutation, as in the literature, can give a link two events with the same
time; such networks are then less likely than the others. With `simple = true` the
times are instead exchanged between random pairs of events, `swaps` times per event,
rejecting the exchanges that would repeat an event. No repeated event is created
(those of the input disappear as their times are exchanged), and the result is
uniform among the networks without repeated events.
"""
Base.@kwdef struct TimestampShuffling <: ReferenceModel
    simple::Bool = false
    swaps::Int = 10
end

"""
    SequenceShuffling(; active = true, resolution = nothing)

``P[p_\\mathcal{T}(Γ), χ_{\\mathbb{N}^+}(\\mathbf{A})]``, *activity-constrained
sequence shuffling*: the snapshots (events with the same start time) are permuted at
random among the times with at least one event. With `active = false`,
``P[p_\\mathcal{T}(Γ)]``: the snapshots are placed at distinct random times of the time
grid of the observation interval. Keeps every snapshot, destroys the temporal order.
"""
Base.@kwdef struct SequenceShuffling{R} <: ReferenceModel
    active::Bool = true
    resolution::R = nothing
end

"""
    SnapshotShuffling(; directed = true)

``P[\\mathbf{t}]``, *snapshot shuffling*: every snapshot is replaced by a uniformly
random graph with the same number of events (distinct pairs of distinct nodes).
Keeps the times of all events; randomizes the topology.
"""
Base.@kwdef struct SnapshotShuffling <: ReferenceModel
    directed::Bool = true
end

"""
    DegreeSnapshotShuffling(; directed = true, swaps = 10)

``P[\\mathbf{d}]``, *degree-constrained snapshot shuffling*: every snapshot is
randomized keeping the instantaneous degree of every node (in- and out-degrees if
`directed`), with `swaps` double edge swaps per event.

# References

- S. Maslov and K. Sneppen. *Specificity and stability in topology of protein networks.* Science 296(5569), 2002. [DOI](https://doi.org/10.1126/science.1065103), [arXiv](https://arxiv.org/abs/cond-mat/0205380)
"""
Base.@kwdef struct DegreeSnapshotShuffling <: ReferenceModel
    directed::Bool = true
    swaps::Int = 10
end

"""
    IsomorphicSnapshotShuffling()

``P[\\mathrm{iso}(Γ)]``, *isomorphic snapshot shuffling*: the node ids of every
snapshot are permuted at random, independently for every snapshot. Keeps the
structure of every snapshot up to isomorphism.
"""
struct IsomorphicSnapshotShuffling <: ReferenceModel end

"""
    EventShuffling(; directed = true, resolution = nothing)

``P[p(τ)]`` (``P[E]`` for instantaneous events), *event shuffling*: every event gets a
uniformly random pair of distinct nodes and a uniformly random start time, keeping
its duration; all events are distinct. The most random model: it keeps only the
nodes, the observation interval, the number of events and their durations.
"""
Base.@kwdef struct EventShuffling{R} <: ReferenceModel
    directed::Bool = true
    resolution::R = nothing
end

# ----------------------------------------------------------------------------------
# helpers

_link_key(u::Integer, v::Integer, directed::Bool) = directed ? (Int(u), Int(v)) : minmax(Int(u), Int(v))

# The links of es and the indices of their events in time order.
function _timelines(es::Vector{<:TemporalEdge}, directed::Bool)
    index = Dict{Tuple{Int,Int},Int}()
    links = Tuple{Int,Int}[]
    lines = Vector{Vector{Int}}()
    for (i, e) in enumerate(es)
        k = _link_key(e.u, e.v, directed)
        j = get!(index, k) do
            push!(links, k)
            push!(lines, Int[])
            length(links)
        end
        push!(lines[j], i)
    end
    return links, lines
end

# The snapshots of es (sorted by time): ranges of equal start times.
function _snapshots(es::Vector{<:TemporalEdge})
    out = UnitRange{Int}[]
    lo = 1
    for i in eachindex(es)
        if i == length(es) || es[i+1].t != es[i].t
            push!(out, lo:i)
            lo = i + 1
        end
    end
    return out
end

# Time grid (a, b, Δ) of the observation interval; Δ === nothing for continuous times.
function _time_grid(g::OrderedEdgeList{V,T}, resolution) where {V,T}
    a, b = time_interval(g)
    resolution === nothing || return (a, b, convert(T, resolution))
    T <: Integer || return (a, b, nothing)
    Δ = zero(T)
    for e in edges(g)
        Δ = gcd(Δ, e.t - a)
    end
    return (a, b, Δ == 0 ? one(T) : Δ)
end

# Uniformly random start time of an event of duration tt.
function _random_time(rng, (a, b, Δ), tt)
    if Δ === nothing
        return a + rand(rng) * max(b - tt - a, zero(b - a))
    end
    K = max(fld(b - tt - a, Δ), 0)
    return a + Δ * rand(rng, 0:K)
end

# Uniformly random link of n nodes (ordered pair of distinct nodes, sorted if undirected).
function _random_link(rng, n::Int, directed::Bool)
    u = rand(rng, 1:n)
    v = rand(rng, 1:n-1)
    v >= u && (v += 1)
    return directed ? (u, v) : minmax(u, v)
end

function _check_pairs(n::Int, L::Int, directed::Bool)
    n >= 2 || L == 0 || throw(ArgumentError("at least two nodes are needed"))
    L <= (directed ? n * (n - 1) : n * (n - 1) ÷ 2) ||
        throw(ArgumentError("more links than pairs of nodes"))
    return nothing
end

# Event e moved to link (x, y): if undirected, keep its orientation relative to the
# order of the endpoints.
function _move(e::TemporalEdge{V,T}, (x, y)::Tuple{Int,Int}, directed::Bool) where {V,T}
    (directed || e.u <= e.v) && return TemporalEdge{V,T}(x, y, e.t, e.tt)
    return TemporalEdge{V,T}(y, x, e.t, e.tt)
end

_rebuilt(g::OrderedEdgeList, es) = OrderedEdgeList(num_nodes(g), es, time_interval(g); original_ids=g.original_ids)

# Timelines moved to new links: link j of `links` receives the timeline perm[j].
function _move_timelines(g::OrderedEdgeList{V,T}, lines, links, perm, directed) where {V,T}
    es = edges(g)
    out = Vector{TemporalEdge{V,T}}(undef, length(es))
    x = 0
    for (j, l) in enumerate(links), i in lines[perm[j]]
        out[x+=1] = _move(es[i], l, directed)
    end
    return _rebuilt(g, out)
end

# Double edge swaps on a list of links keeping (in- and out-) degrees and simplicity:
# (a, b), (c, d) → (a, d), (c, b), or (a, c), (b, d) for undirected links.
function _swap_links!(rng, links::Vector{Tuple{Int,Int}}, nswaps::Int, directed::Bool,
                      present=Dict{Tuple{Int,Int},Int}(l => 1 for l in links))
    L = length(links)
    L < 2 && return links
    for _ in 1:nswaps
        i, j = rand(rng, 1:L), rand(rng, 1:L)
        i == j && continue
        (a, b), (c, d) = links[i], links[j]
        if !directed && rand(rng, Bool)
            c, d = d, c
        end
        (a == d || c == b) && continue
        n1, n2 = _link_key(a, d, directed), _link_key(c, b, directed)
        (get(present, n1, 0) > 0 || get(present, n2, 0) > 0 || n1 == n2) && continue
        for (old, new) in ((links[i], n1), (links[j], n2))
            present[old] -= 1
            present[new] = get(present, new, 0) + 1
        end
        links[i], links[j] = n1, n2
    end
    return links
end

# ----------------------------------------------------------------------------------
# shufflings

function _randomize(rng::AbstractRNG, g::OrderedEdgeList, m::LinkShuffling)
    links, lines = _timelines(edges(g), m.directed)
    n, L = num_nodes(g), length(links)
    _check_pairs(n, L, m.directed)
    new = Set{Tuple{Int,Int}}()
    while length(new) < L
        push!(new, _random_link(rng, n, m.directed))
    end
    return _move_timelines(g, lines, collect(new), Random.randperm(rng, L), m.directed)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList, m::DegreeLinkShuffling)
    links, lines = _timelines(edges(g), m.directed)
    new = _swap_links!(rng, copy(links), m.swaps * length(links), m.directed)
    return _move_timelines(g, lines, new, Random.randperm(rng, length(links)), m.directed)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList, m::TopologyLinkShuffling)
    links, lines = _timelines(edges(g), m.directed)
    return _move_timelines(g, lines, links, Random.randperm(rng, length(links)), m.directed)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList, m::WeightLinkShuffling)
    links, lines = _timelines(edges(g), m.directed)
    perm = collect(eachindex(links))
    byweight = Dict{Int,Vector{Int}}()
    for j in eachindex(links)
        push!(get!(byweight, length(lines[j]), Int[]), j)
    end
    for js in values(byweight)
        perm[js] = js[Random.randperm(rng, length(js))]
    end
    return _move_timelines(g, lines, links, perm, m.directed)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, m::RandomTimes) where {V,T}
    es = edges(g)
    grid = _time_grid(g, m.resolution)
    _, lines = _timelines(es, m.directed)
    out = similar(es)
    x = 0
    for line in lines
        seen = Set{Tuple{T,T}}()
        for i in line
            e = es[i]
            for attempt in 1:1000
                t = _random_time(rng, grid, e.tt)
                if !((t, e.tt) in seen)
                    push!(seen, (t, e.tt))
                    out[x+=1] = TemporalEdge{V,T}(e.u, e.v, t, e.tt)
                    break
                end
                attempt == 1000 && throw(ArgumentError("a link has more events than times in the observation interval"))
            end
        end
    end
    return _rebuilt(g, out)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, m::InterEventShuffling) where {V,T}
    es = edges(g)
    _, lines = _timelines(es, m.directed)
    out = similar(es)
    x = 0
    for line in lines
        gaps = [es[line[r+1]].t - es[line[r]].t for r in 1:length(line)-1]
        Random.shuffle!(rng, gaps)
        t = es[line[1]].t
        for (r, i) in enumerate(line)
            r > 1 && (t += gaps[r-1])
            out[x+=1] = TemporalEdge{V,T}(es[i].u, es[i].v, t, es[i].tt)
        end
    end
    return _rebuilt(g, out)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, m::TimelineShifting) where {V,T}
    es = edges(g)
    a, b, Δ = _time_grid(g, m.resolution)
    # period: the length of the interval, or the number of grid points times Δ
    P = Δ === nothing ? b - a : (fld(b - a, Δ) + 1) * Δ
    _, lines = _timelines(es, m.directed)
    out = similar(es)
    x = 0
    for line in lines
        s = Δ === nothing ? rand(rng) * P : Δ * rand(rng, 0:fld(P, Δ)-1)
        for i in line
            e = es[i]
            out[x+=1] = TemporalEdge{V,T}(e.u, e.v, a + mod(e.t - a + s, P), e.tt)
        end
    end
    return _rebuilt(g, out)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, m::TimestampShuffling) where {V,T}
    es = edges(g)
    E = length(es)
    if !m.simple
        p = Random.randperm(rng, E)
        return _rebuilt(g, [TemporalEdge{V,T}(e.u, e.v, es[p[i]].t, es[p[i]].tt) for (i, e) in enumerate(es)])
    end
    present = Dict{TemporalEdge{V,T},Int}()  # multiplicities
    for e in es
        present[e] = get(present, e, 0) + 1
    end
    out = copy(es)
    E < 2 && return _rebuilt(g, out)
    for _ in 1:m.swaps*E
        i, j = rand(rng, 1:E), rand(rng, 1:E)
        a, b = out[i], out[j]
        (a.t == b.t && a.tt == b.tt) && continue
        a2 = TemporalEdge{V,T}(a.u, a.v, b.t, b.tt)
        b2 = TemporalEdge{V,T}(b.u, b.v, a.t, a.tt)
        (get(present, a2, 0) > 0 || get(present, b2, 0) > 0) && continue
        present[a] -= 1
        present[b] -= 1
        present[a2] = 1
        present[b2] = 1
        out[i], out[j] = a2, b2
    end
    return _rebuilt(g, out)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, m::SequenceShuffling) where {V,T}
    es = edges(g)
    snaps = _snapshots(es)
    if m.active
        times = [es[first(r)].t for r in snaps][Random.randperm(rng, length(snaps))]
    else
        a, b, Δ = _time_grid(g, m.resolution)
        Δ === nothing && throw(ArgumentError("SequenceShuffling(active = false) needs a time resolution"))
        K = fld(b - a, Δ) + 1
        K >= length(snaps) || throw(ArgumentError("more snapshots than times in the observation interval"))
        slots = Set{Int}()
        while length(slots) < length(snaps)
            push!(slots, rand(rng, 0:K-1))
        end
        times = [a + Δ * k for k in Random.shuffle!(rng, collect(slots))]
    end
    out = similar(es)
    for (r, t) in zip(snaps, times), i in r
        out[i] = TemporalEdge{V,T}(es[i].u, es[i].v, t, es[i].tt)
    end
    return _rebuilt(g, out)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, m::SnapshotShuffling) where {V,T}
    es = edges(g)
    n = num_nodes(g)
    out = similar(es)
    for r in _snapshots(es)
        _check_pairs(n, length(r), m.directed)
        seen = Set{Tuple{Int,Int}}()
        for i in r
            l = _random_link(rng, n, m.directed)
            while l in seen
                l = _random_link(rng, n, m.directed)
            end
            push!(seen, l)
            out[i] = TemporalEdge{V,T}(l[1], l[2], es[i].t, es[i].tt)
        end
    end
    return _rebuilt(g, out)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, m::DegreeSnapshotShuffling) where {V,T}
    es = edges(g)
    out = similar(es)
    for r in _snapshots(es)
        links = [_link_key(es[i].u, es[i].v, m.directed) for i in r]
        present = Dict{Tuple{Int,Int},Int}()
        for l in links
            present[l] = get(present, l, 0) + 1
        end
        _swap_links!(rng, links, m.swaps * length(r), m.directed, present)
        for (x, i) in enumerate(r)
            out[i] = _move(es[i], links[x], m.directed)
        end
    end
    return _rebuilt(g, out)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, ::IsomorphicSnapshotShuffling) where {V,T}
    es = edges(g)
    n = num_nodes(g)
    out = similar(es)
    image = Dict{Int,Int}()
    used = Set{Int}()
    for r in _snapshots(es)
        empty!(image)
        empty!(used)
        # a uniformly random injection of the nodes of the snapshot into 1:n
        node(u) = get!(image, Int(u)) do
            x = rand(rng, 1:n)
            while x in used
                x = rand(rng, 1:n)
            end
            push!(used, x)
            x
        end
        for i in r
            out[i] = TemporalEdge{V,T}(node(es[i].u), node(es[i].v), es[i].t, es[i].tt)
        end
    end
    return _rebuilt(g, out)
end

function _randomize(rng::AbstractRNG, g::OrderedEdgeList{V,T}, m::EventShuffling) where {V,T}
    es = edges(g)
    n = num_nodes(g)
    _check_pairs(n, isempty(es) ? 0 : 1, m.directed)
    grid = _time_grid(g, m.resolution)
    seen = Set{TemporalEdge{V,T}}()
    out = similar(es)
    for (i, e) in enumerate(es)
        for attempt in 1:1000
            l = _random_link(rng, n, m.directed)
            f = TemporalEdge{V,T}(l[1], l[2], _random_time(rng, grid, e.tt), e.tt)
            if !(f in seen)
                push!(seen, f)
                out[i] = f
                break
            end
            attempt == 1000 && throw(ArgumentError("too many events for the number of pairs and times"))
        end
    end
    return _rebuilt(g, out)
end

# ----------------------------------------------------------------------------------

"""
    randomize(g::OrderedEdgeList, models::ReferenceModel...; rng = Random.default_rng())

A random temporal network drawn from the microcanonical randomized reference model
`model` applied to `g` (Gauvin et al., SIAM Review 2022). With several models they
are applied in sequence: `randomize(g, A, B)` is the composition of `A` followed by
`B`. A link shuffling followed by a timeline shuffling, or a sequence shuffling
followed by a snapshot shuffling, is again a microcanonical model (Propositions V.5
and V.6 of the paper), e.g. `randomize(g, TopologyLinkShuffling(), TimestampShuffling())`
is ``P[\\mathcal{L}, p(\\mathbf{w}), \\mathbf{t}]``.

| model | canonical name | keeps |
|---|---|---|
| [`LinkShuffling`](@ref) | ``P[p_\\mathcal{L}(Θ)]`` | timelines (not their links), number of links |
| [`DegreeLinkShuffling`](@ref) | ``P[\\mathbf{k}, p_\\mathcal{L}(Θ)]`` | timelines, static degrees |
| [`TopologyLinkShuffling`](@ref) | ``P[\\mathcal{L}, p_\\mathcal{L}(Θ)]`` | timelines, static graph |
| [`WeightLinkShuffling`](@ref) | ``P[\\mathbf{w}, p_\\mathcal{L}(Θ)]`` | timelines, static graph, link weights |
| [`RandomTimes`](@ref) | ``P[\\mathbf{w}]`` | static graph, link weights |
| [`InterEventShuffling`](@ref) | ``P[π_\\mathcal{L}(Δτ), \\mathbf{t}_1]`` | static graph, inter-event times and first time of every link |
| [`TimelineShifting`](@ref) | ``P[\\mathrm{per}(Θ)]`` | timelines up to a periodic shift |
| [`TimestampShuffling`](@ref) | ``P[\\mathbf{w}, \\mathbf{t}]`` | static graph, link weights, times of the events |
| [`SequenceShuffling`](@ref) | ``P[p_\\mathcal{T}(Γ), χ_{\\mathbb{N}^+}(\\mathbf{A})]`` | snapshots (not their times) |
| [`SnapshotShuffling`](@ref) | ``P[\\mathbf{t}]`` | times of the events |
| [`DegreeSnapshotShuffling`](@ref) | ``P[\\mathbf{d}]`` | times, instantaneous degrees |
| [`IsomorphicSnapshotShuffling`](@ref) | ``P[\\mathrm{iso}(Γ)]`` | snapshots up to isomorphism |
| [`EventShuffling`](@ref) | ``P[p(τ)]`` | number and durations of the events |

Links are ordered pairs of nodes; for undirected networks store every contact once
and pass `directed = false` to the models (then use [`make_undirected`](@ref) on
the result if needed). The result has the same nodes and time interval as `g`.

# References

- L. Gauvin, M. Génois, M. Karsai, M. Kivelä, T. Takaguchi, E. Valdano, and C. L. Vestergaard. *Randomized reference models for temporal networks.* SIAM Review 64(4), 2022. [DOI](https://doi.org/10.1137/19M1242252), [arXiv](https://arxiv.org/abs/1806.04032)
"""
function randomize(g::OrderedEdgeList, models::ReferenceModel...; rng::AbstractRNG=Random.default_rng())
    for m in models
        g = _randomize(rng, g, m)
    end
    return g
end

"""
    reference_samples(f, g::OrderedEdgeList, models::ReferenceModel...; samples = 100,
                      rng = Random.default_rng())

The values `f(h)` for `samples` independent random networks `h = randomize(g,
models...)`, computed in parallel (each sample has its own random number generator
seeded from `rng`, so the result is reproducible). Comparing `f(g)` with these values
tells how much of `f` is explained by the features kept by the model, e.g. the
z-scores of motif counts:

```julia
null = reference_samples(h -> temporal_motif_counts(h, 3600), g, TimestampShuffling(); samples = 100)
z = (temporal_motif_counts(g, 3600) .- mean(null)) ./ std(null)
```
"""
function reference_samples(f, g::OrderedEdgeList, models::ReferenceModel...; samples::Integer=100,
                           rng::AbstractRNG=Random.default_rng())
    samples >= 1 || throw(ArgumentError("at least one sample is needed"))
    seeds = rand(rng, UInt64, samples)
    return tmap(s -> f(randomize(g, models...; rng=Random.Xoshiro(s))), seeds)
end
