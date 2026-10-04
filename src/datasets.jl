# Catalog of public temporal network datasets, downloaded on demand from their
# original sources (SNAP and SocioPatterns) and cached locally. The data are not
# redistributed: their licenses and the papers to cite are listed in the catalog.

"""
    TemporalDataset

Metadata of a dataset of the catalog returned by [`temporal_datasets`](@ref): `name`,
`source`, `url`, `description`, whether its edges are `directed`, its `license` and
the `citation` requested by the authors.
"""
struct TemporalDataset
    name::String
    source::String
    url::String
    file::String                 # data file after decompression
    description::String
    directed::Bool
    columns::NTuple{3,Int}       # columns of u, v and t
    delim::Union{Nothing,Char}   # nothing: whitespace
    header::Int                  # lines to skip
    license::String
    citation::String
end

function Base.show(io::IO, d::TemporalDataset)
    print(io, "TemporalDataset(\"", d.name, "\", ", d.source, ", ", d.directed ? "directed" : "undirected", ")")
end

function Base.show(io::IO, ::MIME"text/plain", d::TemporalDataset)
    println(io, d.name, " (", d.source, ", ", d.directed ? "directed" : "undirected", ")")
    println(io, "  ", d.description)
    println(io, "  url:      ", d.url)
    println(io, "  license:  ", d.license)
    print(io, "  cite:     ", d.citation)
end

const _SNAP = "https://snap.stanford.edu/data/"
const _SP = "http://www.sociopatterns.org/assets/data/"
const _SNAP_LICENSE = "see https://snap.stanford.edu/data/ (free for research use)"
const _SP_CC0 = "CC0 (Creative Commons Public Domain Dedication)"
const _SP_BYNCSA = "CC BY-NC-SA (Creative Commons Attribution-NonCommercial-ShareAlike)"

function _snap(name, description, citation; file=name * ".txt", columns=(1, 2, 3), delim=nothing)
    url = _SNAP * file * ".gz"
    return TemporalDataset(name, "SNAP", url, file, description, true, columns, delim, 0, _SNAP_LICENSE, citation)
end

function _sp(name, archive, file, description, license, citation; delim=nothing)
    return TemporalDataset(name, "SocioPatterns", _SP * archive, file, description, false, (2, 3, 1), delim, 0, license, citation)
end

const _CITE_PANZARASA = "P. Panzarasa, T. Opsahl, and K. M. Carley. Patterns and dynamics of users' behavior and interaction: network analysis of an online community. JASIST 60(5), 2009."
const _CITE_PBL = "A. Paranjape, A. R. Benson, and J. Leskovec. Motifs in temporal networks. WSDM, 2017."
const _CITE_KUMAR = "S. Kumar, F. Spezzano, V. S. Subrahmanian, and C. Faloutsos. Edge weight prediction in weighted signed networks. ICDM, 2016. And: S. Kumar, B. Hooi, D. Makhija, M. Kumar, C. Faloutsos, and V. S. Subrahmanian. REV2: fraudulent user prediction in rating platforms. WSDM, 2018."
const _CITE_WIKI = "J. Sun, J. Kunegis, and S. Staab. Predicting user roles in social networks using transfer learning with feature transformation. ICDM Workshops, 2016. And: J. Leskovec, D. Huttenlocher, and J. Kleinberg. Governance in social media: a case study of the Wikipedia promotion process. ICWSM, 2010."

function _catalog()
    d = TemporalDataset[]
    push!(d, _snap("CollegeMsg", "Private messages on an online social network at the University of California, Irvine (1,899 nodes, 59,835 edges, 193 days).", _CITE_PANZARASA))
    push!(d, _snap("email-Eu-core-temporal", "Emails between members of a European research institution (986 nodes, 332,334 edges, 803 days).", _CITE_PBL))
    for k in 1:4
        push!(d, _snap("email-Eu-core-temporal-Dept$k", "Emails inside department $k of a European research institution.", _CITE_PBL))
    end
    for site in ("mathoverflow", "askubuntu", "superuser", "stackoverflow"), (suffix, what) in (("", "all interactions"),
                 ("-a2q", "answers to questions"), ("-c2q", "comments to questions"), ("-c2a", "comments to answers"))
        push!(d, _snap("sx-$site$suffix", "Interactions on the Stack Exchange site $site, $what.", _CITE_PBL))
    end
    push!(d, _snap("wiki-talk-temporal", "Edits of user talk pages of the English Wikipedia (1,140,149 nodes, 7,833,140 edges).", _CITE_WIKI))
    push!(d, _snap("soc-sign-bitcoin-otc", "Who-trusts-whom network of the Bitcoin OTC platform (5,881 nodes, 35,592 rated edges; ratings are ignored).",
                   _CITE_KUMAR; file="soc-sign-bitcoinotc.csv", columns=(1, 2, 4), delim=','))
    push!(d, _snap("soc-sign-bitcoin-alpha", "Who-trusts-whom network of the Bitcoin Alpha platform (3,783 nodes, 24,186 rated edges; ratings are ignored).",
                   _CITE_KUMAR; file="soc-sign-bitcoinalpha.csv", columns=(1, 2, 4), delim=','))
    push!(d, _sp("sp-workplace-2013", "workplace_InVS_tij.dat.zip", "tij_InVS.dat",
                 "Face-to-face contacts in an office building in France, 2013 (92 nodes, 20 s resolution).", _SP_CC0,
                 "M. Génois, C. L. Vestergaard, J. Fournet, A. Panisson, I. Bonmarin, and A. Barrat. Data on face-to-face contacts in an office building suggest a low-cost vaccination strategy based on community linkers. Network Science 3(3), 2015."))
    push!(d, _sp("sp-workplace-2015", "workplace_InVS15_tij.dat.gz", "workplace_InVS15_tij.dat",
                 "Face-to-face contacts in an office building in France, 2015 (217 nodes, 20 s resolution).", _SP_CC0,
                 "M. Génois and A. Barrat. Can co-location be used as a proxy for face-to-face contacts? EPJ Data Science 7, 2018."))
    push!(d, _sp("sp-sfhh", "SFHH_tij.dat.gz", "SFHH_tij.dat",
                 "Face-to-face contacts at the SFHH conference in Nice, 2009 (403 nodes, 20 s resolution).", _SP_CC0,
                 "M. Génois and A. Barrat. Can co-location be used as a proxy for face-to-face contacts? EPJ Data Science 7, 2018."))
    push!(d, _sp("sp-primary-school", "primaryschool.csv.gz", "primaryschool.csv",
                 "Face-to-face contacts in a primary school in Lyon, 2009 (242 nodes, 20 s resolution).", _SP_BYNCSA,
                 "J. Stehlé et al. High-resolution measurements of face-to-face contact patterns in a primary school. PLoS ONE 6(8), 2011.";
                 delim='\t'))
    for (year, nodes) in (("2011", 126), ("2012", 180))
        push!(d, _sp("sp-high-school-$year", "highschool_$year.csv.gz", "highschool_$year.csv",
                     "Face-to-face contacts in a high school in Marseille, $year ($nodes nodes).", _SP_BYNCSA,
                     "J. Fournet and A. Barrat. Contact patterns among high school students. PLoS ONE 9(9), 2014.";
                     delim='\t'))
    end
    push!(d, _sp("sp-hospital", "hospital_lyon_contacts.dat.gz", "hospital_lyon_contacts.dat",
                 "Face-to-face contacts in a hospital ward in Lyon, 2010 (75 nodes, 20 s resolution).", _SP_BYNCSA,
                 "P. Vanhems et al. Estimating potential infection transmission routes in hospital wards using wearable proximity sensors. PLoS ONE 8(9), 2013.";
                 delim='\t'))
    push!(d, _sp("sp-hypertext-2009", "ht2009_contact_list.dat.gz", "ht2009_contact_list.dat",
                 "Face-to-face contacts at the ACM Hypertext 2009 conference (113 nodes, 20 s resolution).", _SP_BYNCSA,
                 "L. Isella et al. What's in a crowd? Analysis of face-to-face behavioral networks. Journal of Theoretical Biology 271(1), 2011."))
    return d
end

const _DATASETS = _catalog()

"""
    temporal_datasets()

The catalog of public temporal network datasets that [`load_dataset`](@ref) can
download: SNAP temporal networks (messages, e-mails, Stack Exchange interactions,
Wikipedia talk pages, Bitcoin trust networks) and SocioPatterns face-to-face contact
networks. Every entry is a [`TemporalDataset`](@ref) with its source, license and
the paper to cite.
"""
temporal_datasets() = copy(_DATASETS)

function _dataset(name::AbstractString)
    i = findfirst(d -> d.name == name, _DATASETS)
    i === nothing && throw(ArgumentError("unknown dataset \"$name\"; see temporal_datasets()"))
    return _DATASETS[i]
end

"""
    dataset_dir()

The directory where [`load_dataset`](@ref) caches the downloaded files: the
environment variable `TEMPORALGRAPHS_DATA` if set, otherwise a directory inside the
first Julia depot.
"""
dataset_dir() = get(ENV, "TEMPORALGRAPHS_DATA",
                    joinpath(first(DEPOT_PATH), "scratchspaces", "82aa5c3e-d57f-4b0f-ac00-57e9d2a8eea2", "datasets"))

# Path of the data file of d, downloading and decompressing it if needed.
function _dataset_file(d::TemporalDataset, dir::AbstractString)
    folder = joinpath(dir, d.name)
    path = joinpath(folder, d.file)
    isfile(path) && return path
    mkpath(folder)
    archive = joinpath(folder, basename(d.url))
    @info "Downloading $(d.name) from $(d.source) ($(d.url)). License: $(d.license). Please cite: $(d.citation)"
    Downloads.download(d.url, archive)
    exe = p7zip_jll.p7zip()
    if endswith(archive, ".gz")  # a single compressed file, whatever its stored name
        run(pipeline(`$exe e -so $archive`; stdout=path, stderr=devnull))
    else
        before = Set(readdir(folder))
        run(pipeline(`$exe x $archive -o$folder -y`; stdout=devnull, stderr=devnull))
        new = setdiff(readdir(folder), before)
        (!isfile(path) && length(new) == 1) && mv(joinpath(folder, only(new)), path)
    end
    isfile(path) || throw(ErrorException("$(d.file) not found in the downloaded archive $archive"))
    rm(archive; force=true)
    return path
end

# Columns u, v, t of a delimited text file (comment lines starting with % or # skipped).
function _read_columns(path::AbstractString, d::TemporalDataset)
    us, vs, ts = Int[], Int[], Int[]
    cu, cv, ct = d.columns
    k = max(cu, cv, ct)
    for (i, line) in enumerate(eachline(path))
        i <= d.header && continue
        s = strip(line)
        (isempty(s) || startswith(s, '%') || startswith(s, '#')) && continue
        fields = d.delim === nothing ? split(s) : split(s, d.delim)
        length(fields) >= k || throw(ArgumentError("line $i of $path has fewer than $k fields"))
        push!(us, parse(Int, fields[cu]))
        push!(vs, parse(Int, fields[cv]))
        x = fields[ct]
        push!(ts, something(tryparse(Int, x), round(Int, parse(Float64, x))))
    end
    return (u=us, v=vs, t=ts)
end

"""
    load_dataset(name; dir = dataset_dir(), directed = nothing, transition_time = 1,
                 node_type = Int32)

Load the temporal network `name` of the catalog [`temporal_datasets`](@ref) as an
[`OrderedEdgeList`](@ref), downloading it from its original source on first use (the
license and the paper to cite are printed then) and caching it in `dir`. Node ids
are renumbered `1:n`; [`original_ids`](@ref) gives the ids of the source. Every edge
gets the transition time `transition_time`. SocioPatterns contacts are undirected and
loaded in both directions unless `directed = true`.

```julia
g = load_dataset("CollegeMsg")
h = load_dataset("sp-workplace-2015")   # face-to-face contacts, 20 s resolution
```
"""
function load_dataset(name::AbstractString; dir::AbstractString=dataset_dir(), directed::Union{Nothing,Bool}=nothing,
                      transition_time::Real=1, node_type::Type{<:Integer}=Int32)
    d = _dataset(name)
    cols = _read_columns(_dataset_file(d, dir), d)
    return OrderedEdgeList(cols; u=:u, v=:v, t=:t, directed=something(directed, d.directed),
                           default_tt=transition_time, node_type=node_type)
end
