# Drawing temporal graphs as SVG images: the static graph with every edge labelled by
# its (t, tt), and one timeline per node. The images use grey strokes and dark text
# with a white halo, so they are readable on light and dark backgrounds.

"""
    TemporalGraphDrawing

An SVG drawing of a temporal graph, returned by [`draw_graph`](@ref) and
[`draw_timelines`](@ref). Notebooks, the VS Code plot pane and Documenter display it as
an image (`MIME"image/svg+xml"`); `write("file.svg", d)` saves it and `String(d)`
returns its SVG source.
"""
struct TemporalGraphDrawing
    svg::String
end

Base.String(d::TemporalGraphDrawing) = d.svg
Base.write(io::IO, d::TemporalGraphDrawing) = write(io, d.svg)
Base.show(io::IO, ::MIME"image/svg+xml", d::TemporalGraphDrawing) = print(io, d.svg)
Base.show(io::IO, d::TemporalGraphDrawing) = print(io, "TemporalGraphDrawing (", sizeof(d.svg), " bytes of SVG)")

const _NODE_COLORS = ["#4063D8", "#CB3C33", "#389826", "#9558B2"]  # the Julia colors
const _EDGE_COLOR = "#707070"
const _LINE_COLOR = "#8C8C8C"
const _FONT = "font-family=\"Lato, Helvetica, Arial, sans-serif\""
const _HALO = "$_FONT fill=\"#1A1A1A\" stroke=\"#FFFFFF\" stroke-width=\"4\" stroke-linejoin=\"round\" paint-order=\"stroke\""

_num(x::Real) = (y = round(Float64(x); digits=1); isinteger(y) ? string(Int(y)) : string(y))
_label(x::Real) = isinteger(x) ? string(Int(x)) : string(x)
_escape(s::AbstractString) = replace(s, "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")

function _check_nodes(n, node_labels, node_colors)
    length(node_labels) == n || throw(ArgumentError("node_labels must have length $n"))
    isempty(node_colors) && throw(ArgumentError("node_colors must not be empty"))
    return nothing
end

function _svg_start(io, w, h, title)
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 ", _num(w), " ", _num(h),
            "\" width=\"", _num(w), "\" height=\"", _num(h), "\">")
    println(io, "  <title>", _escape(title), "</title>")
    println(io, "  <defs>\n    <marker id=\"head\" viewBox=\"0 0 10 10\" refX=\"10\" refY=\"5\" ",
            "markerWidth=\"4.4\" markerHeight=\"4.4\" orient=\"auto\">\n      <path d=\"M0,0 L10,5 L0,10 Z\" ",
            "fill=\"", _EDGE_COLOR, "\"/>\n    </marker>\n  </defs>")
end

function _node(io, x, y, r, label, color)
    println(io, "  <circle cx=\"", _num(x), "\" cy=\"", _num(y), "\" r=\"", _num(r), "\" fill=\"", color, "\"/>")
    println(io, "  <text x=\"", _num(x), "\" y=\"", _num(y), "\" text-anchor=\"middle\" dominant-baseline=\"central\" ",
            "font-size=\"", _num(0.9r), "\" font-weight=\"bold\" fill=\"#FFFFFF\" ", _FONT, ">", _escape(string(label)), "</text>")
end

# Text centred at (x, y), rotated by `angle` degrees.
function _text(io, x, y, s; angle=0.0, anchor="middle", size=15)
    rot = abs(angle) > 0.01 ? " transform=\"rotate($(_num(angle)) $(_num(x)) $(_num(y)))\"" : ""
    println(io, "  <text x=\"", _num(x), "\" y=\"", _num(y), "\" text-anchor=\"", anchor, "\" dominant-baseline=\"central\" ",
            "font-size=\"", size, "\" ", _HALO, rot, ">", _escape(s), "</text>")
end

# Default positions of draw_graph: the nodes on a circle, clockwise, the first one on
# the top left (with 4 nodes they are the corners of a square). y points up.
function _circle_positions(n)
    n == 1 && return [(0.0, 0.0)]
    return [(cos(θ), sin(θ)) for θ in (π / 2 + π / n - 2π * (i - 1) / n for i in 1:n)]
end

"""
    draw_graph(g; positions = nothing, node_labels = 1:n, node_colors = Julia colors,
               edge_labels = true, width = 420, node_radius = 20)

Draws the temporal graph `g` as its static graph: one arrow from `u` to `v` for every
temporal edge `(u, v, t, tt)`, labelled with `(t, tt)` if `edge_labels`. The edges
between the same two nodes are drawn as separate arcs and self-loops as loops.
Returns a [`TemporalGraphDrawing`](@ref); see [`draw_timelines`](@ref) for a drawing
that shows the times.

The nodes are placed on a circle, unless `positions[i] = (x, y)` gives the position of
node `i` (with `y` pointing up, in any unit: the drawing is scaled to `width` pixels).
`node_colors` is cycled over the nodes. Meant for small graphs: every edge is drawn.
"""
function draw_graph(g::OrderedEdgeList; positions=nothing, node_labels=1:num_nodes(g),
                    node_colors=_NODE_COLORS, edge_labels::Bool=true, width::Real=420, node_radius::Real=20)
    n = num_nodes(g)
    _check_nodes(n, node_labels, node_colors)
    ps = positions === nothing ? _circle_positions(n) : [(Float64(p[1]), Float64(p[2])) for p in positions]
    length(ps) == n || throw(ArgumentError("positions must have length $n"))
    r = Float64(node_radius)
    margin = r + (edge_labels ? 50 : 30)
    # fit the positions into the drawing, flipping y
    xs, ys = first.(ps), last.(ps)
    lo_x, lo_y = n == 0 ? (0.0, 0.0) : (minimum(xs), minimum(ys))
    dx, dy = n == 0 ? (0.0, 0.0) : (maximum(xs) - lo_x, maximum(ys) - lo_y)
    s = max(dx, dy) > 0 ? (width - 2margin) / max(dx, dy) : 0.0
    w, h = Float64(width), dy * s + 2margin
    pos = [(margin + (x - lo_x) * s + (dx < dy ? (dy - dx) * s / 2 : 0.0), margin + (lo_y + dy - y) * s) for (x, y) in ps]
    center = (w / 2, h / 2)

    io = IOBuffer()
    _svg_start(io, w, h, "Temporal graph")
    labels = Tuple{Float64,Float64,String,Float64,String}[]  # x, y, text, angle, anchor
    println(io, "  <g stroke=\"", _EDGE_COLOR, "\" stroke-width=\"2.5\" fill=\"none\" marker-end=\"url(#head)\">")
    # the edges between the same two nodes, in chronological order
    groups = Dict{Tuple{Int,Int},Vector{Int}}()
    for (k, e) in enumerate(edges(g))
        push!(get!(groups, minmax(Int(e.u), Int(e.v)), Int[]), k)
    end
    for k in sort!(collect(keys(groups)))
        ids = groups[k]
        a, b = k
        m = length(ids)
        if a == b  # self-loops: loops of growing size away from the centre
            (x, y) = pos[a]
            φ = hypot(x - center[1], y - center[2]) < 1 ? -π / 2 : atan(y - center[2], x - center[1])
            for (i, id) in enumerate(ids)
                ℓ = r * (2 + 1.2i)
                p1 = (x + r * cos(φ - 0.45), y + r * sin(φ - 0.45))
                p2 = (x + r * cos(φ + 0.45), y + r * sin(φ + 0.45))
                c1 = (x + ℓ * cos(φ - 0.6), y + ℓ * sin(φ - 0.6))
                c2 = (x + ℓ * cos(φ + 0.6), y + ℓ * sin(φ + 0.6))
                println(io, "    <path d=\"M", _num(p1[1]), ",", _num(p1[2]), " C", _num(c1[1]), ",", _num(c1[2]), " ",
                        _num(c2[1]), ",", _num(c2[2]), " ", _num(p2[1]), ",", _num(p2[2]), "\"/>")
                e = edges(g)[id]
                d = 0.75ℓ + 14
                cosφ = cos(φ)
                push!(labels, (x + d * cosφ, y + d * sin(φ), "($(_label(e.t)), $(_label(e.tt)))", 0.0,
                               abs(cosφ) < 0.4 ? "middle" : cosφ > 0 ? "start" : "end"))
            end
            continue
        end
        (xa, ya), (xb, yb) = pos[a], pos[b]
        L = hypot(xb - xa, yb - ya)
        L > 0 || continue
        ux, uy = (xb - xa) / L, (yb - ya) / L
        nx, ny = -uy, ux                      # normal of the direction a → b
        mx, my = (xa + xb) / 2, (ya + yb) / 2
        # side of the labels of a straight edge: away from the centre, else up
        out = nx * (mx - center[1]) + ny * (my - center[2])
        side0 = abs(out) > 1 ? sign(out) : (ny <= 0 ? 1.0 : -1.0)
        for (i, id) in enumerate(ids)
            e = edges(g)[id]
            off = (i - (m + 1) / 2) * 34.0      # distance of the arc from the segment
            cx, cy = mx + 2off * nx, my + 2off * ny   # control point of the arc
            (xu, yu), (xv, yv) = pos[e.u], pos[e.v]
            du, dv = hypot(cx - xu, cy - yu), hypot(cx - xv, cy - yv)
            p1 = (xu + r * (cx - xu) / du, yu + r * (cy - yu) / du)
            p2 = (xv + r * (cx - xv) / dv, yv + r * (cy - yv) / dv)
            println(io, "    <path d=\"M", _num(p1[1]), ",", _num(p1[2]), " Q", _num(cx), ",", _num(cy), " ",
                    _num(p2[1]), ",", _num(p2[2]), "\"/>")
            side = off == 0 ? side0 : sign(off)
            # the label is beyond the middle of the arc
            ax, ay = (p1[1] + p2[1]) / 4 + cx / 2, (p1[2] + p2[2]) / 4 + cy / 2
            lx, ly = ax + side * 13 * nx, ay + side * 13 * ny
            angle = rad2deg(atan(uy, ux))
            angle > 90 && (angle -= 180)
            angle < -90 && (angle += 180)
            text = "($(_label(e.t)), $(_label(e.tt)))"
            if abs(angle) > 60  # nearly vertical: horizontal text beside the arrow
                push!(labels, (lx + side * nx * 8, ly, text, 0.0, side * nx > 0 ? "start" : "end"))
            else
                push!(labels, (lx, ly, text, angle, "middle"))
            end
        end
    end
    println(io, "  </g>")
    if edge_labels
        for (x, y, text, angle, anchor) in labels
            _text(io, x, y, text; angle, anchor)
        end
    end
    for i in 1:n
        _node(io, pos[i]..., r, node_labels[i], node_colors[mod1(i, length(node_colors))])
    end
    println(io, "</svg>")
    return TemporalGraphDrawing(String(take!(io)))
end

# Tick marks of the time axis: every time stamp if there are few integer ones, else
# about 8 round values.
function _ticks(a::Real, b::Real)
    a == b && return [a]
    isinteger(a) && isinteger(b) && b - a <= 20 && return collect(a:b)
    raw = (b - a) / 8
    mag = 10.0^floor(log10(raw))
    f = raw / mag
    step = mag * (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10)
    return collect(ceil(a / step) * step:step:b)
end

"""
    draw_timelines(g, ti = time_interval(g); order = 1:n, node_labels = 1:n,
                   node_colors = Julia colors, width = 600, row_height = 70)

Draws the temporal graph `g` with one horizontal timeline per node and time on the
horizontal axis: every temporal edge `(u, v, t, tt)` is an arrow from the timeline of
`u` at time `t` to the timeline of `v` at time `t + tt`, and the dots mark the
departures and arrivals. Returns a [`TemporalGraphDrawing`](@ref); see
[`draw_graph`](@ref) for the static graph.

`ti = (a, b)` is the time range of the axis; only the edges inside it are drawn.
`order` lists the nodes from top to bottom (a subset of the nodes draws only the edges
between them); changing it can avoid arrows that cross the dots of other nodes.
`node_colors` is cycled over the nodes. Meant for small graphs: every edge is drawn.
"""
function draw_timelines(g::OrderedEdgeList{V,T}, ti=time_interval(g); order=1:num_nodes(g),
                        node_labels=1:num_nodes(g), node_colors=_NODE_COLORS,
                        width::Real=600, row_height::Real=70) where {V,T}
    n = num_nodes(g)
    _check_nodes(n, node_labels, node_colors)
    all(i -> 1 <= i <= n, order) && allunique(order) ||
        throw(ArgumentError("order must list distinct nodes of 1:$n"))
    a, b = _interval(T, ti)
    row = zeros(Int, n)
    for (k, i) in enumerate(order)
        row[i] = k
    end
    top, left, right = 30.0, 75.0, 50.0
    w = Float64(width)
    yof(i) = top + (row[i] - 1) * row_height
    xof(t) = b > a ? left + (t - a) / (b - a) * (w - left - right) : (left + w - right) / 2
    bottom = top + (max(length(order), 1) - 1) * row_height + 20
    h = bottom + 50
    color(i) = node_colors[mod1(i, length(node_colors))]

    io = IOBuffer()
    _svg_start(io, w, h, "Timelines of a temporal graph")
    ticks = _ticks(a, b)
    println(io, "  <g stroke=\"", _LINE_COLOR, "\" stroke-width=\"1\" stroke-opacity=\"0.45\">")
    for t in ticks
        println(io, "    <line x1=\"", _num(xof(t)), "\" y1=\"", _num(top - 20), "\" x2=\"", _num(xof(t)), "\" y2=\"", _num(bottom), "\"/>")
    end
    println(io, "  </g>")
    println(io, "  <path d=\"M", _num(left - 30), ",", _num(bottom + 20), " H", _num(w - 30), "\" stroke=\"", _EDGE_COLOR,
            "\" stroke-width=\"2\" fill=\"none\" marker-end=\"url(#head)\"/>")
    for t in ticks
        _text(io, xof(t), bottom + 40, _label(t); size=14)
    end
    _text(io, w - 22, bottom + 20, "t"; anchor="start", size=15)
    println(io, "  <g stroke=\"", _LINE_COLOR, "\" stroke-width=\"2\" stroke-dasharray=\"4 4\" stroke-linecap=\"round\">")
    for i in order
        println(io, "    <line x1=\"", _num(left - 30), "\" y1=\"", _num(yof(i)), "\" x2=\"", _num(w - right + 20),
                "\" y2=\"", _num(yof(i)), "\"/>")
    end
    println(io, "  </g>")
    # the edges inside ti between the drawn nodes, and their departures and arrivals
    dot = 5.0
    events = Set{Tuple{T,Int}}()
    println(io, "  <g stroke=\"", _EDGE_COLOR, "\" stroke-width=\"2.5\" fill=\"none\" marker-end=\"url(#head)\">")
    for e in edges(g)
        u, v = Int(e.u), Int(e.v)
        (row[u] > 0 && row[v] > 0 && a <= e.t && e.t + e.tt <= b) || continue
        push!(events, (e.t, u), (e.t + e.tt, v))
        x1, y1, x2, y2 = xof(e.t), yof(u), xof(e.t + e.tt), yof(v)
        L = hypot(x2 - x1, y2 - y1)
        L > 2dot || continue
        if u == v  # a self-loop: an arc above the timeline
            cy = y1 - 0.6row_height
            println(io, "    <path d=\"M", _num(x1), ",", _num(y1 - dot), " Q", _num((x1 + x2) / 2), ",", _num(cy), " ",
                    _num(x2), ",", _num(y2 - dot), "\"/>")
        else
            ux, uy = (x2 - x1) / L, (y2 - y1) / L
            println(io, "    <path d=\"M", _num(x1 + dot * ux), ",", _num(y1 + dot * uy), " L", _num(x2 - dot * ux), ",",
                    _num(y2 - dot * uy), "\"/>")
        end
    end
    println(io, "  </g>")
    for (t, i) in sort!(collect(events))
        println(io, "  <circle cx=\"", _num(xof(t)), "\" cy=\"", _num(yof(i)), "\" r=\"", _num(dot), "\" fill=\"", color(i), "\"/>")
    end
    for i in order
        _node(io, 25.0, yof(i), 16.0, node_labels[i], color(i))
    end
    println(io, "</svg>")
    return TemporalGraphDrawing(String(take!(io)))
end

draw_graph(g; kwargs...) = draw_graph(to_ordered_edge_list(g); kwargs...)
draw_timelines(g, ti...; kwargs...) = draw_timelines(to_ordered_edge_list(g), ti...; kwargs...)
