# SVG drawings of temporal graphs.

@testset "drawing" begin
    g = paper_file()
    arrows(svg) = length(collect(eachmatch(r"<path d=\"M", svg))) - 1  # minus the arrow head marker
    for d in (draw_graph(g), draw_timelines(g), draw_graph(to_incident_lists(g)), draw_timelines(to_incident_lists(g)))
        svg = String(d)
        @test startswith(svg, "<svg") && endswith(svg, "</svg>\n")
        @test arrows(svg) == num_edges(g) + (occursin("Timelines", svg) ? 1 : 0)  # + the time axis
        @test repr(MIME"image/svg+xml"(), d) == svg
        @test occursin("bytes of SVG", repr(d))
        path = tempname() * ".svg"
        write(path, d)
        @test read(path, String) == svg
    end
    # one label (t, tt) per edge, the nodes are labelled with node_labels
    svg = String(draw_graph(g; node_labels=["a", "b", "c", "<d>"]))
    @test all(e -> occursin("($(e.t), $(e.tt))", svg), edges(g))
    @test occursin("&lt;d&gt;", svg) && !occursin("<d>", svg)
    @test !occursin("(1, 5)", String(draw_graph(g; edge_labels=false)))
    @test occursin("viewBox=\"0 0 300 ", String(draw_graph(g; width=300)))
    # self-loops, opposite edges, floating point times, a single node, no nodes
    h = OrderedEdgeList(3, [(1, 1, 0.5, 1.5), (1, 2, 1.0, 0.0), (2, 1, 1.0, 2.0), (2, 3, 4.0, 1.0)])
    @test arrows(String(draw_graph(h; positions=[(0, 0), (1, 0), (2, 1)]))) == 4
    @test occursin("(0.5, 1.5)", String(draw_graph(h)))
    @test arrows(String(draw_timelines(h))) == 5
    @test arrows(String(draw_timelines(h, (0.5, 3.0)))) == 4  # (2 3 4 1) is outside
    @test arrows(String(draw_timelines(h; order=[2, 3]))) == 2
    @test occursin("<circle", String(draw_graph(OrderedEdgeList(1, TemporalEdge{Int32,Int}[]))))
    @test !occursin("<circle", String(draw_graph(OrderedEdgeList(0, TemporalEdge{Int32,Int}[]))))
    @test_throws ArgumentError draw_graph(g; node_labels=1:3)
    @test_throws ArgumentError draw_graph(g; positions=[(0, 0)])
    @test_throws ArgumentError draw_graph(g; node_colors=String[])
    @test_throws ArgumentError draw_timelines(g; order=[1, 1])
    @test_throws ArgumentError draw_timelines(g; order=[5])
end
