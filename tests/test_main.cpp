#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "graph.hpp"
#include "vtk_data.hpp"
#include "random_walk.hpp"

// Defined in graph.cpp (not in a header).
Vec2 intersection_helper(const Graph& g, std::pair<int, int> l1, std::pair<int, int> l2, bool& valid);

namespace {

// ---------- helpers for small hand-built datasets ----------

struct TempFile {
    std::filesystem::path path;
    explicit TempFile(const std::string& contents) {
        path = std::filesystem::temp_directory_path() / "braindiff_graph_test.vtk";
        std::ofstream(path) << contents;
    }
    ~TempFile() { std::filesystem::remove(path); }
    std::string str() const { return path.string(); }
};

using Pts = std::vector<std::array<double, 3>>;
using Lns = std::vector<std::vector<int>>;

std::string vtk_text(const Pts& pts, const Lns& lines) {
    std::ostringstream s;
    s << "# vtk DataFile Version 3.0\ntest\nASCII\nDATASET POLYDATA\n";
    s << "POINTS " << pts.size() << " float\n";
    for (const auto& p : pts) s << p[0] << " " << p[1] << " " << p[2] << "\n";
    size_t total = 0;
    for (const auto& l : lines) total += l.size() + 1;
    s << "LINES " << lines.size() << " " << total << "\n";
    for (const auto& l : lines) {
        s << l.size();
        for (int i : l) s << " " << i;
        s << "\n";
    }
    return s.str();
}

Graph graph_from(const Pts& pts, const Lns& lines) {
    TempFile f(vtk_text(pts, lines));
    VTKdata v(f.str());
    return Graph(v);
}

size_t edge_count(Graph& g) {
    size_t sum = 0;
    for (size_t i = 0; i < g.length(); ++i) sum += g.neighbors(int(i)).size();
    return sum / 2;
}

bool is_symmetric(Graph& g) {
    for (int i = 0; i < int(g.length()); ++i)
        for (int j : g.neighbors(i))
            if (!g.are_neighbors(j, i)) return false;
    return true;
}

}  // namespace

// ---------- existing tests ----------

TEST_CASE("VTKdata matches expected counts: data/001k.vtk") {
    const std::string path = "data/001k.vtk";

    VTKdata v(path);
    CHECK(v.points().size() == 163991);
    CHECK(v.lines().size() == 1000);

    for (const auto& line : v.lines())
        for (int idx : line)
            CHECK((idx >= 0 && size_t(idx) < v.points().size()));
}

TEST_CASE("VTKdata after filter pass matches expected values") {
    const std::string path = "data/001k.vtk";

    VTKdata unfiltered_v(path);
    VTKdata v = unfiltered_v.filter(zfilter<-10.0, 10.0>, 100);

    CHECK(v.points().size() == 15448);
    CHECK(v.lines().size() == 99);

    for (const auto& line : v.lines())
        for (int idx : line)
            CHECK((idx >= 0 && size_t(idx) < v.points().size()));
}

TEST_CASE("Verifying expected graph behavior without intersections and filtering") {
    VTKdata v("data/001k.vtk");
    Graph g(v);
    for (int i = 0; i < int(g.length()); i++) {
        CHECK((g.neighbors(i).size() == 1 || g.neighbors(i).size() == 2));
    }
}

TEST_CASE("Verifying expected graph behavior after filtering") {
    VTKdata v1("data/001k.vtk");
    VTKdata v2 = v1.filter(zfilter<-10.0, 10.0>, 100);
    Graph g(v2);
    for (int i = 0; i < int(g.length()); i++) {
        CHECK((g.neighbors(i).size() == 1 || g.neighbors(i).size() == 2));
    }

    for (int i = 0; i < int(g.length()); i++) {
        for (int j : g.neighbors(i)) {
            CHECK(g.are_neighbors(i, j) == g.are_neighbors(j, i));
        }
    }
}


TEST_CASE("Add intersections check") {
    VTKdata v("data/001k.vtk");
    Graph g(v.filter(zfilter<-10.0, 10.0>, 100));
    g.add_intersections();

    CHECK(g.length() == 17966);
}

// ---------- new tests: small hand-built graphs ----------

TEST_CASE("convenience data on a simple chain") {
    // 3 collinear points, edges of length 1 and 2
    Graph g = graph_from({{{0, 0, 0}}, {{1, 0, 0}}, {{3, 0, 0}}}, {{0, 1, 2}});
    CHECK(g.length() == 3);
    CHECK(g.min_x == doctest::Approx(0.0));
    CHECK(g.max_x == doctest::Approx(3.0));
    CHECK(g.min_y == doctest::Approx(0.0));
    CHECK(g.max_y == doctest::Approx(0.0));
    CHECK(g.min_norm == doctest::Approx(1.0));
    CHECK(g.max_norm == doctest::Approx(2.0));
    CHECK(g.avg_norm == doctest::Approx(1.5));
}

TEST_CASE("adjacency of a chain is symmetric with correct degrees") {
    Graph g = graph_from({{{0, 0, 0}}, {{1, 0, 0}}, {{2, 0, 0}}, {{3, 0, 0}}}, {{0, 1, 2, 3}});
    CHECK(g.neighbors(0).size() == 1);
    CHECK(g.neighbors(1).size() == 2);
    CHECK(g.neighbors(2).size() == 2);
    CHECK(g.neighbors(3).size() == 1);
    CHECK(g.are_neighbors(0, 1));
    CHECK(g.are_neighbors(1, 0));
    CHECK_FALSE(g.are_neighbors(0, 2));
    CHECK(is_symmetric(g));
}

TEST_CASE("duplicate edges from two lines are deduplicated") {
    // both lines traverse segment 0-1; the second continues on to 2
    Graph g = graph_from({{{0, 0, 0}}, {{1, 0, 0}}, {{2, 0, 0}}}, {{0, 1}, {0, 1, 2}});
    CHECK(g.neighbors(0).size() == 1);   // just 1, not listed twice
    CHECK(g.neighbors(1).size() == 2);   // 0 and 2
    CHECK(g.neighbors(2).size() == 1);
    CHECK(edge_count(g) == 2);
    CHECK(is_symmetric(g));
}

TEST_CASE("intersection_helper: crossing, parallel, short, shared endpoint") {
    Graph g = graph_from(
        {{{0, 0, 0}}, {{2, 2, 0}},      // 0,1  diagonal
         {{0, 2, 0}}, {{2, 0, 0}},      // 2,3  anti-diagonal (crosses 0-1)
         {{0, 1, 0}}, {{2, 1, 0}},      // 4,5  horizontal through the middle
         {{5, 5, 0}}, {{6, 6, 0}}},     // 6,7  far away
        {{0, 1}, {2, 3}, {4, 5}, {6, 7}});

    bool ok = false;
    Vec2 p = intersection_helper(g, {0, 1}, {2, 3}, ok);
    CHECK(ok);
    CHECK(p.getX() == doctest::Approx(1.0));
    CHECK(p.getY() == doctest::Approx(1.0));

    // symmetric in its arguments
    Vec2 q = intersection_helper(g, {2, 3}, {0, 1}, ok);
    CHECK(ok);
    CHECK(q.getX() == doctest::Approx(1.0));
    CHECK(q.getY() == doctest::Approx(1.0));

    // would cross if extended, but stops short
    intersection_helper(g, {0, 1}, {6, 7}, ok);
    CHECK_FALSE(ok);

    // sharing an endpoint is not a new crossing
    intersection_helper(g, {0, 1}, {1, 0}, ok);
    CHECK_FALSE(ok);
}

TEST_CASE("intersection_helper: parallel segments never cross") {
    Graph g = graph_from({{{0, 0, 0}}, {{2, 0, 0}}, {{0, 1, 0}}, {{2, 1, 0}}}, {{0, 1}, {2, 3}});
    bool ok = true;
    intersection_helper(g, {0, 1}, {2, 3}, ok);
    CHECK_FALSE(ok);
}

TEST_CASE("add_intersections: X shape gets one degree-4 vertex at the center") {
    Graph g = graph_from({{{0, 0, 0}}, {{2, 2, 0}}, {{0, 2, 0}}, {{2, 0, 0}}}, {{0, 1}, {2, 3}});
    g.add_intersections();

    REQUIRE(g.length() == 5);
    CHECK(g[4].getX() == doctest::Approx(1.0));
    CHECK(g[4].getY() == doctest::Approx(1.0));
    CHECK(g.neighbors(4).size() == 4);

    // the original crossing edges are gone, replaced by edges to the center
    CHECK_FALSE(g.are_neighbors(0, 1));
    CHECK_FALSE(g.are_neighbors(2, 3));
    for (int corner = 0; corner < 4; ++corner) {
        CHECK(g.are_neighbors(corner, 4));
        CHECK(g.neighbors(corner).size() == 1);
    }
    CHECK(is_symmetric(g));
}

TEST_CASE("add_intersections: parallel lines are untouched") {
    Graph g = graph_from({{{0, 0, 0}}, {{2, 0, 0}}, {{0, 1, 0}}, {{2, 1, 0}}}, {{0, 1}, {2, 3}});
    g.add_intersections();
    CHECK(g.length() == 4);
    CHECK(edge_count(g) == 2);
}

TEST_CASE("add_intersections: lines meeting at a shared vertex are untouched") {
    // a V shape: one line, so its edges share vertex 1
    Graph g = graph_from({{{0, 0, 0}}, {{1, 1, 0}}, {{2, 0, 0}}}, {{0, 1, 2}});
    g.add_intersections();
    CHECK(g.length() == 3);
    CHECK(edge_count(g) == 2);
}

TEST_CASE("add_intersections: one long edge crossed twice") {
    // horizontal edge 0-1 is crossed by two vertical edges.
    // The second crossing is only findable after the first split.
    Graph g = graph_from(
        {{{0, 0, 0}}, {{10, 0, 0}},     // 0,1
         {{3, -1, 0}}, {{3, 1, 0}},     // 2,3
         {{7, -1, 0}}, {{7, 1, 0}}},    // 4,5
        {{0, 1}, {2, 3}, {4, 5}});
    g.add_intersections();

    REQUIRE(g.length() == 8);
    CHECK(edge_count(g) == 3 + 2 * 2);   // each crossing replaces 2 edges with 4

    // the two new vertices sit at (3,0) and (7,0) and have degree 4
    std::vector<double> xs;
    for (int i = 6; i < 8; ++i) {
        CHECK(g.neighbors(i).size() == 4);
        CHECK(g[i].getY() == doctest::Approx(0.0));
        xs.push_back(g[i].getX());
    }
    std::sort(xs.begin(), xs.end());
    CHECK(xs[0] == doctest::Approx(3.0));
    CHECK(xs[1] == doctest::Approx(7.0));

    // the original long edge is now a chain 0 - x(3) - x(7) - 1
    CHECK_FALSE(g.are_neighbors(0, 1));
    CHECK(is_symmetric(g));
}

TEST_CASE("add_intersections is idempotent") {
    Graph g = graph_from({{{0, 0, 0}}, {{2, 2, 0}}, {{0, 2, 0}}, {{2, 0, 0}}}, {{0, 1}, {2, 3}});
    g.add_intersections();
    size_t after_first = g.length();
    g.add_intersections();
    CHECK(g.length() == after_first);
}

// ---------- new tests: real data invariants ----------

TEST_CASE("add_intersections on real data preserves graph invariants") {
    VTKdata v("data/001k.vtk");
    Graph g(v.filter(zfilter<-10.0, 10.0>, 100));

    const size_t before = g.length();
    const size_t edges_before = edge_count(g);
    std::vector<size_t> degree_before(before);
    for (size_t i = 0; i < before; ++i) degree_before[i] = g.neighbors(int(i)).size();

    g.add_intersections();

    const size_t added = g.length() - before;
    CHECK(added == 17966 - before);
    CHECK(added > 0);

    // every crossing replaces 2 edges with 4: net +2 edges per new vertex
    CHECK(edge_count(g) == edges_before + 2 * added);

    // original vertices keep their degree (a neighbor is swapped, never added)
    for (size_t i = 0; i < before; ++i)
        CHECK(g.neighbors(int(i)).size() == degree_before[i]);

    // new vertices are all degree 4
    for (size_t i = before; i < g.length(); ++i)
        CHECK(g.neighbors(int(i)).size() == 4);

    CHECK(is_symmetric(g));

    // no self loops, no duplicate neighbors
    for (int i = 0; i < int(g.length()); ++i) {
        auto nb = g.neighbors(i);
        CHECK_FALSE(std::ranges::contains(nb, i));
        std::sort(nb.begin(), nb.end());
        CHECK(std::adjacent_find(nb.begin(), nb.end()) == nb.end());
    }

    // idempotent on real data too
    const size_t after = g.length();
    g.add_intersections();
    CHECK(g.length() == after);
}

TEST_CASE("no proper crossings remain after add_intersections (sampled brute force)") {
    VTKdata v("data/001k.vtk");
    Graph g(v.filter(zfilter<-10.0, 10.0>, 100));
    g.add_intersections();

    std::vector<std::pair<int, int>> edges;
    for (int a = 0; a < int(g.length()); ++a)
        for (int b : g.neighbors(a))
            if (a < b) edges.push_back({a, b});

    // test every 20th edge against all edges
    size_t remaining = 0;
    for (size_t i = 0; i < edges.size(); i += 20) {
        for (size_t j = 0; j < edges.size(); ++j) {
            if (i == j) continue;
            bool ok = false;
            intersection_helper(g, edges[i], edges[j], ok);
            if (ok) ++remaining;
        }
    }
    CHECK(remaining == 0);
}

TEST_CASE("PlanarGraphView puts every vertex in exactly one cell") {
    VTKdata v("data/001k.vtk");
    Graph g(v.filter(zfilter<-10.0, 10.0>, 100));
    g.refresh_convenience_data();
    PlanarGraphView view(g, g.max_norm * 2);

    std::vector<int> seen(g.length(), 0);
    for (size_t x = 0; x < view.w; ++x)
        for (size_t y = 0; y < view.h; ++y)
            for (int id : view[{x, y}]) {
                REQUIRE((id >= 0 && size_t(id) < g.length()));
                ++seen[id];
            }
    for (int count : seen) CHECK(count == 1);
}

TEST_CASE("PlanarGraphView to_cell puts each vertex in its own cell") {
    VTKdata v("data/001k.vtk");
    Graph g(v.filter(zfilter<-10.0, 10.0>, 100));
    g.refresh_convenience_data();
    PlanarGraphView view(g, g.max_norm * 2);

    for (int i = 0; i < int(g.length()); ++i) {
        auto cell = view.to_cell(g[i]);
        CHECK(cell.first < view.w);
        CHECK(cell.second < view.h);
        CHECK(std::ranges::contains(view[cell], i));
    }
}

TEST_CASE("timing: full pipeline on data/010k.vtk"*doctest::skip()) {
    using clock = std::chrono::steady_clock;
    auto seconds = [](clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    constexpr bool APPLY_FILTER = false;

    auto t0 = clock::now();
    VTKdata raw("data/010k.vtk");
    auto t1 = clock::now();

    VTKdata data = APPLY_FILTER ? raw.filter(zfilter<-35.0,45.0>, 100) : raw;
    auto t2 = clock::now();

    Graph g(data);
    auto t3 = clock::now();

    const size_t before = g.length();
    g.add_intersections();
    auto t4 = clock::now();

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "010k pipeline (" << (APPLY_FILTER ? "filtered" : "unfiltered") << ")\n"
              << "  file points / lines : " << raw.points().size() << " / " << raw.lines().size() << "\n"
              << "  graph vertices      : " << before << "\n"
              << "  intersections added : " << (g.length() - before) << "\n"
              << "  load                : " << seconds(t0, t1) << " s\n"
              << "  filter              : " << seconds(t1, t2) << " s\n"
              << "  graph build         : " << seconds(t2, t3) << " s\n"
              << "  add_intersections   : " << seconds(t3, t4) << " s\n"
              << "  total               : " << seconds(t0, t4) << " s\n";

    CHECK(g.length() >= before);
}

namespace {

// X-shaped graph: after add_intersections, vertex 4 is the center (degree 4).
Graph make_x_graph() {
    Graph g = graph_from({{{0, 0, 0}}, {{2, 2, 0}}, {{0, 2, 0}}, {{2, 0, 0}}}, {{0, 1}, {2, 3}});
    g.add_intersections();
    return g;
}

double total_mass(const ParticleAggregate& agg, const PlanarGraphView& view) {
    double sum = 0;
    for (size_t x = 0; x < view.w; ++x)
        for (size_t y = 0; y < view.h; ++y) sum += agg[int(x), int(y)];
    return sum;
}

size_t nonzero_cells(const ParticleAggregate& agg, const PlanarGraphView& view) {
    size_t n = 0;
    for (size_t x = 0; x < view.w; ++x)
        for (size_t y = 0; y < view.h; ++y) n += agg[int(x), int(y)] > 0;
    return n;
}

}  // namespace

// ---------- Particle ----------

TEST_CASE("Particle: stays on an edge when dist is shorter than the edge") {
    Graph g = graph_from({{{0, 0, 0}}, {{1, 0, 0}}}, {{0, 1}});
    Particle p; p.a = 0; p.b = 1; p.s = 0;
    p.step(g, 0.25);
    CHECK(p.a == 0);
    CHECK(p.b == 1);
    CHECK(p.s == doctest::Approx(0.25));
}

TEST_CASE("Particle: reflects at a degree-1 vertex (2.5 along a unit edge)") {
    // 0->1 (1.0), 1->0 (1.0), then 0.5 into 0->1 again: fully deterministic.
    Graph g = graph_from({{{0, 0, 0}}, {{1, 0, 0}}}, {{0, 1}});
    Particle p; p.a = 0; p.b = 1; p.s = 0;
    p.step(g, 2.5);
    CHECK(p.a == 0);
    CHECK(p.b == 1);
    CHECK(p.s == doctest::Approx(0.5));

    Vec2 pos = p.position(g);
    CHECK(pos.getX() == doctest::Approx(0.5));
    CHECK(pos.getY() == doctest::Approx(0.0));
}

TEST_CASE("Particle: dist exactly equal to the edge length lands on the vertex") {
    Graph g = graph_from({{{0, 0, 0}}, {{1, 0, 0}}}, {{0, 1}});
    Particle p; p.a = 0; p.b = 1; p.s = 0;
    p.step(g, 1.0);
    CHECK(p.a == 1);
    CHECK(p.b == 0);
    CHECK(p.s == doctest::Approx(0.0));
}

TEST_CASE("Particle: zero distance changes nothing") {
    Graph g = graph_from({{{0, 0, 0}}, {{1, 0, 0}}}, {{0, 1}});
    Particle p; p.a = 0; p.b = 1; p.s = 0.3;
    p.step(g, 0.0);
    CHECK(p.a == 0);
    CHECK(p.b == 1);
    CHECK(p.s == doctest::Approx(0.3));
}

TEST_CASE("Particle: always on a real edge, with 0 <= s <= edge length") {
    Graph g = make_x_graph();
    rand_particle = Lcg(7);
    Particle p; p.a = p.b = 4; p.s = 0;
    p.choose_rand_nb(g);

    for (int i = 0; i < 20000; ++i) {
        p.step(g, 0.37);
        REQUIRE(g.are_neighbors(p.a, p.b));
        double len = (g[p.b] - g[p.a]).norm();
        REQUIRE(p.s >= 0.0);
        REQUIRE(p.s <= len + 1e-9);
    }
}

TEST_CASE("Particle: same seed gives the same trajectory") {
    Graph g = make_x_graph();
    auto run = [&](uint64_t seed) {
        rand_particle = Lcg(seed);
        Particle p; p.a = p.b = 4; p.s = 0;
        p.choose_rand_nb(g);
        for (int i = 0; i < 1000; ++i) p.step(g, 0.5);
        return std::make_tuple(p.a, p.b, p.s);
    };
    CHECK(run(42) == run(42));
}

TEST_CASE("Particle: picks all four spokes from a degree-4 vertex") {
    Graph g = make_x_graph();
    rand_particle = Lcg(99);
    std::array<int, 4> hits{};
    for (int i = 0; i < 4000; ++i) {
        Particle p; p.a = p.b = 4; p.s = 0;
        p.choose_rand_nb(g);
        ++hits[p.b];                 // spoke endpoints are vertices 0..3
    }
    for (int c : hits) {
        CHECK(c > 800);              // ~1000 expected each; loose bound
        CHECK(c < 1200);
    }
}

// ---------- ParticleAggregate ----------

TEST_CASE("ParticleAggregate: all mass starts in the start vertex's cell") {
    Graph g = make_x_graph();
    PlanarGraphView view(g, size_t(8), size_t(8));
    rand_particle = Lcg(1);
    ParticleAggregate agg(view, 4, 500);

    auto [cx, cy] = view.to_cell(g[4]);
    CHECK(agg[int(cx), int(cy)] == doctest::Approx(1.0));
    CHECK(total_mass(agg, view) == doctest::Approx(1.0));
    CHECK(nonzero_cells(agg, view) == 1);
}

TEST_CASE("ParticleAggregate: mass is conserved while diffusing") {
    Graph g = make_x_graph();
    PlanarGraphView view(g, size_t(8), size_t(8));
    rand_particle = Lcg(2);
    ParticleAggregate agg(view, 4, 2000);

    for (int step = 0; step < 50; ++step) {
        agg.step_all(0.2);
        if (step % 10 == 0) CHECK(total_mass(agg, view) == doctest::Approx(1.0));
    }
    CHECK(total_mass(agg, view) == doctest::Approx(1.0));
}

TEST_CASE("ParticleAggregate: mass spreads out from the start cell") {
    Graph g = make_x_graph();
    PlanarGraphView view(g, size_t(8), size_t(8));
    rand_particle = Lcg(3);
    ParticleAggregate agg(view, 4, 2000);

    for (int step = 0; step < 100; ++step) agg.step_all(0.2);

    auto [cx, cy] = view.to_cell(g[4]);
    CHECK(nonzero_cells(agg, view) > 1);
    CHECK(agg[int(cx), int(cy)] < 1.0);
}

TEST_CASE("ParticleAggregate: reading twice gives identical values (lazy refresh is stable)") {
    Graph g = make_x_graph();
    PlanarGraphView view(g, size_t(8), size_t(8));
    rand_particle = Lcg(4);
    ParticleAggregate agg(view, 4, 1000);
    agg.step_all(0.3);

    std::vector<double> first, second;
    for (size_t x = 0; x < view.w; ++x)
        for (size_t y = 0; y < view.h; ++y) first.push_back(agg[int(x), int(y)]);
    for (size_t x = 0; x < view.w; ++x)
        for (size_t y = 0; y < view.h; ++y) second.push_back(agg[int(x), int(y)]);
    CHECK(first == second);
}

TEST_CASE("ParticleAggregate: many steps without reading, then one read, is still correct") {
    Graph g = make_x_graph();
    PlanarGraphView view(g, size_t(8), size_t(8));
    rand_particle = Lcg(5);
    ParticleAggregate agg(view, 4, 1000);

    for (int i = 0; i < 200; ++i) agg.step_all(0.25);   // never read in between
    CHECK(total_mass(agg, view) == doctest::Approx(1.0));
    CHECK(nonzero_cells(agg, view) > 1);                // the field really was refreshed
}

TEST_CASE("ParticleAggregate: same seed gives the same field") {
    Graph g = make_x_graph();
    PlanarGraphView view(g, size_t(8), size_t(8));
    auto run = [&](uint64_t seed) {
        rand_particle = Lcg(seed);
        ParticleAggregate agg(view, 4, 500);
        for (int i = 0; i < 30; ++i) agg.step_all(0.3);
        std::vector<double> field;
        for (size_t x = 0; x < view.w; ++x)
            for (size_t y = 0; y < view.h; ++y) field.push_back(agg[int(x), int(y)]);
        return field;
    };
    CHECK(run(11) == run(11));
}

TEST_CASE("ParticleAggregate: point_mass resets the distribution") {
    Graph g = make_x_graph();
    PlanarGraphView view(g, size_t(8), size_t(8));
    rand_particle = Lcg(6);
    ParticleAggregate agg(view, 4, 1000);
    for (int i = 0; i < 100; ++i) agg.step_all(0.3);
    REQUIRE(nonzero_cells(agg, view) > 1);

    agg.point_mass_reset(4, 1000);
    CHECK(nonzero_cells(agg, view) == 1);
    CHECK(total_mass(agg, view) == doctest::Approx(1.0));
}
TEST_CASE("timing: corner fill census on data/010k.vtk" * doctest::skip()) {
    constexpr double Z = 10.0;
    VTKdata raw("data/010k.vtk");
    VTKdata data = raw.filter([=](const Vec3& p) { return p.getZ() >= -Z && p.getZ() <= Z; }, 100);
    Graph g(data);
    const double step = g.avg_norm;     // read before add_intersections shrinks it
    g.add_intersections();

    // Corner pattern bits: SW = 1, SE = 2, NW = 4, NE = 8.
    enum Case { ALL4, THREE, HORIZ, VERT, DIAG, ANTIDIAG, ONE, N_CASES };
    const char* NAME[N_CASES] = {"all four", "three active", "horizontal pair (SW,SE / NW,NE)",
                                 "vertical pair (SW,NW / SE,NE)", "diagonal pair (SW,NE)",
                                 "anti-diagonal pair (SE,NW)", "one active"};
    auto classify = [](int m) {
        switch (m) {
            case 15: return ALL4;
            case 7: case 11: case 13: case 14: return THREE;
            case 3: case 12: return HORIZ;
            case 5: case 10: return VERT;
            case 9: return DIAG;
            case 6: return ANTIDIAG;
            default: return ONE;        // 1, 2, 4, 8 (0 is filtered out before this)
        }
    };

    for (int mode = 0; mode < 2; ++mode) {
        std::cout << "\n=== active = " << (mode == 0 ? "cell contains a vertex" : "cell touched by an edge")
                  << " ===\n";

        for (double mult : {0.25, 0.5, 1.0, 2.0, 4.0}) {          // cell size in units of the average edge
            const double cell = mult * step;
            PlanarGraphView view(g, cell);
            const size_t W = view.w, H = view.h;
            if (W < 2 || H < 2) continue;

            auto idx = [&](Vec2 p) { auto [x, y] = view.to_cell(p); return y * W + x; };

            std::vector<char> active(W * H, 0);
            for (int i = 0; i < int(g.length()); ++i) active[idx(g[i])] = 1;

            if (mode == 1) {                                       // also mark cells an edge passes through
                for (int i = 0; i < int(g.length()); ++i)
                    for (int j : g.neighbors(i)) {
                        if (j < i) continue;
                        Vec2 a = g[i], b = g[j];
                        const int n = int((b - a).norm() / (cell / 4));
                        for (int k = 1; k < n; ++k) active[idx(a + (b - a) * (double(k) / n))] = 1;
                    }
            }

            size_t counts[N_CASES] = {}, touched = 0;
            for (size_t j = 0; j + 1 < H; ++j)
                for (size_t i = 0; i + 1 < W; ++i) {
                    const size_t sw = j * W + i;
                    const int m = active[sw] | (active[sw + 1] << 1) | (active[sw + W] << 2) | (active[sw + W + 1] << 3);
                    if (m == 0) continue;                          // empty corners cost nothing
                    ++counts[classify(m)];
                    ++touched;
                }

            std::cout << std::fixed << std::setprecision(1) << "\ncell = " << mult << " x edge   grid "
                      << W << " x " << H << "   corners touching a filled cell: " << touched << "\n";
            for (int c = 0; c < N_CASES; ++c)
                std::cout << "  " << std::left << std::setw(34) << NAME[c] << std::right
                          << std::setw(9) << counts[c] << "  " << std::setw(5)
                          << (touched ? 100.0 * counts[c] / touched : 0.0) << " %\n";

            size_t sum = 0;
            for (int c = 0; c < N_CASES; ++c) sum += counts[c];
            CHECK(sum == touched);
        }
    }
}

// ---------- benchmark (skipped by default; run with `make bench`) ----------

TEST_CASE("timing: particle walk on data/010k.vtk" * doctest::skip()) {
    using clock = std::chrono::steady_clock;
    auto secs = [](clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    constexpr double Z = 10.0;             // slab half-width; raise for more fibers
    constexpr size_t N_PARTICLES = 100000;
    constexpr int    N_STEPS = 50;
    constexpr size_t GRID = 512;

    auto t0 = clock::now();
    VTKdata raw("data/010k.vtk");
    VTKdata data = raw.filter([=](const Vec3& p) { return p.getZ() >= -Z && p.getZ() <= Z; }, 100);
    Graph g(data);
    auto t1 = clock::now();

    g.add_intersections();
    auto t2 = clock::now();

    PlanarGraphView view(g, GRID, GRID);
    auto t3 = clock::now();

    // start at the degree>=2 vertex closest to the middle of the extent
    const double mx = (g.min_x + g.max_x) / 2, my = (g.min_y + g.max_y) / 2;
    int start = -1;
    double best = 1e300;
    for (int i = 0; i < int(g.length()); ++i) {
        if (g.neighbors(i).size() < 2) continue;
        double d = (g[i] - Vec2(mx, my)).norm2();
        if (d < best) { best = d; start = i; }
    }
    REQUIRE(start >= 0);

    std::cout << std::fixed << std::setprecision(3)
              << "vertices " << g.length() << "  avg edge " << g.avg_norm << "\n"
              << "  load+filter+graph   : " << secs(t0, t1) << " s\n"
              << "  add_intersections   : " << secs(t1, t2) << " s\n"
              << "  PlanarGraphView     : " << secs(t2, t3) << " s\n";

    rand_particle = Lcg(1234);
    ParticleAggregate agg(view, start, int(N_PARTICLES));

    for (double mult : {0.25, 1.0, 4.0}) {          // step length relative to the average edge
        const double dist = mult * g.avg_norm;
        agg.point_mass_reset(start, int(N_PARTICLES));

        auto s0 = clock::now();
        for (int i = 0; i < N_STEPS; ++i) agg.step_all(dist);
        auto s1 = clock::now();

        double first_read = agg[0, 0];              // first read triggers the lazy refresh
        auto s2 = clock::now();
        (void)first_read;

        const double particle_steps = double(N_PARTICLES) * N_STEPS;
        std::cout << "  dist = " << mult << " x avg edge (" << dist << ")\n"
                  << "    walk              : " << secs(s0, s1) << " s  ("
                  << particle_steps / secs(s0, s1) / 1e6 << " M particle-steps/s)\n"
                  << "    field refresh     : " << secs(s1, s2) << " s\n";

        CHECK(total_mass(agg, view) == doctest::Approx(1.0));
    }
}
