#include "graph.hpp"
#include <algorithm>
#include <numeric>


Graph::Graph(const VTKdata& data)
    : points(data.points()), adjacency_list(data.points().size()) {
    for (const auto& line : data.lines()) {
        for (size_t j = 0; j + 1 < line.size(); ++j) {
            int a = line[j], b = line[j + 1];
            if (a == b) continue;
            adjacency_list[a].push_back(b);
            adjacency_list[b].push_back(a);
        }
    }
    for (auto& nb : adjacency_list) {
        std::sort(nb.begin(), nb.end());
        nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    }
    refresh_convenience_data();
}

void Graph::refresh_convenience_data(){
    min_x = (*this)[0].getX(), min_y = (*this)[0].getY(), max_x = (*this)[0].getX(), max_y = (*this)[0].getY();
    min_norm = std::numeric_limits<double>::max(), max_norm = std::numeric_limits<double>::lowest();
    double norm_sum = 0;
    int edge_count = 0;
    for (int i = 0; i < length(); i++){
        min_x = std::min((*this)[i].getX(), min_x);
        min_y = std::min((*this)[i].getY(), min_y);
        max_x = std::max((*this)[i].getX(), max_x);
        max_y = std::max((*this)[i].getY(), max_y);
        for (int j : neighbors(i)){
            min_norm = std::min(((*this)[i] - (*this)[j]).norm(), min_norm);
            max_norm = std::max(((*this)[i] - (*this)[j]).norm(), max_norm);
            norm_sum += ((*this)[i] - (*this)[j]).norm();
            edge_count += 1;
        }
    }
    avg_norm = norm_sum/edge_count;
}


constexpr int INTERSECTION_Z = 378;
void Graph::add_intersection(std::pair<int, int> l1, std::pair<int, int> l2, Vec2 v){
    std::ranges::replace(neighbors(l1.first), l1.second, points.size());
    std::ranges::replace(neighbors(l1.second), l1.first, points.size());
    std::ranges::replace(neighbors(l2.first), l2.second, points.size());
    std::ranges::replace(neighbors(l2.second), l2.first, points.size());

    //THE PUSH BACK MUST COME SECOND
    points.push_back(Vec3(v.getX(), v.getY(), INTERSECTION_Z));
    adjacency_list.push_back(std::vector({l1.first, l1.second, l2.first, l2.second}));
}

//HELPERS FOR ADD INTERSECTION
//-------------------------------
double cross(const Vec2& a, const Vec2& b, const Vec2& c) {
    Vec2 ab = b - a, ac = c - a;
    return ab.getX() * ac.getY() - ab.getY() * ac.getX();
}

Vec2 intersection_helper(const Graph& g, std::pair<int, int> l1, std::pair<int, int> l2, bool& valid) {
    valid = false;

    // Shared endpoint: already connected at a vertex.
    if (l1.first == l2.first || l1.first == l2.second ||
        l1.second == l2.first || l1.second == l2.second)
        return Vec2();

    Vec2 a = g[l1.first], b = g[l1.second];
    Vec2 c = g[l2.first], d = g[l2.second];

    const double len_ab = (b - a).norm();
    const double len_cd = (d - c).norm();
    if (len_ab == 0 || len_cd == 0) return Vec2();

    double d1 = cross(a, b, c), d2 = cross(a, b, d);
    double d3 = cross(c, d, a), d4 = cross(c, d, b);

    constexpr double DIST_EPS = 1e-6;   
    if (std::abs(d1) <= DIST_EPS * len_ab || std::abs(d2) <= DIST_EPS * len_ab ||
        std::abs(d3) <= DIST_EPS * len_cd || std::abs(d4) <= DIST_EPS * len_cd)
        return Vec2();

    bool proper = ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
                  ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0));
    if (!proper) return Vec2();

    valid = true;
    return a + (b - a) * (d3 / (d3 - d4));
}
//------------------------------------
namespace {
// Fraction of the way from p to q at which x lies (x on segment pq).
double param_along(const Vec2& p, const Vec2& q, const Vec2& x) {
    return (x - p).dot(q - p) / (q - p).norm2();
}
}

void Graph::add_intersections() {
    refresh_convenience_data();
    PlanarGraphView view(*this, avg_norm / 4);   // cell size is only a speed knob
    const size_t W = view.w, H = view.h;

    struct CellBox { size_t x0, x1, y0, y1; };
    std::vector<std::pair<int,int>> edges;
    std::vector<CellBox> boxes;
    std::vector<std::vector<int>> bins(W * H);   // edge ids

    for (size_t i = 0; i < length(); i++)
        for (int j : neighbors(i)) {
            if (int(i) > j) continue;
            auto [ax, ay] = view.to_cell((*this)[i]);
            auto [bx, by] = view.to_cell((*this)[j]);
            CellBox box{std::min(ax, bx), std::max(ax, bx),
                        std::min(ay, by), std::max(ay, by)};
            int id = int(edges.size());
            edges.push_back({int(i), j});
            boxes.push_back(box);
            for (size_t x = box.x0; x <= box.x1; x++)
                for (size_t y = box.y0; y <= box.y1; y++)
                    bins[y * W + x].push_back(id);
        }

    // 1. Collect every crossing once. The graph is not modified here.
    struct Crossing { int e, f; Vec2 p; };
    std::vector<Crossing> found;
    for (size_t x = 0; x < W; x++)
      for (size_t y = 0; y < H; y++) {
        const auto& cell = bins[y * W + x];
        for (size_t m = 0; m < cell.size(); m++)
          for (size_t n = m + 1; n < cell.size(); n++) {
            int e = cell[m], f = cell[n];
            if (x != std::max(boxes[e].x0, boxes[f].x0) ||
                y != std::max(boxes[e].y0, boxes[f].y0)) continue;   // test each pair once
            bool ok;
            Vec2 p = intersection_helper(*this, edges[e], edges[f], ok);
            if (ok) found.push_back({e, f, p});
          }
      }

    // 2. Create the new vertices, recording where each lies along its two edges.
    std::vector<std::vector<std::pair<double,int>>> splits(edges.size());
    for (const Crossing& k : found) {
        int x = int(points.size());
        points.push_back(Vec3(k.p.getX(), k.p.getY(), INTERSECTION_Z));
        adjacency_list.emplace_back();
        auto [a, b] = edges[k.e];
        auto [c, d] = edges[k.f];
        splits[k.e].push_back({param_along((*this)[a], (*this)[b], k.p), x});
        splits[k.f].push_back({param_along((*this)[c], (*this)[d], k.p), x});
    }

    // 3. Rebuild each split edge once, as a, x1, x2, ..., b.
    for (size_t id = 0; id < edges.size(); id++) {
        auto& list = splits[id];
        if (list.empty()) continue;
        auto [a, b] = edges[id];
        std::erase(adjacency_list[a], b);
        std::erase(adjacency_list[b], a);
        std::ranges::sort(list);
        int prev = a;
        for (auto& [t, x] : list) {
            adjacency_list[prev].push_back(x);
            adjacency_list[x].push_back(prev);
            prev = x;
        }
        adjacency_list[prev].push_back(b);
        adjacency_list[b].push_back(prev);
    }

    refresh_convenience_data();
}

Vec2 PlanarGraphView::to_world(size_t ix, size_t iy) const {
    return Vec2(g.min_x + ix * cell_w_, g.min_y + iy * cell_h_);   
}


std::pair<size_t, size_t> PlanarGraphView::to_cell(Vec2 p) const {
    double fx = (p.getX() - g.min_x) / cell_w_;
    double fy = (p.getY() - g.min_y) / cell_h_;
    size_t ix = fx <= 0 ? 0 : std::min(size_t(fx), w - 1);
    size_t iy = fy <= 0 ? 0 : std::min(size_t(fy), h - 1);
    return {ix, iy};
}


//must have recent convenience data
PlanarGraphView::PlanarGraphView(Graph& _graph, size_t _w, size_t _h)
    : g(_graph), w(_w), h(_h),
      spacially_organized_points(_w*_h, std::vector<int>()),
      view(spacially_organized_points.data(), h, w){


    cell_w_ = (g.max_x - g.min_x) / w;
    cell_h_ = (g.max_y - g.min_y) / h;

    for (int i = 0; i < g.length(); i++){
        (*this)[to_cell(g[i])].push_back(i);
    }
}

//must have recent convenience data
PlanarGraphView::PlanarGraphView(Graph& g, double resolution)
: PlanarGraphView(g, (((size_t)((g.max_x - g.min_x)/resolution)) + 1), ((size_t)((g.max_y - g.min_y)/resolution)) + 1){}



