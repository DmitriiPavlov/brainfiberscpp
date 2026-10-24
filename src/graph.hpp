#pragma once
#include "vtk_data.hpp"
#include <mdspan>
#include <algorithm>

using SpacialView = std::mdspan<std::vector<int>, std::dextents<size_t, 2>>;
class Graph{
public:
    //builds an adjacency list out of the graph
    explicit Graph(const VTKdata& data);

    //methods delegating to underlying vector
    Vec2 operator[](int i) const{ return Vec2(points[i].getX(), points[i].getY()); }
    size_t length() const {return points.size();}


    //methods for neighbor look ups
    std::vector<int>& neighbors(int point_index) {return adjacency_list[point_index];}
    bool are_neighbors(int i1, int i2) {return std::ranges::contains(neighbors(i1), i2);}
    
    //adds virtual points to the adjacency and point lists
    void add_intersections();


    //some convenient data to know
    double min_norm;
    double max_norm;
    double avg_norm;
    double min_x;
    double min_y;
    double max_x;
    double max_y;

    void refresh_convenience_data();



private:
    Graph() = default;
    
    std::vector<Vec3> points;
    std::vector<std::vector<int>> adjacency_list;

    void add_intersection(std::pair<int, int> l1, std::pair<int,int> l2, Vec2 v);

};


class PlanarGraphView{
public:
    size_t w, h;

    //both require recently updated convenience data
    explicit PlanarGraphView(Graph& graph, size_t w, size_t h);
    //best faith resolution, in reality each cell is going to be a little smaller
    explicit PlanarGraphView(Graph& graph, double resolution);

    Vec2 operator[](int i) const{return  g[i];}
    std::vector<int>& operator[](std::pair<size_t, size_t> cell) const {return  view[cell.second, cell.first];}
    std::vector<int>& operator[](size_t x, size_t y) const {return  view[y,x];}


    std::vector<int>& neighbors(int point_index) {return g.neighbors(point_index);}
    size_t length() const {return g.length();}

    //should make sure that g has refreshed data
    std::pair<size_t, size_t> to_cell(Vec2 v) const;
    Vec2 to_world(size_t x, size_t y) const;

private:
    double cell_w_, cell_h_;
    Graph& g;
    std::vector<std::vector<int>> spacially_organized_points;
    const SpacialView view;
};
