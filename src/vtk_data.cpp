#include "vtk_data.hpp"
#include <fstream>
#include <stdexcept>

VTKdata::VTKdata(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("cannot open " + filename);

    std::string s;
    for (int i = 0; i < 4; ++i)
        if (!std::getline(in, s)) throw std::runtime_error("truncated VTK header");

    int n_points;
    if (!(in >> s) || s != "POINTS" || !(in >> n_points >> s) || n_points < 0)
        throw std::runtime_error("bad POINTS header");

    points_.reserve(n_points);
    for (int i = 0; i < n_points; ++i) {
        double x, y, z;
        if (!(in >> x >> y >> z)) throw std::runtime_error("bad point data");
        points_.push_back(Vec3(x, y, z));
    }

    int n_lines, total;
    if (!(in >> s) || s != "LINES" || !(in >> n_lines >> total) || n_lines < 0)
        throw std::runtime_error("bad LINES header");

    lines_.reserve(n_lines);
    for (int i = 0; i < n_lines; ++i) {
        int len;
        if (!(in >> len) || len < 0) throw std::runtime_error("bad line length");
        std::vector<int> line(len);
        for (int& ref : line) {
            if (!(in >> ref) || ref < 0 || ref >= n_points)
                throw std::runtime_error("bad point index in line");
        }
        lines_.push_back(std::move(line));
    }
}

VTKdata VTKdata::filter(const VTKfilter& keep, size_t min_line_length) const {
    std::vector<char> pass(points_.size());
    for (size_t i = 0; i < points_.size(); ++i) pass[i] = keep(points_[i]);

    VTKdata out;
    std::vector<int> remap(points_.size(), -1);   // old id -> new id

    auto flush = [&](std::vector<int>& run) {
        if (run.size() >= min_line_length && !run.empty()) {
            for (int& old_id : run) {
                if (remap[old_id] < 0) {
                    remap[old_id] = int(out.points_.size());
                    out.points_.push_back(points_[old_id]);
                }
                old_id = remap[old_id];
            }
            out.lines_.push_back(std::move(run));
        }
        run.clear();
    };

    std::vector<int> run;
    for (const auto& line : lines_) {
        for (int id : line) {
            if (pass[id]) run.push_back(id);
            else flush(run);          
        }
        flush(run);
    }
    return out;
}
