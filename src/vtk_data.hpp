#pragma once
#include "linear_algebra.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>




using VTKfilter = std::function<bool(const Vec3&)>;

// True iff minz <= p.z <= maxz. Usage: zfilter<-1.0, 1.0>
template <double minz, double maxz>
bool zfilter(const Vec3& p) {
    return p.getZ() >= minz && p.getZ() <= maxz;
}

class VTKdata {
public:
    explicit VTKdata(const std::string& filename);

    //keeps only points where keep(v) = true, and throws away any lines below minimum length
    VTKdata filter(const VTKfilter& keep, size_t min_line_length) const;

    const std::vector<Vec3>& points() const { return points_; }
    const std::vector<std::vector<int>>& lines() const { return lines_; }

private:
    VTKdata() = default;   // used by filter()

    std::vector<Vec3> points_;
    std::vector<std::vector<int>> lines_;
};
