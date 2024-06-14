#pragma once

#include <cstddef>
#include <cstdint>
#include <spatialindex/SpatialIndex.h>
std::vector<uint64_t> rstar_res;
int total_cnt=0;
namespace rs
{
    template <class KeyType>
    struct Coord
    {
        KeyType x;
        double y;
    };

    struct FitParam
    {
        double slope;
        double intercept;
    };

} // namespace rs
