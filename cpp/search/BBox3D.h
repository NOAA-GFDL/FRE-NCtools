//
// Created by mzuniga on 3/9/23.
//

#ifndef FREGRID_BBOX3D_H
#define FREGRID_BBOX3D_H

#include <array>
#include <iostream>

namespace nct {
    class BBox3D {
    public:
        BBox3D(const std::array<float, 3> &lo, const std::array<float, 3> &hi) : lo(lo), hi(hi) {}
        static bool intersect(const BBox3D& A, const BBox3D& B) {
            if (A.lo[0] > B.hi[0] || A.lo[1] > B.hi[1] || A.lo[2] > B.hi[2] ||
                A.hi[0] < B.lo[0] || A.hi[1] < B.lo[1] || A.hi[2] < B.lo[2]) {
                return false;
            } else {
                return true;
            }
        }
    private:
        std::array<float,3> lo;
        std::array<float,3> hi;
    };

    class BoxPair {
    public:
        int id;
        BBox3D& bb;
        BoxPair(int id, BBox3D &bb) : id(id), bb(bb) {}
    };

} // nct

#endif //FREGRID_BBOX3D_H
