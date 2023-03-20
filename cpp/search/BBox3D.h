
#ifndef FREGRID_BBOX3D_H
#define FREGRID_BBOX3D_H

#include <array>
#include <iostream>
#include <limits>
#include <vector>
#include <ranges>
#include <algorithm>
#include <cmath>

#include "Point3D.h"
#include "Polygon.h"

namespace nct {
    using std::ranges::minmax_element;
    using std::vector;

    class BBox3D {
    private:
        //bounding boxes are always made up of floats.
        std::array<float, 3> lo;
        std::array<float, 3> hi;
    public:
        //Constructor from a mesh of polygons, which in turn is really
        //  just a vector of pointers to points.
        //  This is the constructor for use with legacy nctools, where
        explicit BBox3D(MeshPolygon<double> &poly)  {
            unsigned int dim;
            auto comp {
                    [&dim](const Point3D<double> *a, const Point3D<double> *b) {
                return (a->p[dim] < b->p[dim]);}};

            for (dim = 0; dim < 3; ++dim) {
                auto [min, max] = minmax_element(poly.points, comp);
                lo[dim] = static_cast<float>((*min)->p[dim]);
                hi[dim] = static_cast<float>((*max)->p[dim]);
            }

            //The hi of the boxes normally need to increased do to casting of double to float.
#ifndef USE_NEXTAFTER
            for (auto &v: hi) {
                v = std::nextafter(v, std::numeric_limits<float>::max());
            }
#endif
        }

        //TODO: investigate using sum of bools. e.g. int r= bool(A.lo[0] > B.hi[0])
        // + bool(A.lo[1] > B.hi[1]) ... check for branch misses w perf stat
        inline static bool intersect(const BBox3D &A, const BBox3D &B) {
            if (A.lo[0] > B.hi[0] || A.lo[1] > B.hi[1] || A.lo[2] > B.hi[2] ||
                A.hi[0] < B.lo[0] || A.hi[1] < B.lo[1] || A.hi[2] < B.lo[2]) {
                return false;
            } else {
                return true;
            }
        }

        //augment box a with extents of box B
        static void augment(BBox3D &A, BBox3D &B) {
            for (int i = 0; i < 3; i++) {
                if (B.lo[i] < A.lo[i])
                    A.lo[i] = B.lo[i];
                if (B.hi[i] > A.hi[i])
                    A.hi[i] = B.hi[i];
            }
        }

        //augment box A with extents of box B in dimension d
        static void augment(BBox3D &A, BBox3D &B, const int d) {
            if (B.lo[d] < A.lo[d])
                A.lo[d] = B.lo[d];
            if (B.hi[d] > A.hi[d])
                A.hi[d] = B.hi[d];
        }
    };

    class BoxPair {
    public:
        int id;
        BBox3D& bb;
        BoxPair(int id, BBox3D &bb) : id(id), bb(bb) {}
    };
} // nct

#endif //FREGRID_BBOX3D_H
