
#ifndef FREGRID_POINT3D_H
#define FREGRID_POINT3D_H

#include <array>

namespace nct {

    template <class T>
    class Point3D {
    public:
        std::array<T,3> p;
        explicit Point3D(const std::array<T,3> &v) : p(v) {}
        explicit Point3D(T x, T y, T z) {
            p[0]=x; p[1]=y; p[2]=z;
        }
    };
} // nct
#endif //FREGRID_POINT3D_H
