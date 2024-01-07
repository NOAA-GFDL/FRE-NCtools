
#ifndef FREGRID_BBOX3D_H
#define FREGRID_BBOX3D_H

#include <array>
#include <iostream>
#include <limits>
#include <vector>
#include <algorithm>
//#include <ranges>
#include <algorithm>
#include <cmath>
#include <format>

#include "DistanceInterval.h"
#include "Point3D.h"
#include "Polygon.h"

namespace nct {
//using namespace oneapi::dpl::experimental::ranges;
    //using std::ranges::minmax_element; //cannot find using oneapi
    using std::vector;

    class BBox3D {
    private:
        //bounding boxes are always made up of floats (i.e. double
        // precision is not needed. See functions expand_for_doubles
        std::array<float, 3> lo;
        std::array<float, 3> hi;
    public:
        inline float getLo(int dim) { return lo[dim]; }
        inline float getHi(int dim) { return hi[dim]; }

      void operator=(const BBox3D& b)
      {
        for (auto i = 0; i < 3; i++) {
          lo[i] = b.lo[i];
          hi[i] = b.hi[i];
        }
      }


        //Constructor from a mesh of polygons, which in turn is really
        //  just a vector of pointers to points.
        //  This is the constructor for use with legacy nctools
        explicit BBox3D(MeshPolygon<double> &poly) {
          unsigned int dim;
          auto comp{
                  [&dim](const Point3D<double> *a, const Point3D<double> *b) {
                      return (a->p[dim] < b->p[dim]);
                  }};

          for (dim = 0; dim < 3; ++dim) {
            auto [min, max] = std::minmax_element(poly.points.begin(), poly.points.end(), comp);
            lo[dim] = static_cast<float>((*min)->p[dim]);
            hi[dim] = static_cast<float>((*max)->p[dim]);
          }
          expand_for_doubles();
        }

        /**
         * Constructor from no polygons. User is expected to incrementally expand the box
         * with the expand method.
         */
        explicit BBox3D() {
          for (auto i = 0; i < 3; ++i) {
            lo[i] = std::numeric_limits<float>::max();
            hi[i] = -std::numeric_limits<float>::max();
          }
        }

        //Expand box so that it contain point p
        template <class T>
        void expand(const std::array<T, 3>& p) {
          for (auto i = 0; i < 3; i++) {
            if (p[i] < lo[i]) lo[i] = p[i];
            if (p[i] > hi[i]) hi[i] = p[i];
          }
        }

        friend std::ostream &operator<<(std::ostream &os, const BBox3D &b) {
          for(int i = 0; i<3; i++){
            auto str = std::format("[ {:16.10e},{:16.10e} ], ", b.lo[i], b.hi[i]);
            os << str << std::endl;
          }
          return os;
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

        inline static bool intersect(const BBox3D &A, DistanceInterval<float> &di, const int dim) {
          if (A.lo[dim] > di.getFar() || A.hi[dim] < di.getNear()) {
            return false;
          } else {
            return true;
          }
        }

        /**
         * Returns true if the  projection of BBOX onto a Z=const plane contains
         * the projection of the point. In this test, containment is strict and
         * points just on boundaries are not considered contained.
         * @param bb the bounding box
         * @param p the point
         * @return
         */
        template<class T>
        inline static bool contains_zk(const BBox3D &bb, const std::array<T, 3> &p) {
          if (bb.lo[0] < p[0] && p[0] < bb.hi[0] &&
              bb.lo[1] < p[1] && p[1] < bb.hi[1]) {
            return true;
          } else {
            return false;
          }
        }

      template<class T>
      inline static bool contains(const BBox3D &bb, const std::array<T, 3> &p) {
        if (bb.lo[0] <= p[0] && p[0] <= bb.hi[0] &&
            bb.lo[1] <= p[1] && p[1] <= bb.hi[1] &&
            bb.lo[2] <= p[2] && p[2] <= bb.hi[2]) {
          return true;
        } else {
          return false;
        }
      }


      // Expand the bbox to take into account the RANGE_CHECK_CRITERIA used in function
      // clip_2dx2d_great_circle. Note that the recommended RANGE_CHECK_CRITERIA may be a
      // magic number, and so this expansion routine is also ad-hock. Its goal is to make
      // bounding boxes large enough to not miss anything.
      inline void expand_for_rcc(float rcc = 0.05) {
        for (int i = 0; i < 3; i++) {
          auto diff = rcc * (hi[i] - lo[i]);
          //diff = std::max(diff, 1000.0f);
          hi[i] += diff;
          lo[i] -= diff;
        }
      }

        //Expand the box by "+- next possible float" in all six interval ends.
        // This is useful as boxes lo[] and hi[] are float; but the boxes are set
        // or expanded with data that is doubles.
        inline void expand_for_doubles(){
          for (int i = 0; i < 3; i++) {
            hi[i] = std::nextafter(hi[i], std::numeric_limits<float>::max());
            lo[i] = std::nextafter(lo[i], -std::numeric_limits<float>::max());
          }
        }

      //Expand the box by "+- next possible float" in all six interval ends.
      // This is useful as boxes lo[] and hi[] are float; but the boxes are set
      // or expanded with data that is doubles.
      // Function actually only does something when compiled with option USE_NEXTAFTER.
      inline void expand_for_doubles_if(){
#ifdef USE_NEXTAFTER
        for (int i = 0; i < 3; i++) {
            hi[i] = std::nextafter(hi[i], std::numeric_limits<float>::max());
            lo[i] = std::nextafter(lo[i], -std::numeric_limits<float>::max());
          }
#endif
      }

    };
} // nct

#endif //FREGRID_BBOX3D_H
