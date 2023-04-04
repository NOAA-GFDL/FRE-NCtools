
#ifndef FREGRID_POLYGON_H
#define FREGRID_POLYGON_H

#include <vector>
#include "BBox3D.h"

namespace nct {
/*
 This MeshPolygon intented for use with the "loose" polygon representation of three vectors
 of NCTools.
 * */
    template<class T>
    class MeshPolygon {
    public:
        std::vector<Point3D<T>*> points;
        //The polygon is a vector of pointers to points in a mesh
        // (which in turn is just an array of points)

        //MeshPolygon(const MeshPolygon& mp) : points{mp.points}{}
        //explicit MeshPolygon(std::vector<Point3D<T>*> const & pv) : points{pv}{}
        MeshPolygon() {}
        ~MeshPolygon(){}
        void push_back(Point3D<T> * pt){
            points.push_back( pt);
        }
        void emplace_back(Point3D<T> * pt){points.emplace_back(pt);}
        //size_t getSize() const { return points.size();  }
    };
}

#endif //FREGRID_POLYGON_H
