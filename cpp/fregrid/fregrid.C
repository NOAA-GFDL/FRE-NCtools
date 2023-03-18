#include <iostream>
#include <vector>
//#include <mdspan> //not yet avail w/ gcc!
#include <array>
#include "BBox3D.h"
#include "Polygon.h"
#include "BruteBoxQuery.h"

using BBox_t = nct::BBox3D;
using BPair_t  = nct::BoxPair;
using Poly_t = nct::MeshPolygon<double>;
using Point_t = nct::Point3D<double>;
using nct::BruteBoxQuery;

const size_t NX{4};
const size_t NY{3};
const size_t NZ{2};

size_t pt_idx(int i, int j, int k) {
    return (k * (NX+1) * (NY+1)  + j * (NX+1)  + i);
}

int main() {
    using namespace std;
    using namespace nct;
    vector<Point_t> points;
    vector<Poly_t> polys;
    vector<BBox_t> boxes;
    vector<BBox_t> qBoxes;
    vector<BPair_t> bpairs;

    const size_t NCells{NX * NY* NZ};
    const size_t NPoints{(NX + 1)* (NY + 1) * (NZ + 1)};//does not wrap around.
    const size_t NPoly {NCells};

    cout << "Hello!" << endl;

    //The legacy fregrid maintains an array of latitudes and an
    //equal array of longitudes, and a polygon is (normally)0 four of these pairs.
    //Here we start with the lat-lon arrays already converted to a point and all
    //stored in an array.
    for (int k = 0; k<=NZ; ++k){
        for(int j = 0; j<=NY; ++j){
            for (int i = 0; i<=NX; ++i){
                points.emplace_back(i, j, k);
            }
        }
    }

    //Make the polygons from points - similar to legacy fregrid
    //mdspan not avail in gcc 12.?
    //auto pts = std::mdspan(points.data(), NX, NY, NZ);
    for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
             for (int i = 0; i < NX; ++i) {
                Poly_t p;
                //polygon is counter clockwise stating from lower-left point
                p.push_back(&points[pt_idx(i, j, k)]);
                p.push_back(&points[pt_idx(i, j + 1, k)]);
                p.push_back(&points[pt_idx(i + 1, j + 1, k)]);
                p.push_back(&points[pt_idx(i + 1, j, k)]);

                //p.push_back(&pts[i, j, k]);
                //p.push_back(&pts[i, j+1, k]);
                //p.push_back(&pts[i + 1, j + 1, k]);
               // p.push_back(&pts[i + 1, j, k]);
                polys.emplace_back(p);
            }
        }
    }

    //Make the boxes and pairs from the polygons
    for(auto i = 0; i< polys.size(); ++i){
        boxes.emplace_back(polys[i]);
        bpairs.emplace_back(i, boxes[i]);
    }

    //Make the search data structure
    BruteBoxQuery bbq{bpairs, boxes};

    //make some query boxes:
    // In this first test they are just copies of those in boxes.
    for(auto i = 0; i<boxes.size(); ++i) {
        qBoxes.push_back(boxes[i]);
    }

    //query the boxes with the query boxes
    vector<vector<size_t>> results;
    results.reserve(qBoxes.size());
    for(int i = 0; i< qBoxes.size(); i++) {
        vector<size_t> v;
        results.push_back(v);
    }
    bbq.search(qBoxes, results);

    cout <<" Bye Bye" << endl;

    //Check the results.
}
