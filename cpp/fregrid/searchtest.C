#include <iostream>
#include <vector>
#include <chrono>
//#include <mdspan> //not yet avail w/ gcc!
#include <array>
#include <algorithm>

#ifdef USE_SYCL
#include <CL/sycl.hpp>
#endif

#include "BBox3D.h"
#include "Polygon.h"

#include "BoxedObj.h"
#include "DITree.h"
#include "BruteBoxQuery.h"

using BBox_t = nct::BBox3D;
using BPair_t  = nct::BoxAndId;
using Poly_t = nct::MeshPolygon<double>;
using Point_t = nct::Point3D<double>;
using nct::BruteBoxQuery;
using nct::fill_value_size_t;

const int NX{100};
const int NY{100};
const int NZ{4};

//constexpr long long size = 1'000'000'000;
//TODO: this will be replaced by unit or functionality tests?

using namespace std;
using namespace std::chrono;

template <typename Func>
void runAndTime(const std::string& name, Func func){
    cout << "Starting " << name << endl;
    const auto start =  high_resolution_clock::now();
    func();
    //const std::chrono::duration<double> dur = std::chrono::steady_clock::now() - sta;
    const auto dt = duration_cast<microseconds>( high_resolution_clock::now() - start);
    cout << "For " << name << " time is : " << (dt.count() / 1.0e6) << " sec. " << endl;
}

void compare_results(std::vector<std::vector<size_t>> &v1,
                     std::vector<std::vector<size_t>> &v2);

size_t pt_idx(int i, int j, int k) {
    return (k * (NX+1) * (NY+1)  + j * (NX+1)  + i);
}

extern void compare_results(const std::vector<std::vector<size_t>> &v1,
                     const std::vector<std::vector<size_t>> &v2);

void create_points_and_polygons(std::vector<Point_t> &points, std::vector<Poly_t> &polys);

int main() {
    using namespace std;
    using namespace std::chrono;
    using namespace nct;
    vector<Point_t> points;
    vector<Poly_t> polys;
    vector<BBox_t> boxes;
    vector<BBox_t> qBoxes;
    vector<BPair_t> bPairs;
    vector<BPair_t> qbPairs;//box pairs  that are used as query boxes

    const size_t NCells{NX * NY* NZ};
   // const size_t NPoints{(NX + 1)* (NY + 1) * (NZ + 1)};//does not wrap around.
    //const size_t NPoly {NCells};

    cout << "Hello from searchtest.x !" << endl;

    //The legacy fregrid maintains an array of latitudes and an
    //equal array of longitudes, and a polygon is (normally)0 four of these pairs.
    //Here we start with the lat-lon arrays already converted to a point and all
    //stored in an array.
    create_points_and_polygons(points, polys);

    // Make the boxes and pairs to be searched.
    //Boxes will b made from the polygons.
    //Note we require the size to be reserved
    boxes.reserve(polys.size());
    bPairs.reserve(polys.size());
    for(size_t  i = 0; i< polys.size(); ++i){
       boxes.emplace_back(polys[i]);
       bPairs.emplace_back(i, &boxes[i]);
    }

    //Make the boxes and pairs that are going to be used
    // as query boxes. In this case just copy those that
    // are to be searched.
    qBoxes.reserve(boxes.size());
    qbPairs.reserve(boxes.size());
    for (const auto & abox : boxes) {
        qBoxes.push_back(abox);
    }
    for(size_t i = 0; i< qBoxes.size(); ++i){
       qbPairs.emplace_back(i, &qBoxes[i]);
    }

    // Make some results vectors
    vector<vector<size_t>> results_t(qBoxes.size(), vector<size_t>());
    vector<vector<size_t>> results_b(qBoxes.size(), 
        vector<size_t>(MAX_NN,fill_value_size_t));

  //Make the BF search data structure


  //std::cout<< "starting box query with subtype="
   // << BruteBoxQuery::search_subtype_str[static_cast<int>(bbq_subtype)] <<std::endl;


  auto start =  high_resolution_clock::now();
  #ifdef USE_SYCL
    std::vector<size_t> mdresults;
  std::vector<size_t> res_count;
  bbq.search_sycl(qBoxes, mdresults,  res_count, 10);
  //put results in the old form
  //for(int i = 0; i< qBoxes.size(); i++) {
   // for (int j = 0; j < res_count[i]; j++) {
      //auto idx_md = BruteBoxQuery::qr_idx(j, i, 10);
     //results_b[i].push_back(mdresults[idx_md]);
    //}
 //}
  #else
    //bbq.search(qBoxes, results_b, bbq_subtype);
  #endif

    BruteBoxQuery bbq{bPairs, boxes};
    auto bbq_subtype = BruteBoxQuery::search_subtype::copy_if;
    auto qname = std::string( BruteBoxQuery::search_subtype_str[static_cast<int>(bbq_subtype)]);
    qname = "BF box queries with subtype=" + qname;
    runAndTime(qname, [=,&bbq, &qBoxes, &results_b]() mutable {
        bbq.search(qBoxes, results_b, bbq_subtype);}
        );

    //Make the tree data structure.
    DITree<BoxAndId> tree (bPairs);
    qname = "DITree box queries";
    runAndTime(qname,  [=,&tree, &qbPairs, &results_t]() mutable{
        tree.search(qbPairs, results_t );}
    );

    //check tree awnsers against brute force.
    compare_results(results_b, results_t);
}

void create_points_and_polygons(std::vector<Point_t> &points, std::vector<Poly_t> &polys) {
    for (auto k = 0; k <= NZ; ++k){
        for(auto j = 0; j<=NY; ++j){
            for (auto i = 0; i<=NX; ++i){
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
}

void compare_results(std::vector<std::vector<size_t>> &v1,
                     std::vector<std::vector<size_t>> &v2) {
    std::cout << "Compare results s1 s2 " << v1.size() << ", " << v2.size() << std::endl;
    for (size_t i = 0; i < v1.size(); ++i) {
       std::sort(v1[i].begin(), v1[i].end());
       std::sort(v2[i].begin(), v2[i].end());
       std::vector<int> intersection;
       intersection.clear();
       if (i == 0) {
            std::cout << "i1 s1 s2 " << v1[i].size() << ", " << v2[i].size() << std::endl;
       }
       if (i == (v1.size() - 1)) {
            std::cout << "1L s1 s2 " << v1[i].size() << ", " << v2[i].size() << std::endl;
       }
       std::set_intersection(v1[i].begin(), v1[i].end(), v2[i].begin(), v2[i].end(),
                             std::back_inserter(intersection));
       if ((intersection.size() != v1[i].size()) || (v1[i].size() != v2[i].size())) {
            std::cout << "answer diff for i = " << i << std::endl;
       }
    }
}
