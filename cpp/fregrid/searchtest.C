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

const size_t NX{1000};
const size_t NY{1000};
const size_t NZ{4};

void compare_results(std::vector<std::vector<size_t>> &v1,
                     std::vector<std::vector<size_t>> &v2);

size_t pt_idx(int i, int j, int k) {
    return (k * (NX+1) * (NY+1)  + j * (NX+1)  + i);
}

extern void compare_results(const std::vector<std::vector<size_t>> &v1,
                     const std::vector<std::vector<size_t>> &v2);

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

    // Make the boxes and pairs to be searched.
    //Boxes will b made from the polygons.
    //Note we require the size to be reserved
    boxes.reserve(polys.size());
    bPairs.reserve(polys.size());
    for(unsigned int  i = 0; i< polys.size(); ++i){
       boxes.emplace_back(polys[i]);
       bPairs.emplace_back(i, &boxes[i]);
    }

    //Make the boxes and pairs that are going to be used
    // as query boxes. In this case just copy those that
    // are to be searched.
    for (const auto & abox : boxes) {
        qBoxes.push_back(abox);
    }
    qbPairs.reserve(bPairs.size());
    for(size_t i = 0; i< bPairs.size(); ++i){
        qbPairs[i] =  bPairs [i];
    }

  //Make some results vectors
    vector<vector<size_t>> results_b;
    vector<vector<size_t>> results_t;
    results_b.reserve(qBoxes.size());
    results_t.reserve(qBoxes.size());
     for(int i = 0; i< qBoxes.size(); i++) {
      vector<size_t> vb;
      vector<size_t> vt;
         results_b.push_back(vb);
         results_t.push_back(vt);
      }

  //Make the BF search data structure
  std::cout<< "starting box query"<<std::endl;
  BruteBoxQuery bbq{bPairs, boxes};

  std::vector<size_t> mdresults;
  std::vector<size_t> res_count;
  auto start =  high_resolution_clock::now();
  #ifdef USE_SYCL
  bbq.search_sycl(qBoxes, mdresults,  res_count, 10);
  #else
  //bbq.search(qBoxes, mdresults,  res_count, 10);
  bbq.search_std_partition(qBoxes[0],results_b[0]);
  #endif

  auto stop =  high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "BBox search_sycl time"
       << duration.count() / 1.0e6 << " seconds" << endl;
  std::cout<< "finished box query"<<std::endl;

  //put results in the old form
  for(int i = 0; i< qBoxes.size(); i++) {
    for (int j = 0; j < res_count[i]; j++) {
      auto idx_md = BruteBoxQuery::qr_idx(j, i, 10);
      results_b[i].push_back(mdresults[idx_md]);
    }
  }

  //Make the tree data structure.
    DITree<BoxAndId> tree (bPairs);

    //std::vector<size_t> stResult;
    //BPair_t bp(0, &(boxes[0]));
    //tree.search(bp, stResult);

    start =  high_resolution_clock::now();
    tree.search(bPairs, results_t );
    stop =  high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    cout << "Tree search time"
       << duration.count() / 1.0e6 << " seconds" << endl;
    std::cout<< "Finished tree search "<<std::endl;


    //check tree awnsers against brute force.
    compare_results(results_b, results_t);
}

void compare_results(std::vector<std::vector<size_t>> &v1,
          std::vector<std::vector<size_t>> &v2) {
    for (int i = 0; i< v1.size() ; ++i){
        std::sort(v1[i].begin(), v1[i].end());
        std::sort(v2[i].begin(), v2[i].end());
        std::vector<int> intersection;
        intersection.clear();
        std::set_intersection(v1[i].begin(), v1[i].end(), v2[i].begin(), v2[i].end(),
                              std::back_inserter(intersection));
        if((intersection.size() != v1[i].size()) || (v1[i].size() != v2[i].size())){
            std::cout << "answer diff for i = " << i << std::endl;
        }
    }
}
