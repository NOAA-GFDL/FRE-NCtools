#ifndef FREGRID_BRUTEBOXQUERY_H
#define FREGRID_BRUTEBOXQUERY_H

#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <ranges>
#include <execution>
#include <numeric>
#include <cmath>
#include <string>

#include <mutex>
#include <thread>

#ifdef USE_SYCL
#include <CL/sycl.hpp>
#endif
//#include <cartesian_product.hpp>

#include "BBox3D.h"
#include "BoxedObj.h"

constexpr int MAX_NN=10;
constexpr int MAX_NNP = MAX_NN + 1;

namespace nct {

using std::array;
using std::string;
using std::vector;
constexpr size_t fill_value_size_t = std::numeric_limits<size_t>::max();

/**
 *  A box query class whos search method is the simple brute (exhaustive) search.
 * The boxes provided in the constructor are the boxes which will later be searched
 * by using the search function. The member bpairs should have a one-to-one correspondence
 * with member boxes Its up to the user to determine what the ids of bpairs are, but
 * usually its just the id index into bpairs.
 **/
class BruteBoxQuery {
 private:
  vector<BoxAndId> &bpairs;
  vector<BBox3D> &boxes;
 public:
  enum class search_subtype
   { simple = 0, parittion, for_each, copy_if,sycl };
   constexpr static std::array  search_subtype_str = 
   {"simple", "parittion", "for_each", "copy_if", "sycl" };

  BruteBoxQuery(vector<BoxAndId> &bpairs, vector<BBox3D> &boxes)
      : bpairs(bpairs), boxes(boxes) {}

  void search(BBox3D &qBox, vector<size_t> &results) {
    // TODO: experiment with iterating directly over bpairs instead.
    for (size_t jb = 0; jb < boxes.size(); jb++) {
      if (BBox3D::intersect(qBox, boxes[jb])) {
        // atomic here ?
        results.push_back(bpairs[jb].id);
      }
    }
  }

  void search(vector<BBox3D> &qboxes, vector<vector<size_t>> &results, search_subtype sstype) {
    switch (sstype) {
      case search_subtype::simple:
        search_simple(qboxes, results);
        break;
      case search_subtype::parittion:
        for (size_t i = 0; i < qboxes.size(); i++) {
          search_std_partition(qboxes[i], results[i]);
        }
        break;
      case search_subtype::for_each:
        for (size_t i = 0; i < qboxes.size(); i++) {
          search_std_for_each(qboxes[i], results[i]);
        }
        break;
      case search_subtype::copy_if:
        search_copy_if(qboxes, results);
        break;
      default:
        std::cout << "Error invalid bute force search subtype" << std::endl;
    }
  }

  void search_simple(vector<BBox3D> &qboxes, vector<vector<size_t>> &results) {
    for (size_t iq = 0; iq < qboxes.size(); iq++) {
      // TODO: experiment with iterating directly over bpairs instead.
      for (size_t jb = 0; jb < boxes.size(); jb++) {
        if (BBox3D::intersect(qboxes[iq], boxes[jb])) {
          // atomic here ?
          results[iq].push_back(bpairs[jb].id);
        }
      }
    }
  }


void search_std_partition(BBox3D &qBox, vector<size_t> &results) {
  // TODO: make copy of bpairs ( vector<BoxAndId> bps) to hav one the perserves
  //  origal order
  auto it = std::partition(
      std::execution::par,
      bpairs.begin(), bpairs.end(),
      [&qBox](BoxAndId bp) { return BBox3D::intersect(*(bp.getBox()), qBox);
  });

  results.clear();
  while (it != bpairs.end()) {
     results.push_back((*it).getId());
     it++;
  }
}

void search_copy_if(vector<BBox3D> &qboxes, vector<vector<size_t>> &results) {
  auto ints = std::views::iota((int)0, (int)boxes.size());
  std::cout << "capacity and size: " << results[0].capacity() << ", " << results[0].size() << std::endl;
  for (int i = 0; i < (int)qboxes.size(); i++) {
     auto nn_iter = results[i].begin();
     std::copy_if(
         std::execution::par,
         // NOTE: nvc++ w/ GPU requires random_access iterators for 3r darg  (result) .
         ints.begin(), ints.end(), nn_iter,
         [=, qboxes = qboxes.data(), boxes = boxes.data()](int j) {
           return BBox3D::intersect(qboxes[i], boxes[j]);
         });
  }

 std::cout << "capacity and size: " << results[0].capacity() << ", " << results[0].size() << std::endl;;
  // remove the fill values at the end of the results vector.
  for (int i = 0; i < (int)qboxes.size(); i++) {
     while ((results[i].size() > 0) &&
            (results[i].end()[-1] == fill_value_size_t)) {
        results[i].pop_back();
     }
  }
}

void search_std_for_each(BBox3D &qBox, vector<size_t> &results) {
  /*
  array<int, MAX_NNP> nns;
  nns[MAX_NN] = 0;
  auto ints = std::views::iota(0, (int)bpairs.size());
  std::mutex r_mutex;
  std::for_each_n(
      std::execution::par,
      ints.begin(), bpairs.size(), [&](int i) {
        if (BBox3D::intersect(boxes[i], qBox)) {
          std::lock_guard<std::mutex> guard(r_mutex);
          assert(nns[MAX_NN] < MAX_NN);
          nns[nns[MAX_NN]] = i;
          ++nns[MAX_NN];
        }
      });

  results.clear();
  for (int i = 0; i < nns[MAX_NN]; i++) {
     results.push_back(nns[i]);
  }
  */
}

#ifdef USE_SYCL
    void
    search_sycl(vector<BBox3D>  & qboxes,  vector<vector<size_t>> & results) {
      using namespace sycl;
      queue q;
      for(auto & elem : results){
        elem.reserve( 10);
        elem[0] = 0;
      }
      buffer dboxes_buf(boxes);
      buffer qboxes_buf(qboxes);
      buffer bpairs_buf(bpairs);
      buffer results_buf(results);
      size_t NQB = qboxes.size();
      size_t NDB = boxes.size();
      q.submit([&](handler &h) {
          accessor a_dboxes(dboxes_buf, h, read_only);
          accessor a_qboxes(qboxes_buf, h, read_only);
          accessor a_bpairs(bpairs_buf, h, read_only);
          accessor a_results(results_buf, h, write_only);

          h.parallel_for(NQB, [=](auto iq) {
              for (auto jb = 0; jb < NDB; jb++) {
                if (BBox3D::intersect(a_qboxes[iq], a_dboxes[jb])) {
                  //a_results[iq].push_back(a_bpairs[jb].id); //sycl exceptions not allowed
                  //need to pre-allocate:
                  //auto no = a_results[iq].size();
                 // a_results[iq][no] = a_bpairs[jb].id;
                  a_results[iq][0] += 1;
                }
              } ; });//end of inner lambda
      }); //end of q.submit lambda
    }

 static int custom_device_selector(const sycl::device &d) {
       using namespace cl::sycl;
      int rating = 0;
      if (d.is_gpu() && (d.get_info<info::device::name>().find("Nvidia") != std::string::npos)) {
        rating = 3;
      } else if (d.is_gpu()) {
        rating = 2;
      } else if (d.is_cpu()) {
        rating = 1;
      }
      return rating;
   }

  static void list_devices(){
    for (auto platform : sycl::platform::get_platforms())
    {
        std::cout << "Platform: "
                  << platform.get_info<sycl::info::platform::name>()
                  << std::endl;

        for (auto device : platform.get_devices())
        {
            std::cout << "\tDevice: "
                      << device.get_info<sycl::info::device::name>()
                      << std::endl;
        }
    }
  }

    void
    search_sycl(vector<BBox3D>  & qboxes,  vector<size_t>  & results,
                vector<size_t>  & res_count ,size_t max_res) {
      using namespace sycl;

      size_t NQB = qboxes.size();
      size_t NDB = boxes.size();

      results.reserve(NQB * max_res);
      res_count.reserve( NQB );

      for(int i=0; i< NQB * max_res; ++i){
        results.push_back(-1);
      }
      for(int i=0; i< NQB; ++i){
        res_count.push_back(0);
      }

      list_devices();

      sycl::device preferred_device { custom_device_selector };

      //queue q;
      queue q(preferred_device);
      std::cout << "Using device: "
                << q.get_device().get_info<sycl::info::device::name>()
                << std::endl;

      buffer dboxes_buf(boxes);
      buffer qboxes_buf(qboxes);
      buffer bpairs_buf(bpairs);
      buffer results_buf(results);
      buffer res_count_buf( res_count);

      q.submit([&](handler &h) {
          accessor a_dboxes(dboxes_buf, h, read_only);
          accessor a_qboxes(qboxes_buf, h, read_only);
          accessor a_bpairs(bpairs_buf, h, read_only);
          accessor a_results(results_buf, h, write_only);
          accessor a_res_count(res_count_buf, h, write_only);

          h.parallel_for(NQB, [=](auto iq) {
              for (auto jb = 0; jb < NDB; jb++) {
                if (BBox3D::intersect(a_qboxes[iq], a_dboxes[jb])) {
                  auto iq_rn = qr_idx( a_res_count[iq], iq, 10);
                  a_results[ iq_rn ] = a_bpairs[jb].id;
                  a_res_count[iq] += 1;
                }
              } ; });//end of inner lambda
      }); //end of q.submit lambda
    }
#endif




    //c is row major order
    inline static size_t qr_idx( size_t ir, size_t iq, size_t RSIZE){
      return (iq * RSIZE + ir);
    }
  };
  }    // namespace nct
#endif //FREGRID_BRUTEBOXQUERY_H
