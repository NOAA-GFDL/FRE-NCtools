#ifndef FREGRID_BRUTEBOXQUERY_H
#define FREGRID_BRUTEBOXQUERY_H

#include <vector>
#include <iostream>
#include <cmath>
#include <CL/sycl.hpp>
#include "BBox3D.h"
#include "BoxedObj.h"

namespace nct{
using  std::vector;
using nct::BBox3D;

/**
 *  A box query class whos search method is the simple brute (exhaustive) search.
 * The boxes provided in the constructor are the boxes which will later be searched
 * by using the search function. The member bpairs should have a one-to-one correspondence
 * with member boxes Its up to the user to determine what the ids of bpairs are, but
 * usually its just the id index into bpairs.
 **/
class BruteBoxQuery {
    private:
        vector<BoxAndId> & bpairs;
        vector<BBox3D> & boxes;
    public:
    BruteBoxQuery(vector <BoxAndId> &bpairs, vector<BBox3D> &boxes)
    : bpairs(bpairs), boxes(boxes) {}

    void search(vector<BBox3D>  & qboxes,  vector<vector<size_t>> & results) {
            for (auto iq = 0; iq<qboxes.size(); iq++){
                //TODO: experiment with iterating directly over bpairs instead.
                for(auto jb = 0; jb < boxes.size(); jb++){
                    if(BBox3D::intersect(qboxes[iq], boxes[jb])) {
                        //atomic here ?
                        results[iq].push_back(bpairs[jb].id);
                    }
            }
        }
    }

    void
    search_sycl(vector<BBox3D>  & qboxes,  vector<vector<size_t>> & results) {
      using namespace cl::sycl;
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


    void
    search_sycl(vector<BBox3D>  & qboxes,  vector<size_t>  & results,
                vector<size_t>  & res_count ,size_t max_res) {
      using namespace cl::sycl;

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


      queue q;
      std::cout << "Using device: "
                << q.get_device().get_info<sycl::info::device::name>()
                << std::endl;;
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



    void search(BBox3D & qBox, vector<size_t> & results) {
        //TODO: experiment with iterating directly over bpairs instead.
        for(auto jb = 0; jb < boxes.size(); jb++){
          if(BBox3D::intersect(qBox, boxes[jb])) {
            //atomic here ?
            results.push_back(bpairs[jb].id);
          }

      }
    }

    //c is row major order
    inline static size_t qr_idx( size_t ir, size_t iq, size_t RSIZE){
      return (iq * RSIZE + ir);
    }
};
} // nct
#endif //FREGRID_BRUTEBOXQUERY_H
