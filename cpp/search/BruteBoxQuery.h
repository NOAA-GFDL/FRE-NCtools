#ifndef FREGRID_BRUTEBOXQUERY_H
#define FREGRID_BRUTEBOXQUERY_H

#include <vector>
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
    void search(BBox3D & qBox, vector<size_t> & results) {
        //TODO: experiment with iterating directly over bpairs instead.
        for(auto jb = 0; jb < boxes.size(); jb++){
          if(BBox3D::intersect(qBox, boxes[jb])) {
            //atomic here ?
            results.push_back(bpairs[jb].id);
          }

      }
    }
};
} // nct
#endif //FREGRID_BRUTEBOXQUERY_H
