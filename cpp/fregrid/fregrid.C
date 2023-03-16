#include <iostream>
#include <vector>
#include <array>
#include "BBox3D.h"

using BoxVec_t  = std::vector<nct::BBox3D>;
using BPairVec_t  = std::vector<nct::BoxPair>;

void  initBoxesAndPairs(BoxVec_t & boxes, BPairVec_t & bps);

int main() {
    using namespace std;
    using namespace nct;

    using BoxVec_t  = std::vector<BBox3D>;

    size_t npoly {100};

    cout << "Hello!" << endl;

    vector<BoxPair> bpairs;
    BoxVec_t boxes;

    bpairs.reserve(npoly);
    boxes.reserve(npoly);
    initBoxesAndPairs(boxes, bpairs);

    cout << "Goodbye!" << endl;
    return 0;
}

void initBoxesAndPairs(BoxVec_t & boxes, BPairVec_t & bps){
    int i = 0;
    std::array<float, 3> lo;
    std::array<float, 3> hi;
    for (size_t j = 0; j<boxes.capacity(); j++){
        for(auto & item : lo) item=j;
        for(auto & item : hi ) item = j+1;
        boxes.push_back(nct::BBox3D{lo,hi});
    }
    for (auto& box : boxes){
        bps.push_back(nct::BoxPair{i,box});
        ++i;
    }
    std::cout << "initBoxesAndPairs i="<<i<<std::endl;
}
