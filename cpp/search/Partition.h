#ifndef FREGRID_PARTITION_H
#define FREGRID_PARTITION_H

namespace nct {

}
#include <algorithm>
#include <vector>

#include "DistanceInterval.h"
#include "DITreeNode.h"
#include "Comparators.h"

    namespace nct{
    enum class PartType {OM, BOM, DMR };

// Abstract base class
    template < typename N, typename  NodeItr>
    class PartitionFunction {
        using DI = nct::DistanceInterval<float>;
    public:
        PartitionFunction() = default;;
        virtual void operator() (const NodeItr firstC, NodeItr& median, const NodeItr end) = 0;
    };

// Add two doubles
    template < typename N, typename  NodeItr>
    class BOM : public PartitionFunction<N,NodeItr> {
    public:
        BOM() {};
        virtual NodeItr operator() (const NodeItr firstC, NodeItr& median, const NodeItr end, int dim) {
            std::nth_element(firstC, median, end, LessThanBoxMed<N>(dim));
            return;
        }
    };

    template < typename N, typename  NodeItr, typename T>
    class DMR : public PartitionFunction<N,NodeItr>{
    public:
        DMR(){};
        NodeItr operator()(const NodeItr firstC, const NodeItr end, int dim){
            auto minElem = std::min_element(firstC, end, LessThanBoxLo<N>(dim));
            auto maxElem = std::max_element(firstC, end, LessThanBoxHi<N>(dim));
            double medDist = 0.5 * ((*maxElem)->obj->getBox()->getLo(dim) + (*minElem)->obj->getBox()->getHi(dim));

            LessThanBoxMed<N> ltm(dim, medDist);
            auto median = std::stable_partition(firstC, end, ltm);
            // One way to take care of degenerate case:
            if (median == firstC || median == end){
                median = firstC + (end - firstC) / 2;
            }
            return median;
        }
    };



} // nct

#endif //FREGRID_PARTITION_H
