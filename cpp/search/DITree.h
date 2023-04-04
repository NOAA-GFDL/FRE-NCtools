#ifndef FREGRID_DITREE_H
#define FREGRID_DITREE_H

#include "BBox3D.h"
#include "BoxedObj.h"
#include "DistanceInterval.h"
#include "DITreeNode.h"
#include "Comparators.h"
#include "Partition.h"

using DI=nct::DistanceInterval<float>;

/* The DITree is a version of the bounding volume hierarchy where the nodes store
 * one dimensional bounding (interval) information of instead of 3D bounding box
 * information. The dimension of the node is the partitioning dimention - the one
 * which was used to partition the node during tree construction, which partitioned
 * (or divided) the elements under that node into two sets - one for the left
 * child and one for the right child, that are not significantly overlapping in that
 * dimension. In this first version:
 *  1. the bounding intervals of both the left and child nodes are stored at each node.
 *  2. each node (both internal and leaf) has a reference to one object (with a bounding
 *  volume)
 * */

namespace nct {
    template<typename T>
    class DITree {
        using Node = DITreeNode<T>;
        using NodePtr = Node *;
        using Box = BBox3D;
        using BoxPtr = Box *;
        using NodeItr = typename std::vector<NodePtr>::iterator;
    private:
        Node* root;
        std::vector < Node*> nodes;
        PartType  partitionType{ PartType::DMR };
    protected:
        Node* buildTree(const NodeItr begin, const NodeItr end, int dim) {
            if (begin == end) {
                return nullptr;
            }
            int size = end - begin;

            if (size == 1) {
                (*begin)->setLeaf();
                return *begin;
            }

            // partition all objects into two sets
            auto median = partition(begin, end, dim);

            //Make the begin iterator now point to the max elem from the LHS;
            // this one becomes the nodes element
            auto maxElemLHS = std::max_element(begin, median, LessThanBoxLo<Node>(dim));
            std::iter_swap(begin, maxElemLHS); //swaps the objects pointed to by two iterators

            auto firstC = begin + 1;
            calculateDI(begin, firstC, median, end, dim);
            (*begin)->left = this->buildTree(firstC, median, (dim + 1) % 3);
            (*begin)->right = this->buildTree(median, end, (dim + 1) % 3);

            return *begin;
        }

        //DMR partitioning
        auto partition(const NodeItr firstC, const NodeItr end,  unsigned int dim) {
            auto minElem = std::min_element(firstC, end, LessThanBoxLo<Node>(dim));
            auto maxElem = std::max_element(firstC, end, LessThanBoxHi<Node>(dim));
            double medDist = 0.5 * ((*minElem)->obj->getBox()->getLo(dim) +
                    (*maxElem)->obj->getBox()->getHi(dim));

            LessThanValBoxMed<Node> ltm(dim, medDist);
            auto median = std::stable_partition(firstC, end, ltm);
            // One way to take care of degenerate case:
            if (median == firstC || median == end) {
                median = firstC + (end - firstC) / 2;
            }
            return median;
        }


        void calculateDI(NodeItr ndItr, const NodeItr begin, const NodeItr median,
                         const NodeItr end, int dim ){
            if (begin == end) {
                return; //TODO: dont do anything?
            }
            calculateDI((*(ndItr))->diL, begin, median, dim);
            calculateDI((*(ndItr))->diR, median, end, dim);
        }
        // calculate one distance interval
        void calculateDI(DI & di ,const NodeItr begin, const NodeItr end, int dim ){
            if (begin == end) {
                return; //TODO: dont do anything?
            }
            auto minElem = std::min_element(begin, end, LessThanBoxLo<Node>(dim));
            auto min = (*(minElem))->obj->getBox()->getLo(dim);
            auto maxElem = std::max_element(begin, end, LessThanBoxHi<Node>(dim));
            auto max = (*(minElem))->obj->getBox()->getHi(dim);
            di = DistanceInterval(min, max);
        }

        void search (NodePtr node, BoxAndId& idBox, std::vector<unsigned int>& results) {
            processObj(node->obj, idBox, results);
            if(node->isLeaf() ){
                return;
            }
            BBox3D* box = node->obj->getBox();
            if(node->left != nullptr && intresects (node->diL, idBox.getBox() ,node->dim)) {
                search(node->left, idBox, results);
            }
            if(node->right != nullptr && intresects (node->diR, idBox.getBox() ,node->dim)) {
                search(node->right, idBox, results);
            }
        }
    public:
        DITree(std::vector < T >& objects, PartType partT = PartType::DMR) {
            unsigned int nofDIExpected = objects.size();
            nodes.reserve(objects.size());
            for (auto& object : objects) {
                nodes.push_back(new Node(&object));
            }

            std::cout << "DITree buildTree starting" << std::endl;
            if( partitionType == PartType::DMR ){
                 //this->pFunction = getPartitionFunctor<Node,NodeItr,T,M>( this->partitionType );
                root = buildTree(nodes.begin(), nodes.end(), 0);
            }
            std::cout << "DITree buildTree finished" << std::endl;
        }
        ~DITree() {
            std::cout << "DITree destructor" << std::endl;
            for (const auto &nd: nodes) {
                delete nd;
            }
            nodes.clear();
            //delete pFunction;
        }

        void search (BoxAndId& idBox, std::vector<unsigned int>& results){
          search (root, idBox,results) ;
        }
    };

} // nct

#endif //FREGRID_DITREE_H
