#ifndef FREGRID_DITREE_H
#define FREGRID_DITREE_H

#include "BBox3D.h"
#include "BoxedObj.h"
#include "DistanceInterval.h"
#include "DITreeNode.h"
#include "Comparators.h"
#include "Partition.h"
#include "TreePerfStats.h"

using DI=nct::DistanceInterval<float>;

/* The DITree is a version of the DETree presented in "Ray Queries with Wide Object
 * Isolation and the DE-Tree" M. Zuniga and J. Uhlmann, JGT, 2006. It was chosen for
 * use because it was shown to be as optimal as the BSP tree for ray queries, and yet
 * with a faster setup overhead. It can be seen as a form of the bounding volume hierarchy
 * with these features:
 * a) Nodes store one dimensional extents (or bounding intervals)  instead of 3D bounding boxes.
 * b) A node has a partitioning  dimension  - the one which was used to partition the node
 * during tree construction.
 * c) Each node has a bounding interval for each child, and in the partitioning dimension.
 * d) Node partitioning is designed such that a nodes interval do not significantly overlap.
 * In this version, we also have these features:
 * e) partitioning dimensions are cyclically rotated at ech level of the tree (as introduced in
 * "Multidimensional Binary Search Trees Used for Associative Searching", Jon Louis Bentley, 1975)
 * f) Spatial media (or distance mid-range) is used in partitioning. Other partitioning schemes are possible.
 * g) There is a box (or object) associated with each node - not just the leaf nodes. A leaf-oriented
 *  tree is a possible variation which would be slightly simpler and easier to parallelize if
 *  ever needed - but it is not done here.
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
                calculateDI(begin, begin, end, end, dim);
                return *begin;
            }

            // partition all objects into two sets
            float min = 1.0e10;
            float max = -1.0e10;
            for(auto itr = begin; itr != end; itr++){
              auto boxp = (*itr)->obj->getBox();
              if (boxp->getLo(dim) < min) min = boxp->getLo(dim);
              if (boxp->getHi(dim) > max) max = boxp->getHi(dim);
            }
            auto median = partition(begin, end, dim);

            //Make the begin iterator now point to the max elem from the LHS;
            // this one becomes the nodes element
            auto maxElemLHS = std::max_element(begin, median, LessThanBoxHi<Node>(dim)); //TODO:
            std::iter_swap(begin, maxElemLHS); //swaps the objects pointed to by two iterators

            auto firstC = begin + 1;
            calculateDI(begin, firstC, median, end, dim);
            auto nl = median - firstC;
            auto nr = end - median;
            assert(nl + nr +1 == size);
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

            float  min2 = (*minElem)->obj->getBox()->getLo(dim);
            float max2 =  (*maxElem)->obj->getBox()->getHi(dim);

            LessThanValBoxMed<Node> ltm(dim, medDist);//TODO: experiemnt w/ boxHi .LT. med.
            auto median = std::stable_partition(firstC, end, ltm);
            // One way to take care of degenerate case:
            if (median == firstC || median == end) {
                median = firstC + (end - firstC) / 2;
            }
            return median;
        }


        void calculateDI(NodeItr ndItr, const NodeItr begin, const NodeItr median,
                         const NodeItr end, const int dim ){
            if (begin == end) {
                return; //TODO: dont do anything?
            }
            (*(ndItr))->dim = dim;
            calculateDI((*(ndItr))->diL, begin, median, dim);
            calculateDI((*(ndItr))->diR, median, end, dim);
        }
        // calculate one distance interval
        void calculateDI(DI & di ,const NodeItr begin, const NodeItr end, const int dim ){
            if (begin == end) {
                return; //TODO: dont do anything?
            }
            auto minElem = std::min_element(begin, end, LessThanBoxLo<Node>(dim));
            auto min = (*(minElem))->obj->getBox()->getLo(dim);
            auto maxElem = std::max_element(begin, end, LessThanBoxHi<Node>(dim));
            auto max = (*(maxElem))->obj->getBox()->getHi(dim);
            di = DistanceInterval(min, max);
        }

        void search (NodePtr node, BoxAndId& idBox, std::vector<size_t>& results) {
            perfs.incNodesVisited();
            idBox.addResultIf(node->obj, results);
            if(node->isLeaf() ){
                return;
            }
            BBox3D* box = node->obj->getBox();
            if(node->left != nullptr && BBox3D::intersect (*(idBox.getBox()), node->diL,node->dim)) {
                search(node->left, idBox, results);
            }
            if(node->right != nullptr && BBox3D:: intersect ((*idBox.getBox()) , node->diR, node->dim)) {
                search(node->right, idBox, results);
            }
        }
    public:
        TreePerfStats perfs;
        explicit DITree(std::vector < T >& objects, PartType partT = PartType::DMR) {
            size_t nofDIExpected = objects.size();
            nodes.reserve(objects.size());
            for (auto& object : objects) {
                nodes.push_back(new Node(&object));
            }

            std::cout << "DITree buildTree starting" << std::endl;
            if( partitionType == PartType::DMR ){
                 //this->pFunction = getPartitionFunctor<Node,NodeItr,T,M>( this->partitionType );
                root = buildTree(nodes.begin(), nodes.end(), 0);
            }
            perfs.incObjectsCount( objects.size());
            std::cout << "DITree buildTree finished" << std::endl;
        }
        ~DITree() {
            std::cout << "DITree destructor" << std::endl;
            for (const auto &nd: nodes) {
                delete nd;
            }
            nodes.clear();
        }

        void search (BoxAndId& idBox, std::vector<size_t>& results){
          perfs.incQueriesCount();
          search (root, idBox,results) ;
        }

        void search (std::vector<BoxAndId>& idBoxes, std::vector<std::vector<size_t>>& results) {
            for (int i= 0; i< idBoxes.size(); ++i){
                perfs.incQueriesCount();
                search(idBoxes[i], results[i]);
            }
        }
    };

} // nct

#endif //FREGRID_DITREE_H
