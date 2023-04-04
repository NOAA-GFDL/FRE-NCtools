#ifndef FREGRID_DITREENODE_H
#define FREGRID_DITREENODE_H

#include "DistanceInterval.h"

namespace nct {
    template <class T>
    class DITreeNode {
        public:
            BoxAndId* obj; // Pointer to the nodes object. Objects required to have bounding boxes
            int size;
            int dim; //Component (dimension) of the bounding boxes used in calculation of this node.
            DistanceInterval<float> diL; //Bounding distance interval of Left children of node in the nodes dim
            DistanceInterval<float> diR;  //Bounding distance interval of Right children of node in the nodes dim
            DITreeNode* left;  //left child.
            DITreeNode* right; //right child
            explicit DITreeNode(T* object) : obj{ object } {}
            bool isLeaf() { return ((left == nullptr) && (right == nullptr)); }
            void setLeaf() { left = nullptr;  right = nullptr; }//TODO ; set interval?
        };
} // nct

#endif //FREGRID_DITREENODE_H
