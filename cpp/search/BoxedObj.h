#ifndef FREGRID_BOXEDOBJ_H
#define FREGRID_BOXEDOBJ_H

#include "BBox3D.h"
#include "BruteBoxQuery.h"

//This file defines the objects we can put in the DITree.
//Version 1.0: The DITree will contain (the nodes point to) only
// objects with the same interface as BoxAndId objects. In this version the
// results are ids and stored in a separate container. See the search tree and
// function addResult

namespace nct {
    using ResultIds = std::vector<size_t>;
    class BoxAndId {
        friend class BruteBoxQuery;
    private:
        size_t id;
        BBox3D *box;
    public:
        BoxAndId(size_t id, BBox3D*box ) : id(id), box(box) {}
        BBox3D *getBox() {
            return box;
        }
        size_t getId() const {
            return id;
        }
        /*
         * Check the candidate object cn and add it (add its id) to the results if it intersects
         * this objects box.
         */
        void addResultIf(BoxAndId* cn, std::vector<size_t>& results){
            if(BBox3D::intersect( *(this->getBox()), *(cn->getBox()))){
                results.push_back( cn->getId());
            }
        }
    };
}
#endif //FREGRID_BOXEDOBJ_H
