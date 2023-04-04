#ifndef FREGRID_BOXEDOBJ_H
#define FREGRID_BOXEDOBJ_H

#include "BBox3D.h"
#include "BruteBoxQuery.h"

//This file defines the objects we can put in the DITree.
//Version 1.0: We will only put BoxAndId objects, and the near neighbors search results
// will be vectors of ids (so the objects will need to be looked up elsewhere,
// such as in an oarray of cells

namespace nct {
    using ResultIds = std::vector<unsigned int>;
    class BoxAndId {
        friend class BruteBoxQuery;
    private:
        unsigned int id;
        BBox3D *box;
    public:
        BoxAndId(unsigned int id, BBox3D *box) : id(id), box(box) {}
        BBox3D *getBox() {
            return box;
        }
        unsigned int getId() {
            return id;
        }
    };
}
#endif //FREGRID_BOXEDOBJ_H
