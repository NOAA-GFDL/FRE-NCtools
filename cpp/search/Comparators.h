#ifndef FREGRID_COMPARATORS_H
#define FREGRID_COMPARATORS_H

/*
 * Includes comparators and partition functions for nodes with
 * boxed objects.
 */

#include "BBox3D.h"
namespace nct {

    template <class Nd>
    class LessThanBoxLo {
    private:
        unsigned int dim;
    public:
        explicit LessThanBoxLo(unsigned int d ): dim{d}{ }
        bool operator()(const Nd* left, const Nd* right) const {
            return (left->obj->getBox()->getLo(dim) < right->obj->getBox()->getLo(dim)) ;
        }
    };

    template <class Nd>
    class LessThanBoxHi {
    private:
        unsigned int dim;
    public:
        explicit LessThanBoxHi(unsigned int d ): dim{d}{ }
        bool operator()(const Nd* left, const Nd* right) const {
            return (left->obj->getBox()->getHi(dim) < right->obj->getBox()->getHi(dim)) ;
        }
    };

    template <class Nd>
    class LessThanBoxMed {
    private:
        unsigned int dim;
    public:
        explicit LessThanBoxMed(unsigned int d): dim{d}{ }
        //return true iff the LHS box < RHS box median in the given dimension.
        //In this case the box is actually the box of the object pointed to by the node
        bool operator()(const Nd* left, const Nd* right) const {
            return ((left->obj->getBox()->getLo(dim)  + left->obj->getBox()->getHi(dim)) <
                    (right->obj->getBox()->getLo(dim) + right->obj->getBox()->getHi(dim)));
        }
    };

    template <class Nd>
    class LessThanValBoxMed{
    private:
        unsigned int dim;
         float val;
    public:
        explicit LessThanValBoxMed(unsigned int d, float v): dim{d},val{v}{}
        //return true iff the LHS box < RHS box median in the given dimension.
        //In this case the box is actually the box of the object pointed to by the node
        bool operator()(const Nd* nd) const {
            return ((nd->obj->getBox()->getLo(dim) + nd->obj->getBox()->getHi(dim)) < (2.0 *val));
        }
    };
} // nct

#endif //FREGRID_COMPARATORS_H
