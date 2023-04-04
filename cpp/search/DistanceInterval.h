#ifndef FREGRID_DISTANCEINTERVAL_H
#define FREGRID_DISTANCEINTERVAL_H

namespace nct {

    template <typename T>
    class DistanceInterval {
        private:
            T d[2]; //Also possible: T near, far;
        public:
            DistanceInterval(const T nd, const T fd) {//Note std::nextafter not used in contructor.
                d[0] = nd; d[1] = fd;
            }
            DistanceInterval() {
                d[0] = 0; d[1] = 0;
            }
            inline void setNear(const T nd) {
#ifndef USE_NEXTAFTER
                d[0] = nd;
#else
                d[0] = std::nextafter(nd,  0.0f);
#endif
            }
            inline void setFar(const T nd) {
#ifndef USE_NEXTAFTER
                d[1] = nd;
#else
                d[1] = std::nextafter(nd,  std::numeric_limits<float>::max());
#endif
            }
            inline const T getNear() { return d[0]; }
            inline const T getFar() { return d[1]; }
            inline void reset() { d[0] = 0; d[1] = 0; }
            DistanceInterval& operator = (const DistanceInterval& that) {
                d[0] = that.d[0];
                d[1] = that.d[1];
                return *this;
            }
            inline bool rangeOverlaps(const double pqDistance, const double radius) {
                if ((pqDistance >= getNear() - radius) &&
                    (pqDistance <= getFar() + radius)) {
                    return true;
                } else {
                    return false;
                }
            }
        };
} // nct

#endif //FREGRID_DISTANCEINTERVAL_H
