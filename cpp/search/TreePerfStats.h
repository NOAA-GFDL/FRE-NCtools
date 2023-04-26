//
// Created by mzuniga on 4/25/23.
//

#ifndef FREGRID_TREEPERFSTATS_H
#define FREGRID_TREEPERFSTATS_H

#include <iostream>
#include <cmath>

class TreePerfStats {
private:
    long long int nodesVisited;
    long long int distanceCalls;
    long long int queriesCount;
    long long int objectsCount;
public:
    TreePerfStats() : nodesVisited{ 0 }, distanceCalls{ 0 }, queriesCount {0}, objectsCount{ 0}{}
    [[nodiscard]] long long int getNodesVisited() const { return nodesVisited; }
    [[nodiscard]] long long int getDistanceCalls() const { return distanceCalls; }
    [[nodiscard]] long long int getQueriesCount() const { return queriesCount; }
    [[nodiscard]] long long int getObjectsCount() const {return objectsCount; }
    void incNodesVisited() { nodesVisited++; }
    void incDistanceCalls() { distanceCalls++; }
    void incQueriesCount(){ queriesCount ++; }
    void incNodesVisited(long long int v) { nodesVisited += v;  }
    void incDistanceCalls(long long int v) { distanceCalls += v; }
    void incQueriesCount(long long int v) { queriesCount += v; }
    void incObjectsCount(long long int v) { objectsCount += v; }

    void addStats( TreePerfStats & s) {
      nodesVisited += s.getNodesVisited();
      distanceCalls += s.getDistanceCalls();
      queriesCount += s.getQueriesCount();
      objectsCount += s.getObjectsCount();
    }
    void reset() {
      nodesVisited = 0;
      distanceCalls = 0;
      queriesCount = 0;
      objectsCount = 0;
    }

    friend std::ostream &operator<< (std::ostream & os, const TreePerfStats& ps){
        auto log2size = log2( ps.objectsCount);
        os << "# objects, # queries , #nodes visited, <# nodes visited>, (<#nodes visited>/log2(# obj):"
            << std::endl << ps.objectsCount << ", " << ps.queriesCount << ", " << ps.nodesVisited
            << ", " <<  ps.nodesVisited / ps.queriesCount << ", " <<
              ps.nodesVisited / (log2size * ps.queriesCount );
      return os;
    }

};

#endif //FREGRID_TREEPERFSTATS_H
