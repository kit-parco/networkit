/*
 * UnionFind.cpp
 *
 *  Created on: 11.08.2014
 *      Author: Marcel Radermacher
 *      Changed a bit by Henning Meyerhenke to reflect union by rank and path compression
 *        as taught in "Algorithms 1"
 */

#include <networkit/structures/UnionFind.hpp>

namespace NetworKit {


void UnionFind::allToSingletons() {
    for (index i = 0; i < parent.size(); ++i) {
        parent[i] = i;
    }
}

index UnionFind::find(index u) {
    if (parent[u] == u) {
        return u;
    }
    else {
        // recursion and path compression
        parent[u] = find(parent[u]);
        return parent[u];
    }
}

void UnionFind::merge(index u, index v) {
    index set_u = find(u);
    index set_v = find(v);
    if (set_u == set_v) return;

    if (rank[set_u] < rank[set_v]) {
        parent[set_u] = set_v;
    }
    else {
        parent[set_v] = set_u;
        if (rank[set_u] == rank[set_v]) {
            rank[set_u]++;
        }
    }
}

Partition UnionFind::toPartition() {
    Partition p(parent.size());
    p.setUpperBound(parent.size());
    for (index e = 0; e < parent.size(); ++e) {
        p.addToSubset(find(e), e);
    }
    return p;
}


}
