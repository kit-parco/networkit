/*
 * Cover.cpp
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#include <networkit/structures/Cover.hpp>

#include <algorithm>
#include <iterator>

namespace NetworKit {

Cover::Cover() : z(0), omega(0), data(0) {}

Cover::Cover(index z) : z(z-1), omega(0), data(z) {
}

Cover::Cover(const Partition &p) : z(p.numberOfElements()-1), omega(p.upperBound()-1), data(p.numberOfElements()) {
    p.forEntries([&](index e, index s) {
        if (s != none)
            data[e].insert(s);
    });
}

bool Cover::contains(index e) const {
    return (e <= z) && (! data[e].empty());	// e is in the element index range and the entry is not empty
}

bool Cover::inSameSubset(index e1, index e2) const {
    assert (e1 <= z);
    assert (e2 <= z);
    assert (! data[e1].empty());
    assert (! data[e2].empty()); // elements cannot be unassigned - it may be possible to change this behavior
    std::set<index> intersect;
    std::set_intersection(data[e1].begin(),data[e1].end(),data[e2].begin(), data[e2].end(), std::inserter(intersect,intersect.begin()));
    return (!intersect.empty());
}

std::set<index> Cover::getMembers(const index s) const {
    assert (s <= omega);
    std::set<index> members;
    for (index e = 0; e <= this->z; ++e) {
        for (index t : data[e]) {
            if (t == s) {
                members.insert(e);
            }
        }
    }
    return members;
}

void Cover::addToSubset(index s, index e) {
    assert (e <= z);
    assert (s <= omega);
    data[e].insert(s);
}

void Cover::removeFromSubset(index s, index e) {
    assert (e <= z);
    assert (s <= omega);
    data[e].erase(s);
}


void Cover::moveToSubset(index s, index e) {
    assert (e <= z);
    assert (s <= omega);
    data[e].clear();
    data[e].insert(s);
}

index Cover::toSingleton(index e) {
    assert (e <= z);
    data[e].clear();
    index sid = newSubsetId();
    data[e].insert(sid);
    return sid;
}

void Cover::allToSingletons() {
    for (index e = 0; e <= this->z; ++e) {
        toSingleton(e);
    }
}

void Cover::mergeSubsets(index s, index t) {
    assert (s <= omega);
    assert (t <= omega);
    if ( s != t ) {
        index m = newSubsetId(); // new id for merged set
        for (index e = 0; e <= this->z; ++e) {
            auto its = data[e].find(s);
            auto itt = data[e].find(t);
            if (its != data[e].end()) {
                data[e].erase(its);
                data[e].insert(m);
            }
            // was else if. makes errors, in case an element is in s as well as t
            if (itt != data[e].end()) {
                data[e].erase(itt);
                data[e].insert(m);
            }
        }
    }
}

index Cover::upperBound() const {
    return omega + 1;	// to enable usual loop test x < upperBound()
}

index Cover::lowerBound() const {
    return 0;
}

std::vector<count> Cover::subsetSizes() const {
    std::map<index,count> mapping;
    std::vector<count> sizes;
    count newIndex = 0;
    for (index e = 0; e <= this->z; ++e) { // stores sizes in a vector
        for (index t : data[e]) {
            if (mapping.find(t) == mapping.end()) {
                mapping[t] = newIndex++;
                sizes.push_back(1);
            } else {
                sizes[mapping[t]]++;
            }
        }

    }
    return sizes;
}

std::map<index, count> Cover::subsetSizeMap() const {
    std::map<index,count> sizeMap;
    for (index e = 0; e <= this->z; ++e) { // stores sizes of subsets in a map
        for (index t : data[e]) {
            if (sizeMap.find(t) == sizeMap.end()) {
                sizeMap[t] = 1;
            } else {
                sizeMap[t]++;
            }
        }
    }
    return sizeMap;
}

count Cover::numberOfSubsets() const {
    std::vector<int> exists(upperBound(), 0); // a boolean vector would not be thread-safe

    this->parallelForEntries([&](index, std::set<index> s) {
        if (!s.empty()) {
            for (auto it = s.begin(); it != s.end(); it++) {
                index currentSubset = *it;
                exists[currentSubset] = 1;
            }
        }
    });

    count k = 0; // number of actually existing clusters
    #pragma omp parallel for reduction(+:k)
    for (omp_index i = 0; i < static_cast<omp_index>(upperBound()); ++i) {
        if (exists[i]) {
            k++;
        }
    }

    return k;
}

count Cover::numberOfElements() const {
    return z+1;
}

index Cover::extend() {
    data.emplace_back();
    ++z;
    assert(z + 1 == data.size());
    return z;
}

void Cover::setUpperBound(index upper) {
    this->omega = upper -1;
}

std::set<index> Cover::getSubsetIds() const {
    std::set<index> ids;
    for (std::set<index> subset : data) {
        ids.insert(subset.begin(), subset.end());
    }
    return ids;
}

std::vector<std::set<index>> Cover::getSubsets() const {
    auto ids = getSubsetIds();
    auto idMapping = std::vector<index>(upperBound());
    std::vector<std::set<index>> subsets(ids.size());
    index newSubsetId = 0;
    for (index id : ids) {
        idMapping[id] = newSubsetId++;
    }

    for (index element = 0; element <= this->z; ++element) {
        for (index subsetId : data[element]) {
            subsets[idMapping[subsetId]].insert(element);
        }
    }
    return subsets;
}

void Cover::addSubset(const std::set<index> &subset) {
    index id = newSubsetId();
    for (node e : subset)
        addToSubset(id, e);
}

} /* namespace NetworKit */
