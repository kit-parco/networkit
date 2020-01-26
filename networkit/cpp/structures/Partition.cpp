/*
 * Partition.cpp
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#include <algorithm>
#include <numeric>

#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/structures/Partition.hpp>
#include <tlx/define/likely.hpp>

namespace NetworKit {

Partition::Partition() : z(0), omega(0), data(0) {

}


Partition::Partition(const std::vector<index>& data) : z(data.size()), omega(), data(data) {
    auto max_elem = *std::max_element(data.begin(), data.end());
    this->omega = (max_elem == none) ? 0 : max_elem;
}


Partition::Partition(index z) : z(z), omega(0), data(z, none) {  //z(z-1);data(z,none);

}

Partition::Partition(index z, index defaultValue) : z(z), omega(0), data(z, defaultValue) {  //z(z-1);data(z,none);

}

void Partition::allToSingletons() {
    setUpperBound(numberOfElements());
    std::iota(data.begin(), data.end(), 0);
}

index Partition::mergeSubsets(index s, index t) {
    assert (s <= omega);
    assert (t <= omega);
    if (s != t) {
        index m = newSubsetId(); // new id for merged set
        for (index e = 0; e < this->z; ++e) {
            if (data[e] == s || data[e] == t) {
                data[e] = m;
            }
        }
        return m;
    }
    return none; // no new cluster formed
}

count Partition::numberOfSubsets() const {
    auto n = upperBound();
    std::vector<bool> exists(n);
    count k = 0; // number of actually existing clusters
    this->forEntries([&](index, index s) {
        if (s != none) {
            assert(s < n);
            if (!exists[s]) {
                ++k;
                exists[s] = true;
            }
        }
    });
    return k;
}

void Partition::compact(bool useTurbo) {
    index i = 0;
    if (!useTurbo) {
        std::vector<index> usedIds(data);
        Aux::Parallel::sort(usedIds.begin(), usedIds.end());
        auto last = std::unique(usedIds.begin(), usedIds.end());
        usedIds.erase(last, usedIds.end());
        i = usedIds.size();

        this->forEntries([&](index e, index s){ // replace old SubsetIDs with the new IDs
            if (s != none) {
                data[e] = std::distance(usedIds.begin(), std::lower_bound(usedIds.begin(), usedIds.end(), s));
            }
        });
    } else {
        std::vector<index> compactingMap(this->upperBound(), none);
        for (index e = 0; e < z; ++e) {
            const index cid = data[e];
            if (TLX_LIKELY(cid != none)) {
                if (compactingMap[cid] == none) {
                    compactingMap[cid] = i++;
                }
                data[e] = compactingMap[cid];
            }
        }
    }
    this->setUpperBound(i);
}

std::vector<count> Partition::subsetSizes() const {
    std::vector<count> sizes;
    std::map<index, count> map = this->subsetSizeMap();
    for (auto kv : map) {
        sizes.push_back(kv.second);
    }
    return sizes;
}

std::map<index, count> Partition::subsetSizeMap() const {
    std::map<index, count> subset2size;

    this->forEntries([&](index, index s){
        if (s != none) {
            subset2size[s] += 1;
        }
    });

    return subset2size;
}

std::set<index> Partition::getMembers(const index s) const {
    assert (s <= omega);
    std::set<index> subset;
    for (index e = 0; e < this->z; ++e) {
        if (data[e] == s) {
            subset.insert(e);
        }
    }
    return subset;
}

std::vector<index> Partition::getVector() const {
    return this->data; //FIXME is this appropriate? - why not?
}

std::vector<index> Partition::moveVector() {
    return std::move(this->data);
}


std::set<std::set<index> > Partition::getSubsets() const {
    std::vector<std::set<index> > table(omega+1);
    this->forEntries([&](index e, index s){
        assert(s <= omega);
        table[s].insert(e);
    });

    std::set<std::set<index> > subsets;
    for (auto const &set : table) {
        if (set.size() > 0) {
            subsets.insert(set);
        }
    }
    return subsets;
}

void Partition::allToOnePartition() {
    omega = 0;
    this->parallelForEntries([&](index e, index) {
        this->data[e] = 0;
    }, upperBound() > (1 << 20));
}

std::set<index> Partition::getSubsetIds() const {
    std::set<index> ids;
    for (index id : data) {
        if (id != none) {
            ids.insert(id);
        }
    }
    return ids;
}

} /* namespace NetworKit */
