/*
 * Cover.h
 *
 *  Created on: 03.10.2013
 *      Author: cls
 */

#ifndef NETWORKIT_STRUCTURES_COVER_HPP_
#define NETWORKIT_STRUCTURES_COVER_HPP_

#include <cinttypes>
#include <set>
#include <vector>
#include <map>
#include <cassert>
#include <limits>

#include <networkit/Globals.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup structures
 * Implements a cover of a set, i.e. an assignment of
 * its elements to possibly overlapping subsets.
 */
class Cover {

public:
    /** Default constructor */
    Cover();

    /**
     * Create a new cover data structure for elements up to a maximum element index.
     *
     * @param[in]	z	maximum index
     */
    Cover(index z);

    /**
     * Creates a new cover data structure which contains the given partition.
     *
     * @param[in]	p	The partition to construct the cover from
     */
    Cover(const Partition &p);

    /** Default destructor */
    virtual ~Cover() = default;


    /**
     *  Index operator.
     *
     *  @param[in]	e	an element
     */
    inline std::set<index>& operator [](const index& e) {
        return this->data[e];
    }
    /**
     * Index operator for const instances of this class.
     *
     * @param[in]	e	an element
     */
    inline const std::set<index>& operator [](const index& e) const {
        return this->data[e];
    }

    /**
     * Return the ids of subsets in which the element @a e is contained.
     *
     * @param[in]	e	an element
     * @return A set of subset ids in which @a e is contained.
     */
    inline std::set<index> subsetsOf(index e) const {
        // TODO: assert (e < this->numberOfElements());
        return this->data[e];
    }



    /**
     * Check if cover assigns a valid subset to the element @a e.
     *
     * @param e an element.
     * @return @c true, if @a e is assigned to a valid subset, @c false otherwise.
     */
    bool contains(index e) const;


    /**
     * Check if two elements @a e1 and @a e2 belong to the same subset.
     *
     * @param e1 an element.
     * @param e2 an element.
     * @return @c true, if @a e1 and @a e2 belong to the same subset, @c false otherwise.
     */
    bool inSameSubset(index e1, index e2) const;


    /**
     * Get the members of a specific subset @a s.
     *
     * @return The set of members of subset @a s.
     */
    std::set<index> getMembers(const index s) const;


    /**
     * Add the (previously unassigned) element @a e to the set @a s.
     * @param[in]	s	a subset
     * @param[in]	e	an element
     */
    void addToSubset(index s, index e);


    /**
     * Remove the element @a e from the set @a s.
     * @param[in]	s	a subset
     * @param[in]	e	an element
     */
    void removeFromSubset(index s, index e);


    /**
     * Move the element @a e to subset @a s, i.e. remove it from all
     * other subsets and place it in the subset.
     * 	@param[in]	s	a subset
     *  @param[in]	e	an element
     */
    void moveToSubset(index s, index e);


    /**
     * Creates a singleton set containing the element @a e and returns the index of the new set.
     * @param[in]	e	an element
     * @return The index of the new set.
     */
    index toSingleton(index e);


    /**
     * Assigns every element to a singleton set.
     * Set id is equal to element id.
     */
    void allToSingletons();


    /**
     * Assigns the elements from both sets to a new set.
     * @param[in]	s	a subset
     * @param[in]	t	a subset
     */
    void mergeSubsets(index s, index t);


    /**
     * Get an upper bound for the subset ids that have been assigned.
     * (This is the maximum id + 1.)
     *
     * @return An upper bound.
     */
    index upperBound() const;

    /**
     * Get a lower bound for the subset ids that have been assigned.
     * @return A lower bound.
     */
    index lowerBound() const;


    /**
     * Get a list of subset sizes. Indices do not necessarily correspond to subset ids.
     *
     * @return A list of subset sizes.
     */
    std::vector<count> subsetSizes() const;


    /**
     * Get a map from subset id to size of the subset.
     *
     * @return A map from subset id to size of the subset.
     */
    std::map<index, count> subsetSizeMap() const;


    /**
     * Get the current number of sets in this cover.
     *
     * @return The number of sets in this cover.
     */
    count numberOfSubsets() const;

    /**
     * Get the current number of elements in this cover.
     *
     * @return The current number of elements.
     */
    count numberOfElements() const;

    /**
     * Add an additional element (node).
     */
    index extend();

    /**
     * Get the ids of nonempty subsets.
     *
     * @return A set of ids of nonempty subsets.
     */
    std::set<index> getSubsetIds() const;

    /**
     * Sets an upper bound for the subset ids that CAN be assigned.
     *
     * @param[in]	upper	highest assigned subset ID + 1
     */
    void setUpperBound(index upper);

    /**
     * Returns all subsets. Indicies in the returned vector do not correspond to subset ids.
     * @return The subsets of the cover.
     */
    std::vector<std::set<index>> getSubsets() const;

    /**
     * Add a subset of nodes.
     * @param[in]	subset	a set of nodes
     */
    void addSubset(const std::set<index> &subset);


    /**
     * Iterate over all entries (node, subset ID of node) and execute callback function @a func (lambda closure).
     *
     * @param func Takes parameters <code>(node, index)</code>
     */
    template<typename Callback> void forEntries(Callback func) const;


    /**
     * Iterate over all entries (node, subset ID of node) in parallel and execute callback function @a func (lambda closure).
     *
     * @param func Takes parameters <code>(node, index)</code>
     */
    template<typename Callback> void parallelForEntries(Callback handle) const;


private:

    index z;	//!< maximum element index that can be mapped
    index omega;	//!< maximum subset index ever assigned
    std::vector<std::set<index>> data;	//!< data container, indexed by element id, containing set of subset ids


    /**
     * Allocates and returns a new subset id.
     */
    inline index newSubsetId() {
        omega++;
        index s = omega;
        return s;
    }

};

template<typename Callback>
inline void Cover::forEntries(Callback handle) const {
    for (index e = 0; e <= this->z; e += 1) {
        handle(e, data[e]);
    }
}

template<typename Callback>
inline void Cover::parallelForEntries(Callback handle) const {
    #pragma omp parallel for
    for (omp_index e = 0; e <= static_cast<omp_index>(this->z); e += 1) {
        handle(e, data[e]);
    }
}

} /* namespace NetworKit */

#endif // NETWORKIT_STRUCTURES_COVER_HPP_
