/*
 * PrioQueue.h
 *
 *  Created on: 21.02.2014
 *      Author: Henning
 */

#ifndef PRIOQUEUE_H_
#define PRIOQUEUE_H_

#include <cassert>
#include <set>
#include <vector>
#include <limits>

namespace Aux {


/**
 * Priority queue with extract-min and decrease-key.
 * The type Val takes on integer values between 0 and n-1.
 * O(n log n) for construction, O(log n) for typical operations.
 */
template<class Key, class Val>
class PrioQueue {
	typedef std::pair<Key, Val> ElemType;

protected:
	std::set<ElemType> pqset;
	std::vector<Key> mapValToKey;

public:
	/**
	 * Builds priority queue from the vector @a elems.
	 */
	PrioQueue(const std::vector<ElemType>& elems);

	/**
	 * Builds priority queue from the vector @a keys, values are indices
	 * in @a keys.
	 */
	PrioQueue(std::vector<Key>& keys);

	virtual ~PrioQueue() {}

	/**
	 * Inserts key-value pair stored in @a elem.
	 */
	virtual void insert(Key key, Val value);

	/**
	 * Removes the element with minimum key and returns it.
	 */
	virtual ElemType extractMin();

	/**
	 * Modifies entry with key @a elem.first and sets its value
	 * to @a elem.second. If the corresponding key is not present,
	 * the element will be inserted.
	 */
	virtual void decreaseKey(Key key, Val newValue);

	/**
	 * Removes key-value pair given by @a elem.
	 */
	virtual void remove(const ElemType& elem);

	/**
	 * Removes key-value pair given by value @a val.
	 */
	virtual void remove(const Val& val);

	/**
	 * @return Number of elements in PQ.
	 */
	virtual uint64_t size() const;

	/**
	 * Removes all elements from the PQ.
	 */
	virtual void clear();

	/**
	 * DEBUGGING
	 */
//	virtual void print() {
//		DEBUG("num entries: ", mapValToKey.size());
//		for (uint64_t i = 0; i < mapValToKey.size(); ++i) {
//			std::cout << "key: " << mapValToKey[i] << ", val: " << i << std::endl;
//		}
//	}
};

} /* namespace Aux */


template<class Key, class Val>
Aux::PrioQueue<Key, Val>::PrioQueue(const std::vector<ElemType>& elems) {
	mapValToKey.resize(elems.size());
	for (auto elem: elems) {
		insert(elem.first, elem.second);
	}
}

template<class Key, class Val>
Aux::PrioQueue<Key, Val>::PrioQueue(std::vector<Key>& keys) {
	mapValToKey.resize(keys.size());
	uint64_t index = 0;
	for (auto key: keys) {
		insert(key, index);
		++index;
	}
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::insert(Key key, Val value) {
	if (value >= mapValToKey.size()) {
		uint64_t doubledSize = 2 * mapValToKey.size();
		assert(value < doubledSize);
		mapValToKey.resize(doubledSize);
	}
	pqset.insert(std::make_pair(key, value));
	mapValToKey.at(value) = key;
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::remove(const ElemType& elem) {
	Key key = mapValToKey.at(elem.second);
//	DEBUG("key: ", key);
	if (key != none) {
		pqset.erase(std::make_pair(key, elem.second));
		mapValToKey.at(elem.second) = none;
	}
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::remove(const Val& val) {
	Key key = mapValToKey.at(val);
//	DEBUG("key: ", key);
	if (key != none) {
		pqset.erase(std::make_pair(key, val));
		mapValToKey.at(val) = none;
	}
}

template<class Key, class Val>
std::pair<Key, Val> Aux::PrioQueue<Key, Val>::extractMin() {
	assert(pqset.size() > 0);
	ElemType elem = (* pqset.begin());
	remove(elem);
	mapValToKey.at(elem.second) = none;
	return elem;
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::decreaseKey(Key key, Val newValue) {
	// find and remove element with given key
	remove(std::make_pair(key, newValue));

	// insert element with new value
	insert(key, newValue);
}

template<class Key, class Val>
inline uint64_t Aux::PrioQueue<Key, Val>::size() const {
	return pqset.size();
}

template<class Key, class Val>
inline void Aux::PrioQueue<Key, Val>::clear() {
	pqset.clear();
	mapValToKey.clear();
}


#endif /* PRIOQUEUE_H_ */
