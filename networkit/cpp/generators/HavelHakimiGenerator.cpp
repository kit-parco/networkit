// no-networkit-format
/*
 * HavelHakimiGenerator.cpp
 *
 *  Created on: Dec 10, 2013
 *      Author: Henning
 *      Contributors: Hoske/Weisbarth
 */

#include <list>
#include <stack>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/generators/HavelHakimiGenerator.hpp>

namespace NetworKit {

HavelHakimiGenerator::HavelHakimiGenerator(const std::vector<count> &sequence, bool ignoreIfRealizable) :
        StaticDegreeSequenceGenerator(sequence), ignoreIfRealizable(ignoreIfRealizable) {}


Graph HavelHakimiGenerator::generate() {
    count n = seq.size();

    Graph G(n);
    count numDegVals = (* std::max_element(seq.begin(), seq.end())) + 1;

    // Havel-Hakimi algorithm, adapted with appropriate data structure for linear time
    // bucket data structure: vector of lists
    // a list contains all vertices with the deficit corresponding to that list
    // in each iteration the first element of the highest non-empty list is chosen and removed
    // this vertex is connected to the following elements
    // when done, the other connected vertices are moved to the next list in reverse order

    using DeficitAndNode = std::pair<count, node>;
    using Buckets = std::vector<std::list<DeficitAndNode>>;

    // put nodes in appropriate lists
    Buckets nodesByDeficit(numDegVals);
    for(node v = 0; v < n; v++) {
        nodesByDeficit[seq[v]].push_front(std::make_pair(seq[v], v));
    }

    index maxDeficit = numDegVals - 1;
    while (maxDeficit) {
        // process node in largest bucket
        while(! nodesByDeficit[maxDeficit].empty()) {
            // get element
            std::list<DeficitAndNode>::iterator listIter = nodesByDeficit[maxDeficit].begin();
            count deficit = listIter->first;
            node currentVertex = listIter->second;

            // delete it
            nodesByDeficit[maxDeficit].pop_front();

            // connect corresponding vertex with the following ones
            index currentNeighborList = maxDeficit;
            std::stack<count> numToMove;

            while (deficit) {
                count numDeleteFromCurrentList = 0;

                // search for candidates in current list
                for (auto elem : nodesByDeficit[currentNeighborList]) {
                    // connect
                    node nextNeighbor = elem.second;
                    G.addEdge(currentVertex, nextNeighbor);

                    --deficit;
                    ++numDeleteFromCurrentList;

                    if (deficit == 0) {
                        ++currentNeighborList; // due to -- a few lines below
                        break;
                    }
                }
                numToMove.push(numDeleteFromCurrentList);

                if (currentNeighborList == 1) {
                    // We are at the end of the list of nodes that have degree 1 left, i.e. the current node has already all possible neighbors
                    // This means that the Havel-Hakimi algorithm cannot generate a graph with the given degree sequence which means that
                    // the degree sequence is not realizable.
                    WARN("Degree sequence not realizable, node ", currentVertex, " should have got ", seq[currentVertex], " neighbors but actually only got ", seq[currentVertex] - deficit, " neighbors.");
                    if (!ignoreIfRealizable) throw std::runtime_error("Degree sequence is not realizable");
                    break; // if we break here, we will continue as if currentVertex had already got all the neighbors it needed.
                }

                --currentNeighborList;
            }

            while (! numToMove.empty()) {
                // get head element
                count num = numToMove.top();
                numToMove.pop();

                // move this many items from current list to the next one
                for (index i = 0; i < num; ++i) {
                    DeficitAndNode dan = nodesByDeficit[currentNeighborList].front();
                    dan.first--;
                    nodesByDeficit[currentNeighborList - 1].push_front(dan);
                    nodesByDeficit[currentNeighborList].pop_front();
                }

                ++currentNeighborList;
            }
        }
        maxDeficit--;
    }

    G.shrinkToFit();
    return G;
}

} /* namespace NetworKit */
