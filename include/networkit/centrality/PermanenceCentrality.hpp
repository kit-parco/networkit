// no-networkit-format
#ifndef NETWORKIT_CENTRALITY_PERMANENCE_CENTRALITY_HPP_
#define NETWORKIT_CENTRALITY_PERMANENCE_CENTRALITY_HPP_

#include <networkit/centrality/Centrality.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

// TODO documentation completely missing!
class PermanenceCentrality : public Algorithm {
public:
    PermanenceCentrality(const Graph &G, const Partition &P);
    void run() override;
    double getPermanence(node u);
    double getIntraClustering(node u);
private:
    const Graph &G;
    const Partition &P;
    std::vector<index> inBegin;
    std::vector<node> inEdges;
    std::vector<bool> marker;
};


}

#endif // NETWORKIT_CENTRALITY_PERMANENCE_CENTRALITY_HPP_
