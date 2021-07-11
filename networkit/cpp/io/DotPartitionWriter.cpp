// no-networkit-format
#include <fstream>
#include <set>

#include <networkit/io/DotPartitionWriter.hpp>

namespace NetworKit {

std::map<index, double> DotPartitionWriter::createHueMap(Graph &graph, Partition& zeta) const {
    std::set<index> uniqueIds;

    graph.forNodes([&](node u){
        if (graph.degree(u) == 0) {
            return;
        }
        index c = zeta.subsetOf(u);
        uniqueIds.insert(c);
    });

    std::map<index, double> clusterHueMap;

    int idx = 0;
    double factor = 1.0 / uniqueIds.size();
    for(const auto& value: uniqueIds) {
        clusterHueMap.insert(std::make_pair(value, factor * idx));
        ++idx;
    }

    return clusterHueMap;
}

void DotPartitionWriter::write(Graph &graph, Partition &zeta, const std::string &path) const {
    std::ofstream file{path};

    auto hueMap = this->createHueMap(graph, zeta);

    file << "graph {\n";

    graph.forNodes([&](node u){
        if (graph.degree(u) == 0) {
            return;
        }
        index c = zeta.subsetOf(u);
        double h = hueMap.at(c);
        file << u << " [style=filled, color=\"" << h << ",0.99,0.99\", label=" << c << "];\n";
    });
    graph.forEdges([&](node u, node v){
        file << u << " -- " << v << ";\n";
    });

    file << "}\n";
}

} /* namespace NetworKit */
