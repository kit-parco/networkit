/*
 * IOGTest.cpp
 *
 *  Created on: 12.12.2012
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#include <gtest/gtest.h>

#include <array>
#include <cassert>
#include <chrono>
#include <iostream>
#include <fstream>
#include <limits>
#include <unordered_set>
#include <vector>

#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/METISGraphWriter.hpp>
#include <networkit/io/PartitionWriter.hpp>
#include <networkit/io/PartitionReader.hpp>
#include <networkit/io/GraphIO.hpp>
#include <networkit/io/DotGraphWriter.hpp>
#include <networkit/io/DGSReader.hpp>
#include <networkit/io/EdgeListWriter.hpp>
#include <networkit/io/EdgeListPartitionReader.hpp>
#include <networkit/io/SNAPGraphReader.hpp>
#include <networkit/io/SNAPEdgeListPartitionReader.hpp>
#include <networkit/io/SNAPGraphWriter.hpp>
#include <networkit/io/EdgeListReader.hpp>
#include <networkit/io/KONECTGraphReader.hpp>
#include <networkit/io/GMLGraphWriter.hpp>
#include <networkit/io/EdgeListCoverReader.hpp>
#include <networkit/io/CoverReader.hpp>
#include <networkit/io/CoverWriter.hpp>
#include <networkit/io/GMLGraphReader.hpp>
#include <networkit/io/GraphToolBinaryReader.hpp>
#include <networkit/io/GraphToolBinaryWriter.hpp>
#include <networkit/io/ThrillGraphBinaryWriter.hpp>
#include <networkit/io/ThrillGraphBinaryReader.hpp>
#include <networkit/io/BinaryPartitionWriter.hpp>
#include <networkit/io/BinaryPartitionReader.hpp>
#include <networkit/io/BinaryEdgeListPartitionWriter.hpp>
#include <networkit/io/BinaryEdgeListPartitionReader.hpp>
#include <networkit/io/NetworkitBinaryGraph.hpp>
#include <networkit/io/NetworkitBinaryReader.hpp>
#include <networkit/io/NetworkitBinaryWriter.hpp>
#include <networkit/generators/ErdosRenyiGenerator.hpp>

#include <networkit/community/GraphClusteringTools.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/ClusteringGenerator.hpp>
#include <networkit/structures/Partition.hpp>
#include <networkit/community/Modularity.hpp>
#include <networkit/community/PLP.hpp>
#include <networkit/dynamics/GraphDifference.hpp>

#include <tlx/unused.hpp>

namespace NetworKit {

class IOGTest: public testing::Test {};

TEST_F(IOGTest, testEdgeListWriter){
    ErdosRenyiGenerator graphGen(100, 0.1, true);
    Graph G = graphGen.generate();

    std::string path = "output/edgelist2.txt";
    EdgeListWriter writer(' ', 1, true);
    EXPECT_NO_THROW(writer.write(G, path));

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    file.close();
    EXPECT_TRUE(exists) << "A file should have been created : " << path;

    EdgeListReader reader(' ', 1, "#", true, true);
    Graph G2 = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(),G2.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(),G2.numberOfEdges());
    EXPECT_EQ(G.isDirected(),G2.isDirected());
    EXPECT_EQ(G.isWeighted(),G2.isWeighted());

    // If not continuous, firstNode should be set to 0 automatically
    reader = EdgeListReader(' ', 1, "#", false, true);
    Graph G3 = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), G3.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(), G3.numberOfEdges());
    EXPECT_EQ(G.isDirected(), G3.isDirected());
    EXPECT_EQ(G.isWeighted(), G3.isWeighted());
    std::vector<node> nodes(G.nodeRange().begin(), G.nodeRange().end());
    index i = 0;
    G3.forNodes([&](node u) { EXPECT_EQ(u, nodes[i++]); });
}

TEST_F(IOGTest, testGraphIOEdgeList) {
    ErdosRenyiGenerator graphGen(100, 0.1);
    Graph G = graphGen.generate();

    GraphIO graphio;
    std::string path = "output/edgelist.txt";
    graphio.writeEdgeList(G, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(IOGTest, testGraphIOAdjacencyList) {
    ErdosRenyiGenerator graphGen(100, 0.1);
    Graph G = graphGen.generate();
    GraphIO graphio;
    std::string path = "output/circular.adjlist";
    graphio.writeAdjacencyList(G, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(IOGTest, testGraphIOForIsolatedNodes) {
    Graph G(20);
    GraphIO graphio;
    std::string path = "output/isolated.adjlist";
    graphio.writeAdjacencyList(G, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(IOGTest, testMETISGraphReader) {
    std::string path = "input/jazz.graph";

    METISGraphReader reader;
    Graph G = reader.read(path);
    count n = 198;
    count m = 2742;
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }


    // graph polblogs (has singletons)
    path = "input/polblogs.graph";
    G = reader.read(path);
    n = 1490;
    m = 16715;
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }


    // graph PGPgiantcompo
    path = "input/PGPgiantcompo.graph";
    G = reader.read(path);
    n = 10680;
    m = 24316;
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }
}

TEST_F(IOGTest, testMETISGraphReaderWithTinyGraphs) {
    /* These graphs are from the METIS documentation and cover different settings
        of the fmt flag in the header of a METIS graph file */
    count n = 7;
    count m = 11;
    METISGraphReader reader;

    std::string path = "input/tiny_01.graph";
    Graph G = reader.read(path);
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    path = "input/tiny_02.graph";
    G = reader.read(path);
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    path = "input/tiny_03.graph";
    G = reader.read(path);
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }

    path = "input/tiny_04.graph";
    G = reader.read(path);
    EXPECT_FALSE(G.isEmpty());
    EXPECT_EQ(n, G.numberOfNodes()) << "There are " << n << " nodes in the graph";
    EXPECT_EQ(m, G.numberOfEdges()) << "There are " << m << " edges in the graph";

    for (index v = 0; v < n; ++v) {
        EXPECT_TRUE(G.hasNode(v)) << "Node " << v << " should be there";
    }
}

TEST_F(IOGTest, testMETISGraphWriter) {
    std::string path = "output/jazz1.graph";
    Graph G = Graph(3);
    G.addEdge(0,2);
    G.addEdge(1,1);
    G.addEdge(1,2);
    G.addEdge(2,2);

    METISGraphWriter writer;
    writer.write(G, false, path);
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;

}

TEST_F(IOGTest, testMETISGraphWriterWithWeights) {
    std::string path = "output/jazz2.graph";
    Graph G = Graph(5);
    G.addEdge(0,2);
    G.addEdge(0,1);
    G.addEdge(0,0);
    G.addEdge(1,1);

    METISGraphWriter writer;
    writer.write(G, true, path);
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;

}

TEST_F(IOGTest, testPartitionWriterAndReader) {
    // write clustering first
    std::string path = "output/example.clust";

    count n = 100;
    count k = 3;
    ErdosRenyiGenerator graphGen(n, 0.1);
    Graph G = graphGen.generate();

    ClusteringGenerator clusteringGen;
    Partition zeta = clusteringGen.makeRandomClustering(G, k);

    PartitionWriter writer;
    writer.write(zeta, path);

    // check if file exists
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "clustering file should have been written to: " << path;


    PartitionReader reader;
    Partition read = reader.read(path);

    EXPECT_EQ(n, read.numberOfElements()) << "read clustering should contain n nodes";
    EXPECT_TRUE(GraphClusteringTools::isProperClustering(G, read)) << "read clustering should be proper clustering of G";
    EXPECT_TRUE(GraphClusteringTools::equalClusterings(read, zeta, G)) << "read clustering should be identical to created clustering";
}

TEST_F(IOGTest, testDotGraphWriter) {
    ErdosRenyiGenerator graphGen(100, 0.1);
    Graph G = graphGen.generate();

    std::string path = "output/example.dot";

    DotGraphWriter writer;
    writer.write(G, path);

    // check if file exists
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "graph file should have been written to: " << path;
}

TEST_F(IOGTest, debugDGSReaderOnBigFile) {
    // read example graph
    DGSReader reader;
    Graph G;
    GraphEventProxy Gproxy(G);
    reader.read("/Users/forigem/KIT/NetworKit-CommunityDetection/input/AuthorsGraph.dgs", Gproxy);
}

TEST_F(IOGTest, debugDGSReader) {
    // read example graph
    DGSReader reader;
    Graph G;
    GraphEventProxy Gproxy(G);
    reader.read("input/example2.dgs", Gproxy);

    // get input parameters
    count nodeCount = G.numberOfNodes();
    DEBUG("Number of nodes " , nodeCount);
    EXPECT_EQ(3u, nodeCount);
    count edgeCount = G.numberOfEdges();
    DEBUG("Number of edges " , edgeCount);
    EXPECT_EQ(2u, edgeCount);

    G.forNodes([&](node n) {
        DEBUG("DEGREE OF NODE: " , G.degree(n) , "\n");
        tlx::unused(n);
    });

}

TEST_F(IOGTest, testEdgeListReader) {
    EdgeListReader reader('\t', 1);

    std::string path = "input/network.dat";
    DEBUG("reading file: " , path);
    Graph G = reader.read(path);
    EXPECT_EQ(10u, G.numberOfNodes());
    EXPECT_EQ(10u, G.numberOfEdges());
    EXPECT_TRUE(G.hasEdge(0, 5));
    EXPECT_TRUE(G.hasEdge(2, 9));
    EXPECT_TRUE(G.hasEdge(1, 7));

    path = "input/example.edgelist";
    DEBUG("reading file: " , path);
    EdgeListReader reader2('\t', 1);
    Graph G2 = reader2.read(path);
    EXPECT_EQ(10u, G2.numberOfEdges());
    EXPECT_TRUE(G2.hasEdge(0, 4));

    path = "input/spaceseparated.edgelist";
    DEBUG("reading file: " , path);
    EdgeListReader reader3(' ', 1);
    Graph G3 = reader3.read(path);
    EXPECT_EQ(10u, G3.numberOfEdges());
    EXPECT_TRUE(G3.hasEdge(0, 4));

    path = "input/spaceseparated_weighted.edgelist";
    DEBUG("reading file: " , path);
    Graph G32 = reader3.read(path);
    EXPECT_TRUE(G32.isWeighted());
    EXPECT_EQ(2,G32.weight(0,1));
    EXPECT_EQ(4,G32.weight(0,2));
    EXPECT_EQ(3,G32.weight(1,2));

    path = "input/comments.edgelist";
    DEBUG("reading file: " , path);
    EdgeListReader reader4('\t', 1);
    Graph G4 = reader4.read(path);
    EXPECT_EQ(10u, G4.numberOfEdges());
    EXPECT_TRUE(G4.hasEdge(0, 4));

    path = "input/alphabet.edgelist";
    DEBUG("reading file: " , path);
    EdgeListReader reader5('\t', 0, "#", false, false);
    Graph G5 = reader5.read(path);
    EXPECT_EQ(5u, G5.numberOfNodes());
    EXPECT_EQ(4u, G5.numberOfEdges());
    EXPECT_TRUE(G5.hasEdge(0, 1));
    EXPECT_TRUE(G5.hasEdge(0, 2));
    EXPECT_EQ(5, G5.weight(3,4));
    EXPECT_EQ(1, G5.weight(2,3));
}

TEST_F(IOGTest, testEdgeListPartitionReader) {
    EdgeListPartitionReader reader(1);

    Partition zeta = reader.read("input/community.dat");
    //EXPECT_EQ(10, zeta.size());
    EXPECT_EQ(1u, zeta[0]);
    EXPECT_EQ(3u, zeta[1]);
    EXPECT_EQ(2u, zeta[2]);
    EXPECT_EQ(10u, zeta.numberOfElements());

}

TEST_F(IOGTest, testEdgeListCoverReader) {
    EdgeListCoverReader reader(1);
    EdgeListReader gReader('\t', 1);

    Graph G = gReader.read("input/network_overlapping.dat");
    Cover zeta = reader.read("input/community_overlapping.dat", G);
    EXPECT_EQ(9u, zeta.upperBound());
    EXPECT_EQ(10u, zeta.numberOfElements());
    EXPECT_EQ(1u, zeta[0].count(1));
    EXPECT_EQ(3u, zeta[0].size());
    EXPECT_EQ(1u, zeta[3].size());
}

TEST_F(IOGTest, testCoverReader) {
    CoverReader reader;
    EdgeListReader gReader('\t', 1);

    Graph G = gReader.read("input/network_overlapping.dat");
    Cover zeta = reader.read("input/community_overlapping.cover", G);
    EXPECT_EQ(9u, zeta.upperBound());
    EXPECT_EQ(10u, zeta.numberOfElements());
    EXPECT_EQ(1u, zeta[0].count(1));
    EXPECT_EQ(3u, zeta[0].size());
    EXPECT_EQ(1u, zeta[3].size());
}

TEST_F(IOGTest, testCoverWriter) {
    std::string outpath = "output/coverWriter_test.cover";
    CoverWriter writer;
    CoverReader reader;
    EdgeListReader gReader('\t', 1);

    Graph G = gReader.read("input/network_overlapping.dat");
    Cover zeta = reader.read("input/community_overlapping.cover", G);

    writer.write(zeta, outpath);

    Cover read = reader.read(outpath, G);
    EXPECT_EQ(9u, read.upperBound());
    EXPECT_EQ(10u, read.numberOfElements());
    EXPECT_EQ(1u, read[0].count(1));
    EXPECT_EQ(3u, read[0].size());
    EXPECT_EQ(1u, read[3].size());
}

TEST_F(IOGTest, testMETISGraphReaderForNodeExistence2) {
    METISGraphReader reader;
    Graph G = reader.read("input/jazz.graph");
    EXPECT_TRUE(G.hasNode(0));
    EXPECT_EQ(198u, G.numberOfNodes());
    EXPECT_EQ(2742u, G.numberOfEdges());
}

TEST_F(IOGTest, testMETISGraphReaderWithIsolatedNode) {
    METISGraphReader reader;
    Graph G = reader.read("input/example.graph");
    EXPECT_EQ(4u, G.numberOfNodes());
    EXPECT_EQ(2u, G.numberOfEdges());
    EXPECT_TRUE(G.hasNode(0));
    EXPECT_TRUE(G.hasNode(1));
    EXPECT_TRUE(G.hasNode(2));
    EXPECT_TRUE(G.hasNode(3));
    EXPECT_TRUE(G.hasEdge(0,1));
    EXPECT_TRUE(G.hasEdge(0,3));
}

TEST_F(IOGTest, debugReadingLFR) {
    std::string graphPath;
    std::string clustPath;

    std::cout << "[INPUT] LFR graph file path >" << std::endl;
    std::getline(std::cin, graphPath);

    std::cout << "[INPUT] clustering file path >" << std::endl;
    std::getline(std::cin, clustPath);

    EdgeListReader graphReader('\t',1);
    EdgeListPartitionReader clusteringReader;

    Graph G = graphReader.read(graphPath);
    Partition truth = clusteringReader.read(clustPath);

    PLP PLP(G);
    PLP.run();
    Partition zeta = PLP.getPartition();

    Modularity mod;
    INFO("static clustering quality: " , mod.getQuality(zeta, G));
    INFO("static clustering number of clusters: " , zeta.numberOfSubsets());
    INFO("ground truth quality: " , mod.getQuality(truth, G));
    INFO("ground truth number of clusters: " , truth.numberOfSubsets());

}

TEST_F(IOGTest, debugReadingSNAP) {
    std::string graphPath;

    std::cout << "[INPUT] SNAP graph file path >" << std::endl;
    std::getline(std::cin, graphPath);

    EdgeListReader graphReader(' ', 1);

    Graph G = graphReader.read(graphPath);

    INFO("n = " , G.numberOfNodes());
    INFO("m = " , G.numberOfEdges());

}

TEST_F(IOGTest, testSNAPGraphReader) {
    SNAPGraphReader reader(true);
    Graph G = reader.read("input/wiki-Vote.txt");

    EXPECT_TRUE(G.isDirected());
    EXPECT_EQ(G.numberOfNodes(), 7115);
    EXPECT_EQ(G.numberOfEdges(), 103689);
}

TEST_F(IOGTest, testSNAPGraphWriter) {
    METISGraphReader reader;
    Graph G = reader.read("input/jazz.graph");

    std::string path = """output/SNAPGraphWriter.gr";
    SNAPGraphWriter writer;
    writer.write(G, path);

    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "graph file should have been written to: " << path;
}

TEST_F(IOGTest, debugReadingMETISFile) {
    std::string graphPath;
    std::cout << "[INPUT] graph file path >" << std::endl;
    std::getline(std::cin, graphPath);

    METISGraphReader reader;
    Graph G = reader.read(graphPath);

    EXPECT_TRUE(true);
}

TEST_F(IOGTest, testGMLGraphWriterUndirected) {
    std::string path = "output/jazz2_undirected.gml";
    Graph G = Graph(5);
    G.addEdge(0,2);
    G.addEdge(0,1);
    G.addEdge(0,0);
    G.addEdge(1,1);

    GMLGraphWriter writer;
    writer.write(G,path);
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;


}

TEST_F(IOGTest, testGMLGraphWriterDirected) {
    std::string path = "output/jazz2_directed.gml";
    Graph G = Graph(5,false,true);
    G.addEdge(0,2);
    G.addEdge(0,1);
    G.addEdge(0,0);
    G.addEdge(1,1);

    GMLGraphWriter writer;
    writer.write(G,path);
    bool exists = false;
    std::ifstream file(path);
    if (file) {
        exists = true;
    }
    EXPECT_TRUE(exists) << "A file should have been created : " << path;
}

TEST_F(IOGTest, testGMLGraphReaderUndirected) {
    std::string path = "input/jazz2_undirected.gml";
    GMLGraphReader reader;
    Graph G = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), 5u) << "number of nodes is not correct";
    EXPECT_TRUE(G.hasEdge(0,2));
    EXPECT_TRUE(G.hasEdge(0,1));
    EXPECT_TRUE(G.hasEdge(0,0));
    EXPECT_TRUE(G.hasEdge(1,1));
    EXPECT_FALSE(G.isDirected());
    EXPECT_TRUE(G.hasEdge(2,0));
    EXPECT_TRUE(G.hasEdge(1,0));
}

TEST_F(IOGTest, testGMLGraphReaderDirected) {
    std::string path = "input/jazz2_directed.gml";
    GMLGraphReader reader;
    Graph G = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(), 5u) << "number of nodes is not correct";
    EXPECT_TRUE(G.hasEdge(0,2));
    EXPECT_TRUE(G.hasEdge(0,1));
    EXPECT_TRUE(G.hasEdge(0,0));
    EXPECT_TRUE(G.hasEdge(1,1));
    EXPECT_TRUE(G.isDirected());
    EXPECT_FALSE(G.hasEdge(2,0));
    EXPECT_FALSE(G.hasEdge(1,0));

}

TEST_F(IOGTest, testGraphToolBinaryReader) {
    std::string path = "input/power.gt";
    GraphToolBinaryReader reader;
    Graph G = reader.read(path);
    EXPECT_EQ(4941u,G.numberOfNodes());
    EXPECT_EQ(6594u,G.numberOfEdges());
    EXPECT_FALSE(G.isDirected());
}

TEST_F(IOGTest, testGraphToolBinaryWriter) {
    Graph G(10,false,false);
    G.addEdge(0,1);
    G.addEdge(2,1);
    G.addEdge(2,3);
    G.addEdge(3,4);
    G.addEdge(5,4);
    G.addEdge(5,6);
    G.addEdge(7,6);
    G.addEdge(8,6);
    G.addEdge(7,8);
    G.addEdge(9,8);
    G.addEdge(9,0);
    GraphToolBinaryReader reader;
    GraphToolBinaryWriter writer;
    std::string path = "output/test.gt";
    writer.write(G,path);
    Graph Gread = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(),Gread.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(),Gread.numberOfEdges());
    EXPECT_EQ(G.isDirected(),Gread.isDirected());
    EXPECT_EQ(G.isWeighted(),Gread.isWeighted());
}

TEST_F(IOGTest, testGraphToolBinaryWriterWithDeletedNodes) {
    Graph G(10,false,false);
    G.removeNode(0);
    G.addEdge(2,1);
    G.addEdge(2,3);
    G.removeNode(4);
    G.addEdge(5,6);
    G.addEdge(7,6);
    G.addEdge(8,6);
    G.addEdge(7,8);
    G.removeNode(9);
    GraphToolBinaryReader reader;
    GraphToolBinaryWriter writer;
    std::string path = "output/test.gt";
    writer.write(G,path);
    Graph Gread = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(),Gread.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(),Gread.numberOfEdges());
    EXPECT_EQ(G.isDirected(),Gread.isDirected());
    EXPECT_EQ(G.isWeighted(),Gread.isWeighted());
}

TEST_F(IOGTest, testGraphToolBinaryWriterDirected) {
    Graph G(10,false,true);
    G.addEdge(0,1);
    G.addEdge(2,1);
    G.addEdge(2,3);
    G.addEdge(3,4);
    G.addEdge(5,4);
    G.addEdge(5,6);
    G.addEdge(7,6);
    G.addEdge(8,6);
    G.addEdge(7,8);
    G.addEdge(9,8);
    G.addEdge(9,0);
    GraphToolBinaryReader reader;
    GraphToolBinaryWriter writer;
    std::string path = "output/test.gt";
    writer.write(G,path);
    Graph Gread = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(),Gread.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(),Gread.numberOfEdges());
    EXPECT_EQ(G.isDirected(),Gread.isDirected());
    EXPECT_EQ(G.isWeighted(),Gread.isWeighted());
}

TEST_F(IOGTest, testGraphToolBinaryWriterWithDeletedNodesDirected) {
    Graph G(10,false,true);
    G.removeNode(0);
    G.addEdge(2,1);
    G.addEdge(2,3);
    G.removeNode(4);
    G.addEdge(5,6);
    G.addEdge(7,6);
    G.addEdge(8,6);
    G.addEdge(7,8);
    G.removeNode(9);
    GraphToolBinaryReader reader;
    GraphToolBinaryWriter writer;
    std::string path = "output/test.gt";
    writer.write(G,path);
    Graph Gread = reader.read(path);
    EXPECT_EQ(G.numberOfNodes(),Gread.numberOfNodes());
    EXPECT_EQ(G.numberOfEdges(),Gread.numberOfEdges());
    EXPECT_EQ(G.isDirected(),Gread.isDirected());
    EXPECT_EQ(G.isWeighted(),Gread.isWeighted());
}


TEST_F(IOGTest, testThrillGraphBinaryWriterAndReader) {
    // This test graph has a large maximum degree as degrees smaller than 128
    // do not test the binary writer and reader properly.
    // So we simply use a star.
    count n = 257;
    Graph G(n, false, false);
    node center = 129;

    for (node u = 0; u < n; ++u) {
        if (u != center) {
            G.addEdge(u, center);
        }
    }

    std::string path = "output/test.thrillbin";

    ThrillGraphBinaryReader reader;
    ThrillGraphBinaryWriter writer;

    writer.write(G, path);
    Graph H = reader.read(path);

    GraphDifference diff(G, H);
    diff.run();

    EXPECT_EQ(diff.getEdits().size(), 0);
}

TEST_F(IOGTest, testBinaryPartitionWriterAndReader) {
    Partition P(5);
    P.setUpperBound((1ull<<32));
    P[0] = 0;
    P[1] = 2007;
    P[2] = none;
    P[3] = (1ull << 31 | 1ull << 28);
    P[4] = 3925932491;

    std::string path = "output/partition.bin";

    BinaryPartitionWriter writer(8);
    BinaryPartitionReader reader(8);

    writer.write(P, path);
    Partition Q(reader.read(path));

    EXPECT_EQ(P[0], Q[0]);
    EXPECT_EQ(P[1], Q[1]);
    EXPECT_EQ(P[2], Q[2]);
    EXPECT_EQ(P[3], Q[3]);
    EXPECT_EQ(P[4], Q[4]);
    EXPECT_EQ(Q.upperBound(), P[4]+1);
}

TEST_F(IOGTest, testBinaryEdgeListPartitionWriterAndReader) {
    Partition P(5);
    P.setUpperBound((1ull<<32));
    P[0] = 0;
    P[1] = 2007;
    P[2] = none;
    P[3] = (1ull << 31 | 1ull << 28);
    P[4] = 3925932491;

    std::string path = "output/partition.bin";

    BinaryEdgeListPartitionWriter writer(1, 8);
    BinaryEdgeListPartitionReader reader(1, 8);

    writer.write(P, path);
    Partition Q(reader.read(path));

    EXPECT_EQ(P[0], Q[0]);
    EXPECT_EQ(P[1], Q[1]);
    EXPECT_EQ(P[2], Q[2]);
    EXPECT_EQ(P[3], Q[3]);
    EXPECT_EQ(P[4], Q[4]);
    EXPECT_EQ(Q.upperBound(), P[4]+1);
}

TEST_F(IOGTest, testKONECTGraphReader){
    KONECTGraphReader reader;
    Graph G = reader.read("input/foodweb-baydry.konect");

    ASSERT_TRUE(G.isDirected());
    ASSERT_EQ(G.numberOfEdges() , 2137);
    ASSERT_EQ(G.numberOfNodes() , 128);
    ASSERT_EQ(G.weight(0,1), 1.261404);
    ASSERT_EQ(G.weight(127, 48), 0.03050447);
}
TEST_F(IOGTest, testNetworkitBinaryTiny01) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/tiny_01.graph");
    NetworkitBinaryWriter writer;

    writer.write(G, "output/binary_tiny01");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_tiny01");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryTiny01Indexed) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/tiny_01.graph");
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::autoDetect, true);

    G.indexEdges();
    writer.write(G, "output/binary_tiny01");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_tiny01");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());

    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
            ASSERT_EQ(G2.edgeId(u,v), G.edgeId(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryKonect) {
    KONECTGraphReader reader2;
    Graph G = reader2.read("input/foodweb-baydry.konect");
    NetworkitBinaryWriter writer;

    writer.write(G, "output/binary_konect");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_konect");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), true);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
            ASSERT_EQ(G.weight(u,v), G2.weight(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryKonectIndexed) {
    KONECTGraphReader reader2;
    Graph G = reader2.read("input/foodweb-baydry.konect");
    G.indexEdges();
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::autoDetect, true);
    writer.write(G, "output/binary_konect");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_konect");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), true);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
            ASSERT_EQ(G.weight(u,v), G2.weight(u,v));
            ASSERT_EQ(G.edgeId(u,v), G2.edgeId(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryJazz) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/jazz.graph");

    NetworkitBinaryWriter writer;
    writer.write(G, "output/binary_jazz");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_jazz");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryJazzIndexed) {
    METISGraphReader reader2;
    Graph G = reader2.read("input/jazz.graph");
    G.indexEdges();

    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::autoDetect, true);
    writer.write(G, "output/binary_jazz");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_jazz");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryWiki) {
    SNAPGraphReader reader2(true);
    Graph G = reader2.read("input/wiki-Vote.txt");
    NetworkitBinaryWriter writer;

    writer.write(G, "output/binary_wiki");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_wiki");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryWikiIndexed) {
    SNAPGraphReader reader2(true);
    Graph G = reader2.read("input/wiki-Vote.txt");
    G.indexEdges();
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::autoDetect, true);

    writer.write(G, "output/binary_wiki");
    ASSERT_TRUE(!G.isEmpty());

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binary_wiki");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G2.numberOfEdges(), G.numberOfEdges());
    ASSERT_EQ(G2.numberOfNodes(), G.numberOfNodes());
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinarySignedWeights) {

    Graph G(10, true, false);
    int64_t weight = -1;
    for(count n = 0; n < G.numberOfNodes(); n++) {
        if(n != G.numberOfNodes()-1)
            G.addEdge(n, n+1, weight++);
    }
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::autoDetect);
    writer.write(G, "output/binarySigned");

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binarySigned");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), true);
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
            ASSERT_EQ(G.weight(u,v), G2.weight(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinarySignedWeightsIndexed) {

    Graph G(10, true, false);
    G.indexEdges();
    int64_t weight = -1;
    for(count n = 0; n < G.numberOfNodes(); n++) {
        if(n != G.numberOfNodes()-1)
            G.addEdge(n, n+1, weight++);
    }
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::autoDetect, true);
    writer.write(G, "output/binarySigned");

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binarySigned");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), true);
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
            ASSERT_EQ(G.weight(u,v), G2.weight(u,v));
            ASSERT_EQ(G.edgeId(u,v), G2.edgeId(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryFloatWeights) {

    Graph G(10, true, false);
    float weight = 987.654f;
    for(count n = 0; n < G.numberOfNodes(); n++) {
        if(n != G.numberOfNodes()-1)
            G.addEdge(n, n+1, weight++);
    }
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::autoDetect);
    writer.write(G, "output/binaryFloats");

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binaryFloats");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), true);
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
            ASSERT_EQ(G.weight(u,v), G2.weight(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryFloatWeightsIndexed) {

    Graph G(10, true, false);
    G.indexEdges();
    float weight = 987.654f;
    for(count n = 0; n < G.numberOfNodes(); n++) {
        if(n != G.numberOfNodes()-1)
            G.addEdge(n, n+1, weight++);
    }
    NetworkitBinaryWriter writer(32, NetworkitBinaryWeights::autoDetect, true);
    writer.write(G, "output/binaryFloats");

    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/binaryFloats");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), true);
    G.forNodes([&](node u){
        G.forEdgesOf(u, [&](node v) {
            ASSERT_TRUE(G2.hasEdge(u,v));
            ASSERT_EQ(G.weight(u,v), G2.weight(u,v));
        });
    });
}

TEST_F(IOGTest, testNetworkitBinaryUndirectedSelfLoops) {

    Graph G(5, false, false);
    G.addEdge(0,0);
    G.addEdge(1,1);
    G.addEdge(2,2);
    G.addEdge(3,3);
    G.addEdge(4,4);
    NetworkitBinaryWriter writer;
    writer.write(G, "output/loops");
    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/loops");
    EXPECT_EQ(G2.isDirected(), false);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G.numberOfSelfLoops(), G2.numberOfSelfLoops());
}

TEST_F(IOGTest, testNetworkitBinaryDirectedSelfLoops) {

    Graph G(5, false, true);
    G.addEdge(0,0);
    G.addEdge(1,1);
    G.addEdge(2,2);
    G.addEdge(3,3);
    G.addEdge(4,4);
    NetworkitBinaryWriter writer;
    writer.write(G, "output/loops");
    NetworkitBinaryReader reader;
    Graph G2 = reader.read("output/loops");
    EXPECT_EQ(G2.isDirected(), true);
    EXPECT_EQ(G2.isWeighted(), false);
    ASSERT_EQ(G.numberOfSelfLoops(), G2.numberOfSelfLoops());
}

TEST_F(IOGTest, testNetworkitBinaryVarInt) {
    std::array<uint8_t, 10> buffer;

    // write defined values into buffer
    {
        uint8_t i = 0;
        for(auto& x : buffer)
            x = i++;
    }

    std::mt19937_64 gen{1};

    uint64_t checked_bits = 0;

    for (int bits = 0; bits < 64; ++bits) {
        auto min = uint64_t(1) << bits;
        auto max = 2*min - 1;

        // special cases
        if (bits == 0) {
            min = 0;
            max = 0;
        } else if (bits == 64) {
            max = std::numeric_limits<uint64_t>::max();
        }

        std::uniform_int_distribution<uint64_t> distr{min, max};

        const auto nSamples = std::min<size_t>(10 * max + 2, 1000);

        for (size_t i = 0; i < nSamples; ++i) {
            // first two iterations test min/max values, all other random values
            const auto orig = [&] {
                if (i == 0) return min;
                if (i == 1) return max;
                return distr(gen);
            }();

            uint64_t valueRead;
            const auto bytesWritten = nkbg::varIntEncode(orig, buffer.data());
            const auto bytesRead    = nkbg::varIntDecode(buffer.data(), valueRead);

            ASSERT_EQ(bytesWritten, bytesRead) << "bits=" << bits << ", i=" << i;
            ASSERT_EQ(valueRead, orig)         << "bits=" << bits << ", i=" << i;
            ASSERT_GT(buffer[bytesWritten-1], 0) << "bits=" << bits << ", i=" << i;

            for(size_t j = bytesWritten; j < buffer.size(); ++j)
                ASSERT_EQ(buffer[j], j);

            checked_bits |= orig;
        }
    }

    // make sure we touched each bit at least once
    ASSERT_EQ(checked_bits, std::numeric_limits<uint64_t>::max());

}

TEST_F(IOGTest, testNetworkitBinaryZigzag) {
    std::mt19937_64 gen(1);
    std::uniform_int_distribution<uint64_t> distr(0, (std::numeric_limits<uint64_t>::max() >> 1) - 1);

    for(int i = 0; i < 10000; ++i) {
        auto check = [] (int64_t value) {
            const auto encoded = nkbg::zigzagEncode(value);
            const auto decoded = nkbg::zigzagDecode(encoded);

            ASSERT_EQ(value, decoded);
            ASSERT_LE(encoded, 2u * static_cast<uint64_t>(std::abs(value)));
        };

        const auto x = distr(gen);
        check(x);
        check(-1 * x);
    }
}

} /* namespace NetworKit */
