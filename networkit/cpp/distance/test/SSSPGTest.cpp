// no-networkit-format
/*
 * SSSPGTest.cpp
 *
 *  Created on: 21.07.2014
 *      Author: ebergamini
 */

#include <gtest/gtest.h>

#include <networkit/distance/DynBFS.hpp>
#include <networkit/distance/BFS.hpp>
#include <networkit/distance/DynDijkstra.hpp>
#include <networkit/distance/Dijkstra.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/auxiliary/Log.hpp>

#include <stack>

namespace NetworKit {

class SSSPGTest: public testing::Test{};

TEST_F(SSSPGTest, testDijkstra) {
    /* Graph:
             ______
        /      \
         0    3   6
        \  / \ /
         2    5
        /  \ / \
         1    4   7
    */
    int n = 8;
    Graph G(n, true);

    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);
    G.addEdge(0, 6);

    Dijkstra sssp(G, 5, true, true);
    EXPECT_NO_THROW(sssp.run());
    std::vector<node> stack = sssp.getNodesSortedByDistance();
    for (count i = 0; i < stack.size() - 1; ++i)
        EXPECT_LE(sssp.distance(stack[i]), sssp.distance(stack[i + 1]));
}

TEST_F(SSSPGTest, testShortestPaths) {
    METISGraphReader reader;
    Graph G = reader.read("input/PGPgiantcompo.graph");
    DEBUG("The graph has been read.");
    int source = 2;
    BFS bfs(G, source);
    bfs.run();
    bigfloat max = 0;
    node x;
    G.forNodes([&](node n){
        if(bfs.numberOfPaths(n) > max) {
            max = bfs.numberOfPaths(n);
            x = n;
        }
    });
    count dist = bfs.distance(x);
    std::set<std::vector<node>> paths = bfs.getPaths(x, true);
    count i = 0;
    for (auto path : paths) {
        DEBUG("Path number ", i);
        i ++;
        DEBUG(path);
        EXPECT_EQ(path[0], source);
        EXPECT_EQ(path[dist], x);
    }
    DEBUG("Maximum number of shortest paths: ", bfs.numberOfPaths(x));
    DEBUG("Distance: ", dist);
}

TEST_F(SSSPGTest, testGetAllShortestPaths) {
/* Graph:

       0    3   6   9
        \  / \ / \ /
         2    5   8
        /  \ / \ / \
       1    4   7   10
*/
    int n = 11;
    Graph G(n, false);
    G.addEdge(0, 2);
    G.addEdge(1, 2);
    G.addEdge(2, 3);
    G.addEdge(2, 4);
    G.addEdge(3, 5);
    G.addEdge(4, 5);
    G.addEdge(5, 6);
    G.addEdge(5, 7);
    G.addEdge(6, 8);
    G.addEdge(7, 8);
    G.addEdge(8, 9);
    G.addEdge(8, 10);
    Dijkstra sssp(G, 0);
    EXPECT_NO_THROW(sssp.run());


    std::set<std::vector<node>> paths = sssp.getPaths(9);
    ASSERT_EQ(paths.size(),4);
    std::vector<node> path1 {0,2,3,5,6,8,9};
    std::vector<node> path2 {0,2,3,5,7,8,9};
    std::vector<node> path3 {0,2,4,5,6,8,9};
    std::vector<node> path4 {0,2,4,5,7,8,9};
    std::vector<std::vector<node>> results {path1, path2, path3, path4};
    count i = 0;
    for (auto path : paths) {
        ASSERT_EQ(path,results[i]);
        i++;
    }
}

TEST_F(SSSPGTest, testDirectedBFS) {
/* Graph:
         ________
        /        \.
       0     3.    6
        \. ./  \ ./
          2     .5
        ./  \. / \.
       1     4    7
*/
    int n = 8;
    // G directed unweighted
    Graph G(n, false, true);

    G.addEdge(0, 6);
    G.addEdge(0, 2);
    G.addEdge(3, 2);
    G.addEdge(5, 3);
    G.addEdge(6, 5);
    G.addEdge(5, 7);
    G.addEdge(4, 5);
    G.addEdge(2, 4);
    G.addEdge(2, 1);


    BFS sssp(G, 0);
    sssp.run();
    EXPECT_EQ(sssp.distance(0), 0);
    EXPECT_EQ(sssp.distance(1), 2);
    EXPECT_EQ(sssp.distance(2), 1);
    EXPECT_EQ(sssp.distance(3), 3);
    EXPECT_EQ(sssp.distance(4), 2);
    EXPECT_EQ(sssp.distance(5), 2);
    EXPECT_EQ(sssp.distance(6), 1);
    EXPECT_EQ(sssp.distance(7), 3);
}

TEST_F(SSSPGTest, testDirectedDijkstra) {
/* Graph:
         ________
        /        \.
       0     3.    6
        \. ./  \ ./
          2     .5
        ./  \. / \.
       1     4    7
*/
    int n = 8;
    // G directed unweighted
    Graph G(n, false, true);

    G.addEdge(0, 6, 1);
    G.addEdge(0, 2, 1);
    G.addEdge(3, 2, 1);
    G.addEdge(5, 3, 1);
    G.addEdge(6, 5, 1);
    G.addEdge(5, 7, 1);
    G.addEdge(4, 5, 1);
    G.addEdge(2, 4, 1);
    G.addEdge(2, 1, 1);


    Dijkstra sssp(G, 0);
    sssp.run();
    EXPECT_EQ(sssp.distance(0), 0);
    EXPECT_EQ(sssp.distance(1), 2);
    EXPECT_EQ(sssp.distance(2), 1);
    EXPECT_EQ(sssp.distance(3), 3);
    EXPECT_EQ(sssp.distance(4), 2);
    EXPECT_EQ(sssp.distance(5), 2);
    EXPECT_EQ(sssp.distance(6), 1);
    EXPECT_EQ(sssp.distance(7), 3);
}
}
