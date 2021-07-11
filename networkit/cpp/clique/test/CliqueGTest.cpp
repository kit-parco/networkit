// no-networkit-format
/*
 * CliqueGTest.cpp
 *
 *  Created on: 08.12.2014
 *      Author: henningm
 */

#include <gtest/gtest.h>

#include <networkit/clique/MaximalCliques.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/io/EdgeListReader.hpp>

namespace NetworKit {

class CliqueGTest: public testing::Test {};

TEST_F(CliqueGTest, testMaxCliqueOnSmallerGraphs) {
    EdgeListReader r(' ',1,"%");
    Graph gJohnson = r.read("input/johnson8-4-4.edgelist");
    Graph gHamming = r.read("input/hamming6-4.edgelist");

    MaximalCliques mcJohnson(gJohnson,true);
    MaximalCliques mcHamming(gHamming,true);

    mcJohnson.run();
    mcHamming.run();

    std::vector<node> cliqueJohnson = mcJohnson.getCliques()[0];
    count maxCliqueSizeJohnson = cliqueJohnson.size();
    std::vector<node> cliqueHamming = mcHamming.getCliques()[0];
    count maxCliqueSizeHamming = cliqueHamming.size();

    EXPECT_EQ(14u,maxCliqueSizeJohnson) << "maximum clique size on graph johnson8-4-4 is not correct";
    EXPECT_EQ(4u,maxCliqueSizeHamming) << "maximum clique size on graph hamming6-4 is not correct";

    EXPECT_EQ(14u, cliqueJohnson.size());
    EXPECT_EQ(4u, cliqueHamming.size());
}

} /* namespace NetworKit */
