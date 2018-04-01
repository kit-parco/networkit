/*
 * MahmoodyGroupBetweennessGTest.cpp
 *
 *  Created on:23.03.2018 
 *  Author:Marvin Pogoda 
 */

#include "MahmoodyGroupBetweennessGTest.h"
#include "../MahmoodyGroupBetweenness.h"
#include "../../auxiliary/Random.h"


namespace NetworKit {



TEST_F(MahmoodyGroupBetweennessGTest, testMahmoodyGroupBetweennessSmallGraph1) {
/* Graph:
0    3   6
\  / \ /
 2    5
/  \ / \
1    4   7
*/
	Aux::Random::setSeed(2,false);

	Graph g(8,false,false);

	g.addEdge(0, 2);
	g.addEdge(1, 2);
	g.addEdge(2, 3);
	g.addEdge(2, 4);
	g.addEdge(3, 5);
	g.addEdge(4, 5);
	g.addEdge(5, 6);
	g.addEdge(5, 7);

	MahmoodyGroupBetweenness gb(g,5,0.1);
	gb.run();
	

	EXPECT_EQ(gb.groupMaxBetweenness()[0],2);
	EXPECT_EQ(gb.groupMaxBetweenness()[1],5);
	EXPECT_EQ(gb.groupMaxBetweenness()[2],6);
	EXPECT_EQ(gb.groupMaxBetweenness()[3],7);
	EXPECT_EQ(gb.groupMaxBetweenness()[4],4);
	




}


TEST_F(MahmoodyGroupBetweennessGTest, testMahmoodyGroupBetweennessSmallGraph2) {

	Aux::Random::setSeed(2,false);
	Graph g(8, false,false);
	g.addEdge(0,1);
	g.addEdge(0,2);
	g.addEdge(0,3);
	g.addEdge(0,7);
	g.addEdge(1,5);
	g.addEdge(2,5);
	g.addEdge(3,5);
	g.addEdge(3,4);
	g.addEdge(5,6);
	g.addEdge(4,5);
	g.addEdge(4,6);

	MahmoodyGroupBetweenness gb(g,5,0.1);
	gb.run();

	EXPECT_EQ(gb.groupMaxBetweenness()[0],5);
	EXPECT_EQ(gb.groupMaxBetweenness()[1],0);
	EXPECT_EQ(gb.groupMaxBetweenness()[2],4);
	EXPECT_EQ(gb.groupMaxBetweenness()[3],6);
	EXPECT_EQ(gb.groupMaxBetweenness()[4],3);
	





}




TEST_F(MahmoodyGroupBetweennessGTest, testMahmoodyGroupBetweennessSmallGraph3) {

	Aux::Random::setSeed(13,false);
	Graph g(10, false,false);
	g.addEdge(2,4);
	g.addEdge(0,2);
	g.addEdge(1,8);
	g.addEdge(0,9);
	g.addEdge(2,8);
	g.addEdge(3,9);
	g.addEdge(1,5);
	g.addEdge(3,4);
	g.addEdge(1,2);
	g.addEdge(4,7);
	g.addEdge(6,7);

	MahmoodyGroupBetweenness gb(g,3,0.2);
	gb.run();

	EXPECT_EQ(gb.groupMaxBetweenness()[0],2);
	EXPECT_EQ(gb.groupMaxBetweenness()[1],3);
	EXPECT_EQ(gb.groupMaxBetweenness()[2],7);

	





}








} /* namespace NetworKit */
