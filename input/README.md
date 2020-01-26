## Input folder

This folder is used for input files needed in the NetworKit. Keep in mind that this location is not part of the gitignore file and changes are committed to the git.

## Sample graphs

This folder contains many smaller graphs. They are used in the testing suite and can be used by you as a starting point.

| Graph                            | Python format           | Cpp Reader              | Filesize (in byte) | Number of Nodes | Number of Edges | Directed | Weighted | Weblink |
|----------------------------------|-------------------------|-------------------------|--------------------|-----------------|-----------------|----------|----------|---------|
| 10_clusters.graph                | Format.METIS            | METISGraphReader        | 3                  | 100             | 556             | FALSE    | FALSE    | n/a     |
| airfoil1.gi                      | n/a                     | DibapGraphReader        | 149370             | 4253            | 12289           | FALSE    | FALSE    | n/a     |
| airfoil1.graph                   | Format.METIS            | METISGraphReader        | 116538             | 4253            | 12289           | n/a      | n/a      | n/a     |
| astro-ph.graph                   | Format.METIS            | METISGraphReader        | 1259548            | 16706           | 121251          | FALSE    | FALSE    | n/a     |
| caidaRouterLevel.graph           | Format.METIS            | METISGraphReader        | 7490360            | 192244          | 609066          | FALSE    | FALSE    | n/a     |
| celegans_metabolic.graph         | Format.METIS            | METISGraphReader        | 16015              | 453             | 2025            | FALSE    | FALSE    | n/a     |
| celegans_metabolic.thrill        | Format.THRILL           | ThrillGraphReader       | 8554               | 453             | 2025            | FALSE    | FALSE    | n/a     |
| comments.edgelist                | Format.EdgeListTabOne   | EdgeListReader          | 127                | 10              | 10              | FALSE    | FALSE    | n/a     |
| community_overlapping.cover      | n/a                     | CoverReader             | 61                 | n/a             | n/a             | n/a      | n/a      | n/a     |
| community_overlapping.dat        | n/a                     | EdgeListCoverReader     | 83                 | n/a             | n/a             | n/a      | n/a      | n/a     |
| community.dat                    | Format.EdgeListTabOne   | EdgeListPartitionReader | 61                 | 10              | 9               | FALSE    | FALSE    | n/a     |
| dynamicTest.gexf                 | GEXFReader              | n/a                     | 2677               | n/a             | n/a             | n/a      | n/a      | n/a     |
| dynamicTest2.gexf                | GEXFReader              | n/a                     | 74818              | n/a             | n/a             | n/a      | n/a      | n/a     |
| dynamicTest3.gexf                | GEXFReader              | n/a                     | 921                | n/a             | n/a             | n/a      | n/a      | n/a     |
| erdos_renyi_200_0.05.graph       | Format.METIS            | METISGraphReader        | 7                  | 200             | 1002            | FALSE    | FALSE    | n/a     |
| example.edgelist                 | Format.EdgeListTabOne   | EdgeListReader          | 85                 | 10              | 10              | FALSE    | FALSE    | n/a     |
| example.graph                    | Format.METIS            | METISGraphReader        | 16                 | 4               | 2               | FALSE    | FALSE    | n/a     |
| example2.dgs                     | n/a                     | DGSReader               | 161                | n/a             | n/a             | n/a      | n/a      | n/a     |
| foodweb-baydry.konect            | Format.KONECT           | KONECTGraphReader       | 42847              | 128             | 2137            | TRUE     | TRUE     | n/a     |
| foodweb-baydry.networkit         | Format.NETWORKITBINARY  | NetworkitBinaryReader   | 42491              | 128             | 2137            | TRUE     | TRUE     | n/a     |
| grid-5x5-dist-arch.graph         | Format.METIS            | METISGraphReader        | 218                | 25              | 40              | FALSE    | FALSE    | n/a     |
| hamming6-4.edgelist              | Format.EdgeListSpaceOne | EdgeListReader          | 4088               | 64              | 704             | FALSE    | FALSE    | n/a     |
| hep-th.graph                     | Format.METIS            | METISGraphReader        | 157814             | 8361            | 15751           | FALSE    | FALSE    | n/a     |
| jazz.graph                       | Format.METIS            | METISGraphReader        | 19445              | 198             | 2742            | FALSE    | FALSE    | n/a     |
| jazz2_directed.gml               | Format.GML              | GMLGraphReader          | 289                | 5               | 4               | TRUE     | FALSE    | n/a     |
| jazz2_undirected.gml             | Format.GML              | GMLGraphReader          | 276                | 5               | 4               | FALSE    | FALSE    | n/a     |
| jazz2double.graph                | Format.METIS            | METISGraphReader        | 58                 | 5               | 3               | FALSE    | TRUE     | n/a     |
| johnson8-4-4.edgelist            | Format.EdgeListSpaceOne | EdgeListReader          | 10716              | 70              | 1855            | FALSE    | FALSE    | n/a     |
| karate.graph                     | Format.METIS            | METISGraphReader        | 450                | 34              | 78              | FALSE    | FALSE    | n/a     |
| lesmis.graph                     | Format.METIS            | METISGraphReader        | 2630               | 77              | 254             | FALSE    | TRUE     | n/a     |
| lfr_small.graph                  | Format.METIS            | METISGraphReader        | 4                  | 200             | 595             | FALSE    | FALSE    | n/a     |
| looptest1.gml                    | Format.GML              | GMLGraphReader          | 676                | 9               | 12              | FALSE    | FALSE    | n/a     |
| looptest2.gml                    | Format.GML              | GMLGraphReader          | 754                | 9               | 14              | FALSE    | FALSE    | n/a     |
| MIT8.edgelist                    | Format.EdgeListTabZero  | EdgeListReader          | 2425755            | 6440            | 251252          | FALSE    | FALSE    | n/a     |
| network_overlapping.dat          | Format.EdgeListTabOne   | EdgeListReader          | 148                | 10              | 17              | FALSE    | FALSE    | n/a     |
| network.dat                      | Format.EdgeListTabOne   | EdgeListReader          | 96                 | 10              | 10              | FALSE    | FALSE    | n/a     |
| PGPgiantcompo.graph              | Format.METIS            | METISGraphReader        | 249406             | 10680           | 24316           | FALSE    | FALSE    | n/a     |
| polblogs.graph                   | Format.METIS            | METISGraphReader        | 143808             | 1490            | 16715           | FALSE    | FALSE    | n/a     |
| power.graph                      | Format.METIS            | METISGraphReader        | 67968              | 4941            | 6594            | FALSE    | FALSE    | n/a     |
| power.gt                         | Format.GraphToolBinary  | GraphToolBinaryReader   | 172599             | 4941            | 6594            | FALSE    | FALSE    | n/a     |
| spaceseparated_weighted.edgelist | Format.EdgeListSpaceOne | EdgeListReader          | 17                 | 3               | 3               | FALSE    | TRUE     | n/a     |
| spaceseparated.edgelist          | Format.EdgeListSpaceOne | EdgeListReader          | 86                 | 10              | 10              | FALSE    | FALSE    | n/a     |
| staticTest.gexf                  | GEXFReader              | n/a                     | 141043             | n/a             | n/a             | n/a      | n/a      | n/a     |
| tiny_01.graph                    | Format.METIS            | METISGraphReader        | 558                | 7               | 11              | FALSE    | FALSE    | n/a     |
| tiny_02.graph                    | Format.METIS            | METISGraphReader        | 995                | 7               | 11              | FALSE    | TRUE     | n/a     |
| tiny_03.graph                    | Format.METIS            | METISGraphReader        | 1047               | 7               | 11              | FALSE    | TRUE     | n/a     |
| tiny_04.graph                    | Format.METIS            | METISGraphReader        | 1226               | 7               | 11              | FALSE    | FALSE    | n/a     |
| wiki-Vote.txt                    | Format.SNAP             | SNAPGraphReader         | 1095061            | 7115            | 100762          | FALSE    | FALSE    | n/a     |
| wing.graph                       | Format.METIS            | METISGraphReader        | 1482391            | 62032           | 121544          | FALSE    | FALSE    | n/a     |
