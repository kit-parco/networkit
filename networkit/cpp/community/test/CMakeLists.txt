networkit_add_test(community CommunityGTest
    auxiliary generators io overlap scd structures)

networkit_add_benchmark(community CommunityDetectionBenchmark
        auxiliary centrality io structures)

networkit_add_test(community LouvainMapEquationGTest
		auxiliary generators io)

networkit_add_benchmark(community LouvainMapEquationBenchmark
		auxiliary generators io)
