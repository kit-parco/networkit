networkit_add_module(community
    AdjustedRandMeasure.cpp
    ClusteringGenerator.cpp
    CommunityDetectionAlgorithm.cpp
    Conductance.cpp
    CoverHubDominance.cpp
    CoverF1Similarity.cpp
    OverlappingNMIDistance.cpp
    Coverage.cpp
    CutClustering.cpp
    DissimilarityMeasure.cpp
    DynamicNMIDistance.cpp
    EdgeCut.cpp
    GraphClusteringTools.cpp
    GraphStructuralRandMeasure.cpp
    HubDominance.cpp
    IntrapartitionDensity.cpp
    IsolatedInterpartitionConductance.cpp
    IsolatedInterpartitionExpansion.cpp
    JaccardMeasure.cpp
    LFM.cpp
    LPDegreeOrdered.cpp
    LocalCoverEvaluation.cpp
    LocalPartitionEvaluation.cpp
    Modularity.cpp
    NMIDistance.cpp
    NodeStructuralRandMeasure.cpp
    OverlappingCommunityDetectionAlgorithm.cpp
    PLM.cpp
    PLP.cpp
    ParallelAgglomerativeClusterer.cpp
    PartitionFragmentation.cpp
    PartitionHubDominance.cpp
    PartitionIntersection.cpp
    SampledGraphStructuralRandMeasure.cpp
    SampledNodeStructuralRandMeasure.cpp
    StablePartitionNodes.cpp
	LouvainMapEquation.cpp
    )

networkit_module_link_modules(community
    auxiliary base coarsening components flow graph matching scd structures)

add_subdirectory(test)

