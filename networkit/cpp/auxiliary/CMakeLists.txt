networkit_add_module(auxiliary
    BloomFilter.cpp
    BucketPQ.cpp
    Log.cpp
    NumericTools.cpp
    Parallelism.cpp
    Random.cpp
    SignalHandling.cpp
	ParallelTimings.cpp
    ParseString.cpp
    Timer.cpp
    UniformRandomSelector.cpp
    )

add_subdirectory(test)
