// no-networkit-format
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/Log.hpp>

 // OpenMP
#include <omp.h>



void Aux::setNumberOfThreads(int nThreads) {
#ifdef _OPENMP
        omp_set_num_threads(nThreads);
#else
        ERROR("Thread option ignored since OpenMP is deactivated.");
#endif
}


int Aux::getCurrentNumberOfThreads() {
    #ifdef _OPENMP
        return omp_get_num_threads();
    #else
        return 1;
        ERROR("OpenMP is not available");
    #endif
}


int Aux::getMaxNumberOfThreads() {
    #ifdef _OPENMP
        return omp_get_max_threads();
    #else
        return 1;
        ERROR("OpenMP is not available");
    #endif
}
