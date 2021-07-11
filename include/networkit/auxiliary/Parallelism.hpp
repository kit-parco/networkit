// no-networkit-format
/*
 * Parallelism.h
 *
 * Control functions related to parallelism.
 *
 *  Created on: 13.11.2013
 *      Author: Christian Staudt (christian.staudt@kit.edu)
 */

#ifndef NETWORKIT_AUXILIARY_PARALLELISM_HPP_
#define NETWORKIT_AUXILIARY_PARALLELISM_HPP_


namespace Aux {


/**
 *  Set the number of threads available to the program.
 */
void setNumberOfThreads(int nThreads);

/**
 *
 * @return The number of threads currently running.
 */
int getCurrentNumberOfThreads();


/**
 *
 * @return The maximum number of threads available to the program.
 */
int getMaxNumberOfThreads();
}

#endif // NETWORKIT_AUXILIARY_PARALLELISM_HPP_
