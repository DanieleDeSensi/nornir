/*
 * runEM.h
 *
 * Code generation for function 'runEM4'
 *
 * C source code generated on: Thu Jun  5 13:30:30 2014
 *
 */

#ifndef __LEO_H__
#define __LEO_H__
/* Include files */

#include <iostream>
#include <math.h>
//#include <random>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <armadillo>
#include <list>
#include <omp.h>

#define EPSILON 70000
#define ITERATION_LIMIT 100
#define INF 1000000000
/* Function Declarations */

namespace leo{

typedef struct{
    arma::vec predictions;
    double accuracy; // Between 0 and 100
}PredictionResults;

/**
 * @param appId Between [0, N - 1] (N = number of applications in the file).
 * @param dataFile The name of the file containing the data profiles
 *        for all the configurations and for all the applications (not
 *        normalized).
 * @param sampledData A vector with a length equal to the number of
 *        possible configurations. The value in a specific position is 0 if the
 *        corresponding configuration was not sampled. If it is != 0, it is the
 *        data (power or bandwidth) in the corresponding sampled configuration.
 * @param computeError If true computes the model error (comparing with the case 
 *        where all the values for all the configurations are known.
 *        
 **/
PredictionResults compute(uint appId, std::string dataFile,
                          const arma::vec* sampledData,
                          bool perColumnNormalization,
                          bool computeError = false);

}

#endif

