/**
 * This code has been downloaded from https://github.com/LEOAlg/HM and adapted
 *in order to integrate it inside Nornir. The behaviour of this algorithm is
 *described in the paper "A Probabilistic Graphical Model-based Approach for
 *Minimizing Energy Under Performance Constraints" by Mishra, Nikita and Zhang,
 *Huazhe and Lafferty, John D. and Hoffmann, Henry.
 **/

#ifdef ENABLE_ARMADILLO
#include <nornir/leo.h>

using namespace std;
using namespace arma;

#include <cfloat>

namespace leo {

#ifndef typedef_emParam_t
#define typedef_emParam_t
typedef struct {
  vec mu;
  mat C;
  mat T;
  double sigma;
} emParam_t;
#endif // typedef_emParam_t

#ifndef typedef_emValue_t
#define typedef_emValue_t
typedef struct {
  vec w;
  emParam_t em_t;
} emReturn_t;
#endif // typedef_emValue_t

void init_EMParam(emParam_t *Old, int n);
void EM(emParam_t *Old, mat *W, mat *y_em, int i, int ep, int iteration_limit,
        emReturn_t *Appl);
double residualAccuracy(vec A, vec B);
long Time_Difference(struct timeval end, struct timeval start);

void init_EMParam(emParam_t *Old, int n) {
  Old->mu = zeros<vec>(n);
  Old->C = eye<mat>(n, n);
  Old->T = zeros<mat>(n, n);
  Old->sigma = 1;
}

void EM(emParam_t *Old, mat *W, mat *y_em, int i, int ep, int iteration_limit,
        emReturn_t *Appl) {

  int n = y_em->n_rows;
  int m = y_em->n_cols;
  int numSamples = sum(sum(*W));
  double pi = 1;
  double tau = 1;
  // double error = INFINITY;
  mat I = diagmat(*W);
  double sigma = Old->sigma;
  vec mu = Old->mu;
  mat C = Old->C;
  mat wl = zeros<mat>(n, m);

  for (int iterator = 1; iterator <= iteration_limit; iterator++) {
    mat Cinv = inv(C);
    mat Cl0 = inv(I / sigma + Cinv);
    mat Cll = inv(eye<mat>(n, n) / sigma + Cinv);

    wl = (Cll / sigma) * (*y_em) + repmat(Cinv * mu, 1, m);
    wl.col(0) = Cl0 * ((y_em->col(0)) / sigma + Cinv * mu);

    double normSum =
        pow(norm(I * (y_em->col(0) - wl.col(0)), 2), 2) + trace(I * Cl0);
    vec wlSum = sum(wl, 1);
    mat ClSum = eye<mat>(n, n) + Cl0 + (m - 1) * Cll;
    mat wlSumCov = eye<mat>(n, n) + (wl.col(0) - mu) * trans(wl.col(0) - mu);
    for (int l = 1; l < m; l++) {
      wlSumCov = wlSumCov + (wl.col(l) - mu) * trans(wl.col(l) - mu);
      normSum =
          normSum + pow(norm(y_em->col(l) - wl.col(l), 2), 2) + trace(Cll);
    }
    mu = ((1.0 / (double) (pi + m)) * wlSum);
    C = ((1.0 / (double) (tau + m)) *
         (pi * (mu * trans(mu)) + tau * eye<mat>(n, n) + ClSum + wlSumCov));
    sigma = (1.0 / (double) ((m - 1) * n + numSamples)) * normSum;
    // error = norm(Old->mu - mu,"fro")+ norm(Old->C - C,"fro") + abs(Old->sigma
    // - sigma);
    Old->mu = mu;
    Old->C = C;
    Old->sigma = sigma;
  }

  Appl->w = wl.col(0);
  Appl->em_t = *Old;
}

long Time_Difference(struct timeval end, struct timeval start) {
  long mtime, seconds, useconds;
  seconds = end.tv_sec - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;

  mtime = ((seconds) *1000 + useconds / 1000.0) + 0.5;

  return mtime;
}

double residualAccuracy(vec A, vec B) {
  double n = A.n_elem;
  double rss = pow(norm(A - B, 2), 2);
  double tss = pow(norm(A - mean(A), 2), 2);
  double residualsquare = 1 - rss / tss;
  return (1 - (1 - residualsquare) * (n - 1) / (n - 2)) * 100.0;
}

double minNonZero(const mat *data) {
  double min = DBL_MAX;
  for (auto &val : *data) {
    if (val && val < min) {
      min = val;
    }
  }
  return min;
}

mat normalize(const mat *data, double low = 0, double high = 100.0) {
  double min = minNonZero(data);
  double max = data->max();

  double rng = max - min;
  mat copy = *data;

  for (auto &val : copy) {
    if (val) {
      val = high - (((high - low) * (max - val)) / rng);
    }
  }
  return copy;
}

vec denormalize(vec *data, mat *origData, double low = 0, double high = 100.0) {
  double min = minNonZero(origData);
  double max = origData->max();
  double rng = max - min;
  vec copy = *data;

  for (auto &val : copy) {
    val = max + (((val - high) * rng) / (high - low));
  }
  return copy;
}

void loadData(uint appId, mat *data, vec *trueData, mat *W,
              std::string dataFile, const vec *sampledData,
              bool perColumnNormalization, bool computeError = false) {
  int n; //, numSamples = 0;
  mat SupPower, SupPerf;

  // numSamples = 20;

  data->load(dataFile.c_str());
  n = data->n_rows; // # CONFIGURATIONS
  *W = zeros<vec>(n);

  for (size_t i = 0; i < sampledData->size(); i++) {
    if (sampledData->at(i) != 0) {
      //++numSamples;
      (*W)(i) = 1;
    }
  }

  if (computeError) {
    if (perColumnNormalization) {
      mat tmp(data->col(appId));
      *trueData = normalize(&tmp);
    } else {
      *trueData = normalize(data).col(appId);
    }
  }

  data->col(appId) = *sampledData;
  data->swap_cols(0, appId);
}

PredictionResults compute(uint appId, const std::string &dataFile,
                          const arma::vec *sampledData,
                          bool perColumnNormalization, bool computeError) {
  int n;
  mat data, normData, y_em, W;
  vec trueData;
  emParam_t oldData;
  emReturn_t applData;

  // LOAD DATA AND MISSING VALUES, BS is the name of target application
  loadData(appId, &data, &trueData, &W, dataFile, sampledData,
           perColumnNormalization, computeError);

  // Normalization
  if (perColumnNormalization) {
    for (size_t i = 0; i < data.n_cols; i++) {
      mat tmp(data.col(i));
      normData.insert_cols(i, normalize(&tmp));
    }
  } else {
    normData = normalize(&data);
  }

  n = data.n_rows; // # CONFIGURATIONS
  // m = data.n_cols;    // # APPLICATIONS

  init_EMParam(&oldData, n);
  EM(&oldData, &W, &normData, 0, 10000, 5, &applData);

  PredictionResults pr;
  if (computeError) {
    pr.accuracy = residualAccuracy(applData.w, trueData);
  } else {
    pr.accuracy = -1;
  }

  if (perColumnNormalization) {
    mat tmp(data.col(0));
    pr.predictions = denormalize(&(applData.w), &tmp);
  } else {
    pr.predictions = denormalize(&(applData.w), &data);
  }
  return pr;
}
} // namespace leo

#ifdef TESTMAIN
int main(int argc, char **argv) {
  int appId = atoi(argv[1]);
  char *powerFile = argv[2];
  char *bandwidthFile = argv[3];

  vec sampledPower, sampledPerformance;
  mat powerData, perfData;
  powerData.load(powerFile);
  perfData.load(bandwidthFile);
  sampledPower = powerData.col(appId);
  sampledPerformance = perfData.col(appId);

  size_t n = 312; // ATTENTION: This must be set equal to the number of columns
                  // in file (TODO)
  for (size_t i = 0; i < n; i++) {
    if (i % ((n / 20 + 1)) == 0) {
      sampledPower[i] = powerData.col(appId)[i];
      sampledPerformance[i] = perfData.col(appId)[i];
    } else {
      sampledPower[i] = 0;
      sampledPerformance[i] = 0;
    }
  }

  leo::PredictionResults pr =
      leo::compute(appId, bandwidthFile, &sampledPerformance, true, true);
  std::cout << pr.accuracy << std::endl;

  pr = leo::compute(appId, powerFile, &sampledPower, false, true);
  std::cout << pr.accuracy << std::endl;
}
#endif

#endif
