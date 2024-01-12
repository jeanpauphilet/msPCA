#include <RcppEigen.h>
#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <vector>
#include "ConstantArguments.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]

// Absolute value of a double
double absoluteDouble(double aNumber) {
  if (aNumber >= 0) {
    return aNumber;
  }
  return -aNumber;
}

// Compute the value x^T Sigma x
double evaluate(const Eigen::VectorXd& x, const Eigen::MatrixXd& prob_Sigma)
{
  return (x.transpose() * prob_Sigma * x)[0];
}

/*Response:
 I was not familiar with using C++ in R, so I implemented my own sort function to test.
 I then forgot about this issue.
 I have rewritten this function.
 */
// Selects the k indices of x corresponding the k largest coordinates (in absolute value)
Rcpp::NumericVector selectperm2(const Eigen::VectorXd& x, int k)
{
  std::vector<int> indexes(x.size());
  std::iota(indexes.begin(), indexes.end(), 0);
  std::partial_sort(
    indexes.begin(), indexes.begin() + k, indexes.end(),
    [&](int i, int j) -> bool {
      return absoluteDouble(x(i)) > absoluteDouble(x(j));
    });
  indexes.resize(k);
  Rcpp::NumericVector numbers{};
  for (auto index : indexes) {
    numbers.push_back(index);
  }
  return numbers;
}

// Trucating heuristic: Takes a vector origlist, keeps only the k (k=sparsity) largest coordinates (in absolute value), and normalizes the vector
Eigen::VectorXd Hk(const Eigen::VectorXd& origlist, int sparsity, const Rcpp::NumericVector& support)
{
  Eigen::VectorXd list = origlist;
  Eigen::VectorXd kparse = Eigen::VectorXd::Zero(list.size());
  int nbIndicesToKeep = 0;
  for (int s : support)
  {
    nbIndicesToKeep += s == 1;
  }

  double dummyValue = list(0) - 1;
  for (size_t i = 0; i < list.size(); i++)
  {
    if (list(i) - 1 < dummyValue) {
      dummyValue = list(i) - 1;
    }
  }
  for (int i = 0; i < list.size(); i++)
  {
    if (support[i] > -1)
    {
      // kparse[i] = origlist[i]; // Question for Chenkai: I think this initialization is needed
      /* Response: to me, "support" is very mysterious in a way that under no circumstances will it be used.
       In other words, I did not believe "support" could contain any value other than "-1" in the program.
       I almost decided to delete all code related to "support",
       but I kept it because I guess it might have some meaning to you or future developers
       (e.g., as guidance for updates).
       Also, I may be wrong, but I also think "support" is a red herring in Julia code.*/
      list[i] = dummyValue;
    }
  }

  Rcpp::NumericVector newIndices = selectperm2(list, sparsity - nbIndicesToKeep);
  for (auto index : newIndices) {
    kparse[index] = origlist[index];
  }
  kparse = kparse / kparse.norm();
  return kparse;
}

Eigen::VectorXd eigSubset(const Rcpp::NumericVector& support, int k, const Eigen::VectorXd& beta0, const Eigen::MatrixXd& prob_Sigma)
{

  Eigen::VectorXd beta = Hk(beta0, k, support);
  for (int i = 0; i < 100; i++)
  {
    beta = Hk(prob_Sigma * beta, k, support); // Question for Chenkai: Is prob_Sigma * beta a proper matrix-vector multiplication?
    /*Response: I have no idea.
     It is from the original code, and it looks similar to me compared to the Julia code.
     If it is different from the Julia code, let's change it.*/
  }
  return beta;
}

// Note from Jean: support should be a vector of size n with +/- 1. If support[i] = 1, then the i-th coordinate needs to be included. Otherwise, it is free. Flexibility to be removed.
void subset(int k, int timeLimit, Rcpp::NumericVector& support, const Eigen::MatrixXd& prob_Sigma, double& lambda_partial, Eigen::VectorXd& x_output, int countdown = 100)
{
  int n = prob_Sigma.rows();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(prob_Sigma);
  int index;
  solver.eigenvalues().maxCoeff(&index);
  Eigen::VectorXd beta0 = solver.eigenvectors().col(index);

  if (support.size() == 1)
  {
    support[0] = -1;
    for (int i = 0; i < n - 1; i++) {
      support.push_back(-1);
    }
  }

  Eigen::VectorXd bestBeta = eigSubset(support, k, beta0, prob_Sigma);
  double bestObj = evaluate(bestBeta, prob_Sigma);
  time_t start = time(0);
  while (countdown > 0 && difftime(time(0), start) < timeLimit)
  {
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(n);
    for (auto i = 0; i < beta.rows(); i++) {
      for (auto j = 0; j < beta.cols(); j++) {
        beta(i, j) = R::rnorm(0, 1);
      }
    }
    // beta = beta / beta0.norm();
    beta = eigSubset(support, k, beta, prob_Sigma);
    double obj = evaluate(beta, prob_Sigma);
    if (obj > bestObj)
    {
      bestObj = obj;
      bestBeta = bestBeta;
      countdown = 100;
    }
    countdown--;
  }
  lambda_partial = bestObj;
  x_output = bestBeta;
  return;
}

// Computes the orthogonality violation of a family of r vectors x, defined as |x^T x - I_r|
double fnviolation(const Eigen::MatrixXd& x)
{
  double v = 0;
  Eigen::MatrixXd y = x.adjoint() * x;
  for (size_t i = 0; i < y.rows(); i++) {
    for (size_t j = 0; j < y.cols(); j++) {
      if (i == j) {
        v += double(fabs(y(i, j) - 1));
      }
      else {
        v += fabs(y(i, j));
      }
    }
  }
  return v;
}

// [[Rcpp::export]]
List cpp_findmultPCs_deflation(
    Eigen::MatrixXd Sigma,
    int r,
    Rcpp::NumericVector ks, // size r
    int numIters = 200,
    bool verbose = true,
    double violation_tolerance = 1e-4)
{
  int n = Sigma.rows();
  double ofv_best = -1e10;
  double violation_best = n;
  Eigen::MatrixXd x_best = Eigen::MatrixXd::Zero(n, r);
  Eigen::MatrixXd x_current = Eigen::MatrixXd::Zero(n, r);
  double ofv_prev = 0;
  double ofv_overall = 0;

  Eigen::VectorXd weights = Eigen::VectorXd::Zero(r);
  double theLambda = 0;
  double stepSize = 0;
  int slowPeriod = ceil(0.15 * numIters);
  int fastPeriod = ceil(0.75 * numIters);

  if (verbose)
  {
    Rcout << "---- Iterative deflation algorithm for sparse PCA with multiple PCs ---" << std::endl;
    Rcout << "Dimension: " << n << std::endl;
    Rcout << "Number of PCs: " << r << std::endl;
    Rcout << "Sparsity pattern: ";
    for (int t = 0; t < r; t++)
    {
      Rcout << " " << static_cast<int>(ks[t]);
    }
    Rcout << endl;
    Rcout << endl;


    Rcout.width(ConstantArguments::separatorLengthLong + ConstantArguments::wordLengthShort);
    Rcout << "Iteration |";

    Rcout.width(ConstantArguments::separatorLengthLong + ConstantArguments::wordLengthMiddle);
    Rcout << "Objective value |";


    Rcout.width(ConstantArguments::separatorLengthLong + ConstantArguments::wordLengthLong);
    Rcout << "Orthogonality Violation |";

    Rcout.width(ConstantArguments::separatorLengthShort + ConstantArguments::wordLengthShort);
    Rcout << "Time";
    Rcout << endl;
  }

  auto startTime = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd prob_Sigma;
  Eigen::VectorXd x_output;
  double lambda_partial = 0;
  for (int theIter = 1; theIter <= numIters; theIter++)
  {
    theLambda += stepSize;
    for (int t = 0; t < r; t++)
    {
      Eigen::MatrixXd sigma_current = Sigma;
      for (int s = 0; s < r; s++)
      {
        if (s != t)
        {
          sigma_current -= theLambda * weights[s] * x_current.col(s) * x_current.col(s).transpose();
        }
      }
      sigma_current = (sigma_current + sigma_current.transpose()) / 2;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(sigma_current);
      double lambda0 = -solver.eigenvalues().minCoeff() + 1e-4;
      for (int i = 0; i < sigma_current.rows(); i++)
      {
        sigma_current(i, i) += lambda0;
      }
      solver = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(sigma_current);
      prob_Sigma = sigma_current;

      Rcpp::NumericVector support;
      support.push_back(0);
      subset(ks[t], 20, support, prob_Sigma, lambda_partial, x_output);
      x_current.col(t) = x_output;

      if (theIter == 1)
      {
        weights[t] = lambda_partial;
      }
    }
    ofv_prev = ofv_overall;
    ofv_overall = (x_current.transpose() * Sigma * x_current).trace();
    if (theIter == 1)
    {
      ofv_prev = ofv_overall;
    }

    double violation = fnviolation(x_current);
    if (1e-7 > violation) {
      violation = 1e-7;
    }
    stepSize = (theIter < fastPeriod ? 0.01 : 0.05) * (theIter < slowPeriod ? violation : ofv_overall / violation);
    auto stopTime = chrono::high_resolution_clock::now();
    chrono::milliseconds executionTime = chrono::duration_cast<chrono::milliseconds>(stopTime - startTime);
    if (verbose)
    {
      if (numIters <= 25 || theIter % 10 == 0)
      {
        Rcout.width(ConstantArguments::separatorLengthShort + ConstantArguments::wordLengthShort);
        Rcout << theIter << " |";

        Rcout.width(ConstantArguments::separatorLengthShort + ConstantArguments::wordLengthMiddle);
        Rcout << setprecision(ConstantArguments::precisionForObjectiveValue) << ofv_overall / Sigma.trace() << " |";

        Rcout.width(ConstantArguments::separatorLengthShort + ConstantArguments::wordLengthLong);
        Rcout << scientific << setprecision(ConstantArguments::precisionForOrthogonalityViolation)
              << violation << " |" << defaultfloat;

        Rcout.width(ConstantArguments::separatorLengthShort + ConstantArguments::wordLengthShort);
        Rcout << fixed << setprecision(ConstantArguments::precisionForTime)
              << (double)executionTime.count() / ConstantArguments::millisecondsToSeconds << defaultfloat;
        Rcout << endl;
      }
    }

    if (violation < violation_tolerance || (theIter == numIters && ofv_best < 0))
    {
      double ofv_current = (x_current.transpose() * Sigma * x_current).trace();
      if (ofv_best < ofv_current)
      {
        x_best = x_current;
        ofv_best = ofv_current;
      }
    }

    if (fabs(ofv_prev - ofv_overall) < 1e-8 && violation < violation_tolerance)
    {
      if (ofv_best < 0) {
        x_best = x_current;
        ofv_best = (x_current.transpose() * Sigma * x_current).trace();
      }
      break;
    }
  }
  auto stopTime = std::chrono::high_resolution_clock::now();
  std::chrono::milliseconds allExecutionTime = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
  double runtime = (double)allExecutionTime.count() / ConstantArguments::millisecondsToSeconds;
  violation_best = fnviolation(x_best);
  ofv_best = (x_best.transpose() * Sigma * x_best).trace();
  List result = List::create(Named("ofv_best") = ofv_best,
                             Named("violation_best") = violation_best,
                             Named("runtime") = runtime,
                             Named("x_best") = x_best);
  return result;
}
