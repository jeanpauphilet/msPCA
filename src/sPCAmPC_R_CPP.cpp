#include <RcppEigen.h>
#include <Rcpp.h>
#include <chrono>
#include "ConstantArguments.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]

Eigen::MatrixXd prob_data;
Eigen::MatrixXd prob_Sigma;
double ofv_best;
double violation_best;
double lambda_partial;
Eigen::VectorXd x_output;

double absoluteDouble(double aNumber) {
  if (aNumber >= 0) {
    return aNumber;
  }
  return -aNumber;
}

int find_the_index(Rcpp::NumericVector aVector, double aNumber, const Rcpp::NumericVector& blocked)
{
  for (int i = 0; i < aVector.size(); i++) {
    if (absoluteDouble(aVector[i] - aNumber) < 1e-8) {
      int j = 0;
      for (; j < blocked.size(); j++) {
        if (absoluteDouble(blocked[i] - static_cast<double>(i)) < 1e-8) {
          break;
        }
      }
      if (j == blocked.size()) {
        return i;
      }
    }
  }
  return 0;
}

double evaluate(const Eigen::VectorXd& solution)
{
  return (solution.transpose() * prob_Sigma * solution)[0];
}

Rcpp::NumericVector selectperm2(const Eigen::VectorXd& x, int k)
{
  Rcpp::NumericVector numbers{};
  for (int i = 0; i < x.size(); i++) {
    numbers.push_back(absoluteDouble(x(i)));
  }

  Rcpp::NumericVector originalNumbers{};
  for (int i = 0; i < numbers.size(); i++) {
    originalNumbers.push_back(numbers[i]);
  }

  numbers.sort();
  Rcpp::NumericVector indexes{};
  for (int i = 0, j = numbers.size() - 1; i < k; i++, j--) {
    int index = find_the_index(originalNumbers, numbers[j], indexes);
    indexes.push_back(index);
  }
  return indexes;
}

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

Eigen::VectorXd eigSubset(const Rcpp::NumericVector& support, int k, const Eigen::VectorXd& beta0)
{

  Eigen::VectorXd beta = Hk(beta0, k, support);
  for (int i = 0; i < 100; i++)
  {
    beta = Hk(prob_Sigma * beta, k, support);
  }
  return beta;
}

void subset(int k, int timeLimit, Rcpp::NumericVector& support, int countdown = 100)
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

  Eigen::VectorXd bestBeta = eigSubset(support, k, beta0);
  double bestObj = evaluate(bestBeta);
  time_t start = time(0);
  while (countdown > 0 && difftime(time(0), start) < timeLimit)
  {
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(n);
    for (auto i = 0; i < beta.rows(); i++) {
      for (auto j = 0; j < beta.cols(); j++) {
        beta(i, j) = R::rnorm(0, 1);
      }
    }
    beta = beta / beta0.norm();
    beta = eigSubset(support, k, beta);
    double obj = evaluate(beta);
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
Eigen::MatrixXd cpp_findmultPCs_deflation(
    Eigen::MatrixXd Sigma,
    int r,
    Rcpp::NumericVector ks, // size r
    int numIters = 200,
    bool verbose = true,
    double violation_tolerance = 1e-4)
{
  /*Rcpp::NumericVector ks;
  for (int i = 0; i < ks_input.rows(); i++) {
    for (int j = 0; j < ks_input.cols(); j++) {
      ks.push_back(ks_input(i, j));
    }
  }*/
  int n = Sigma.rows();
  ofv_best = -1e10;
  violation_best = n;
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
      prob_data = solver.operatorSqrt();
      prob_Sigma = sigma_current;

      Rcpp::NumericVector support;
      support.push_back(0);
      subset(ks[t], 20, support);

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
  violation_best = fnviolation(x_best);
  ofv_best = (x_best.transpose() * Sigma * x_best).trace();
  return x_best;
}
