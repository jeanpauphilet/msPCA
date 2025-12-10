#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

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

// Compute the value x^T A x -- x is a vector
double evaluate(const Eigen::VectorXd& x, const Eigen::MatrixXd& A)
{
  return (x.transpose() * A * x)[0];
}
// Compute the value tr(x^T A x) -- x is a matrix with r vectors (columns)
double evaluate(const Eigen::MatrixXd& x, const Eigen::MatrixXd& A)
{
  return (x.transpose() * A * x).trace();
}


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

// Trucation operator: Takes a vector origlist, keeps only the k largest coordinates (in absolute value), and normalizes the vector
Eigen::VectorXd truncateVector(const Eigen::VectorXd& origlist, int k)
{
  Eigen::VectorXd list = origlist;
  Eigen::VectorXd kparse = Eigen::VectorXd::Zero(list.size());

  Rcpp::NumericVector newIndices = selectperm2(list, k);
  for (auto index : newIndices) {
    kparse[index] = origlist[index];
  }
  kparse.normalize(); // Normalization in place
  return kparse;
}

// Iterative truncation heuristic from Yuan and Zhang (2013): Looking for a sparse fixed point of x = Sigma*x
Eigen::VectorXd iterativeTruncHeuristic(int k, const Eigen::VectorXd& beta0, const Eigen::MatrixXd& prob_Sigma)
{

  Eigen::VectorXd beta = truncateVector(beta0, k);
  for (int i = 0; i < 100; i++)
  {
    // beta = truncateVector(prob_Sigma * beta, k); 
    beta = prob_Sigma * beta; 
    beta = truncateVector(beta, k); 
  }
  return beta;
}

// Inner routine: sPCA heuristic for a single PC case: Truncated Power Method of Yuan and Zhang (2013) with random restarts 
void singlePCHeuristic(int k, const Eigen::MatrixXd& prob_Sigma, const Eigen::VectorXd& beta0, double& lambda_partial, Eigen::VectorXd& x_output, int maxIter = 100, int timeLimit = 20)
{
  int n = prob_Sigma.rows();

  // Applies the iterative truncation heuristic starting from beta0
  Eigen::VectorXd bestBeta = iterativeTruncHeuristic(k, beta0, prob_Sigma);
  double bestObj = evaluate(bestBeta, prob_Sigma);

  // Applies the iterative truncation heuristic starting from random points
  int countdown = maxIter;
  time_t start = time(0);
  while (countdown > 0 && difftime(time(0), start) < timeLimit)
  {
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(n);
    for (auto i = 0; i < beta.rows(); i++) {
      beta(i) = R::rnorm(0, 1);
      // for (auto j = 0; j < beta.cols(); j++) {
      //   beta(i, j) = R::rnorm(0, 1);
      // }
    }
    beta = beta / beta.norm();

    beta = iterativeTruncHeuristic(k, beta, prob_Sigma);
    double obj = evaluate(beta, prob_Sigma);

    if (obj > bestObj)
    {
      bestObj = obj;
      bestBeta = beta;
      countdown = maxIter; // Reset the countdown if found a better solution via random restart
    }
    countdown--;
  }

  lambda_partial = bestObj;
  x_output = bestBeta;
  return;
}

// Computes the orthogonality violation of a family of r vectors x, defined as |x^T x - I_r|
double fnviolation(const Eigen::MatrixXd& x,
                    const Eigen::MatrixXd& Sigma,
                    int feasibilityConstraintType)
{
  const int r = x.cols();
  double v = 0;
  Eigen::MatrixXd y = x.transpose() * x;
  // for (size_t i = 0; i < r; i++) {
  //   for (size_t j = 0; j < r; j++) {
  //     if (i == j) {
  //       v += std::fabs(y(i, j) - 1);
  //     }
  //     else {
  //       v += std::fabs(y(i, j));
  //     }
  //   }
  // }
  // return v;

  if (feasibilityConstraintType == 0) {
    // Orthogonality: sum of |x_i^T x_j - delta_ij| 
    for (int i = 0; i < r; ++i) {
      for (int j = i; j < r; ++j) {
        if (i == j) {
          v += std::fabs(y(i, j) - 1.0);
        } else {
          v += std::fabs(y(i, j));
        }
      }
    }
  } else {
    // Uncorrelatedness: sum of |x_i^T Sigma x_j| for i > j and ||x_i||^2 - 1 for i == j
    Eigen::MatrixXd C = x.transpose() * Sigma * x; // Xᵀ Σ X

    for (int i = 0; i < r; ++i) {
      for (int j = i; j < r; ++j) {
        if (i == j) {
          v += std::fabs(y(i, j) - 1.0);
        } else {
          v += std::fabs(C(i, j));
        }
      }
    }
  }
  return v;
}

// Main function: Iterative deflation heuristic for sparse PCA with multiple PCs
// [[Rcpp::export]]
List iterativeDeflationHeuristic(
    Eigen::MatrixXd Sigma,
    int r,
    Rcpp::NumericVector ks, // size r
    int maxIter = 200,
    bool verbose = true,
    int feasibilityConstraintType = 0, // NEW v2: 0 = orthogonality constraints, 1 = uncorrelatedness constraints
    double violationTolerance = 1e-4,
    double stallingTolerance = 1e-8,
    int maxIterTPW = 200, 
    int timeLimitTPW = 20)
{
  int n = Sigma.rows();

  if (ks.size() < r) //If sparsity pattern is not fully defined, we complete it with zeros//we reduce the number of PCs to the number of sparsity patterns
  {
    warning("Warning: You requested %i PCs but only provided %i sparsity levels. Ran the algorithm for %i PCs instead.", r, ks.size(), ks.size());
    r = ks.size();
  }

  double ofv_best = -1e10; // Objective value of the best solution found (solution = set of r PCs)
  double violation_best = n; // Orthogonality violation of the best solution found 
  Eigen::MatrixXd x_best = Eigen::MatrixXd::Zero(n, r); // Best solution found

  double ofv_secondbest = -1e10; // Objective value of the second best solution found (solution = set of r PCs)
  Eigen::MatrixXd x_secondbest = Eigen::MatrixXd::Zero(n, r); // Second best solution found

  Eigen::MatrixXd x_current = Eigen::MatrixXd::Zero(n, r); // Current solution (solution = set of r PCs)
  double ofv_prev = 0; // Objective value of the previous solution
  double ofv_overall = 0; // Objective value of the current solution 

  double theLambda = 0; // Penalty parameter on the orthogonality constraint

  Eigen::VectorXd weights = Eigen::VectorXd::Zero(r); // Weights assigned to each PC in the penalization heuristic (initialized through the first iteration of the algorithm)
  
  double stepSize = 0;
  int slowPeriod = ceil(ConstantArguments::slowPeriodRate * maxIter); // Slow phase: smallest step size (for the penalty) in the beginning to encourage exploration
  int fastPeriod = ceil(ConstantArguments::fastPeriodRate * maxIter); // Faster phase: highest step size (for the penalty) in the end to encourage feasibility

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


    Rcout.width(ConstantArguments::separatorLengthLong + ConstantArguments::wordLengthShort);
    Rcout << "Iteration |";

    Rcout.width(ConstantArguments::separatorLengthLong + ConstantArguments::wordLengthMiddle);
    Rcout << "Objective value |";


    Rcout.width(ConstantArguments::separatorLengthLong + ConstantArguments::wordLengthLong);
    Rcout << "Feasibility Violation |";

    Rcout.width(ConstantArguments::separatorLengthShort + ConstantArguments::wordLengthShort);
    Rcout << "Time";
    Rcout << endl;
  }

  auto startTime = std::chrono::high_resolution_clock::now();
  Eigen::MatrixXd sigma_current; // For memory: Current deflated covariance matrix
  Eigen::VectorXd x_output; // For memory: Current PC
  double lambda_partial = 0; // For memory: Fraction of the variance explained by the current PC 

  for (int theIter = 1; theIter <= maxIter; theIter++)
  {
    theLambda += stepSize;
    
    // Iteratively updating each PC
    for (int t = 0; t < r; t++)
    {
      if (theIter == 1 && ks[t] > n) // Verifies that the sparsity level is not higher than the dimension
      {
        warning("Warning: For PC #%i, you requested a sparsity level (%i) that exceeds the dimension (%i). Ran the algorithm with a sparsity level of %i instead.", t+1, ks[t], n, n);
        ks[t] = n;
      }

      // Eigen::MatrixXd 
      sigma_current = Sigma;
      for (int s = 0; s < r; s++)
      {
        if (s != t)
        {
          // sigma_current -= theLambda * weights[s] * x_current.col(s) * x_current.col(s).transpose();

          Eigen::VectorXd w;
          if (feasibilityConstraintType == 0) {
            // Orthogonality: penalty matrix is \sum_s u_s u_s^\top
            w = x_current.col(s);
          } else {
            // Uncorrelatedness: penalty matrix is \sum_s (Σ u_s)(Σ u_s)ᵀ
            w = Sigma * x_current.col(s);
          }
          sigma_current -= theLambda * weights[s] * (w * w.transpose());
        }
      }

      // Ensure sigma_current is symmetric (increase numerical accuracy)
      sigma_current = (sigma_current + sigma_current.transpose()) / 2;
      
      // Make sigma_current PSD
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(sigma_current);
      double lambda0 = -solver.eigenvalues().minCoeff() + 1e-4;
      if (lambda0 > 0) 
      {
        for (int i = 0; i < sigma_current.rows(); i++)
        {
          sigma_current(i, i) += lambda0;
        }
      }

      // Compute largest eigenvector of sigma_current, to start the TPW heuristic
      int index;
      solver.eigenvalues().maxCoeff(&index);
      Eigen::VectorXd beta0 = solver.eigenvectors().col(index); 

      singlePCHeuristic(ks[t], sigma_current, beta0, lambda_partial, x_output, maxIterTPW);

      x_current.col(t) = x_output;

      if (theIter == 1) // Initialize the weights on each PC at the first iteration
      {
        weights[t] = lambda_partial;
      }
    }

    ofv_prev = ofv_overall;
    ofv_overall = evaluate(x_current, Sigma);

    // if (theIter == 1) // TBD if needed or if initialization with -1e10 is enough
    // {
    //   ofv_prev = ofv_overall;
    // }

    double violation = fnviolation(x_current, Sigma, feasibilityConstraintType);
    double violation_forStep = violation; // For the step size: if truncate values that are too small for numerical stability
    if (1e-7 > violation) {
      violation_forStep = 1e-7;
    }

    stepSize = (theIter < fastPeriod ? ConstantArguments::changedRateLow : ConstantArguments::changedRateHigh) * (theIter < slowPeriod ? violation_forStep : ofv_overall / violation_forStep);
    
    auto stopTime = chrono::high_resolution_clock::now();
    chrono::milliseconds executionTime = chrono::duration_cast<chrono::milliseconds>(stopTime - startTime);
    
    bool stopCriterion = (theIter == maxIter) || (std::fabs(ofv_prev - ofv_overall) < stallingTolerance && violation < violationTolerance);
    if (verbose)
    {
      if (maxIter <= 25 || theIter % 10 == 1 || stopCriterion) // Display at every iteration if less than 25 iterations, or every 10 iterations otherwise, or at the last iteration
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

    if (violation < violationTolerance || (theIter == maxIter && ofv_best < 0)) //If current solution is feasible (within tolerance) or if we reached the last iteration and no feasible solution was found (ofv_best still <0)
    {
      // double ofv_current = (x_current.transpose() * Sigma * x_current).trace();
      if (ofv_best < ofv_overall)
      {
        x_secondbest = x_best; // Saving second-best solution
        ofv_secondbest = ofv_best;

        x_best = x_current; // Saving best solution
        ofv_best = ofv_overall;
      }
      if (ofv_best > ofv_overall && ofv_overall > ofv_secondbest)
      {
        x_secondbest = x_current; // Saving second-best solution
        ofv_secondbest = ofv_overall;
      }
    }

    if (std::fabs(ofv_prev - ofv_overall) < stallingTolerance && violation < violationTolerance) //If the algorithm is stalling (in terms of objective value) and the current solution is feasible (within tolerance)
    {
      if (ofv_best < 0) //Safety check: if no feasible solution was found yet, we take the current solution as the best solution
      {
        x_best = x_current;
        ofv_best = ofv_overall; //(x_current.transpose() * Sigma * x_current).trace();
      }
      break; //Stop the algorithm
    }
  }

  auto stopTime = std::chrono::high_resolution_clock::now();
  std::chrono::milliseconds allExecutionTime = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
  double runtime = (double)allExecutionTime.count() / ConstantArguments::millisecondsToSeconds;
  violation_best = fnviolation(x_best, Sigma, feasibilityConstraintType);
  ofv_best = evaluate(x_best, Sigma);
  List result = List::create(Named("objective_value") = ofv_best,
                             Named("feasibility_violation") = violation_best,
                             Named("runtime") = runtime,
                             Named("x_best") = x_best);
  return result;
}

// Main function: Truncated Power Method for sparse PCA with 1 PC
// [[Rcpp::export]]
List truncatedPowerMethod(
    Eigen::MatrixXd Sigma,
    int k, // Sparsity level
    int maxIter = 200,
    bool verbose = true,
    int timeLimit = 10)
{
  int n = Sigma.rows();

  if (k > n) //If sparsity level is higher than dimension, we reduce the sparsity level to the dimension
  {
    warning("You requested a sparsity level (%i) that exceeds the dimension (%i). Ran the algorithm with a sparsity level of %i instead.", k, n, n);
    k = n;
  }
  auto startTime = std::chrono::high_resolution_clock::now();

  Eigen::VectorXd x_output; // For memory: Current PC
  double lambda_partial = 0; // For memory: Fraction of the variance explained by the current PC 

  // Compute largest eigenvector of Sigma
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(Sigma);
  int index;
  solver.eigenvalues().maxCoeff(&index);
  Eigen::VectorXd beta0 = solver.eigenvectors().col(index); 

  singlePCHeuristic(k, Sigma, beta0, lambda_partial, x_output, maxIter, timeLimit);

  auto stopTime = std::chrono::high_resolution_clock::now();
  std::chrono::milliseconds allExecutionTime = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime);
  
  double runtime = (double)allExecutionTime.count() / ConstantArguments::millisecondsToSeconds;

  Eigen::MatrixXd x_best = Eigen::MatrixXd::Zero(n, 1); // Best solution found
  x_best.col(0) = x_output;

  List result = List::create(Named("objective_value") = lambda_partial,
                             Named("runtime") = runtime,
                             Named("x_best") = x_best);
  return result;

}