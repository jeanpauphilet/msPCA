#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#include <algorithm>
#include <chrono>
#include <functional>
#include <numeric>
#include <vector>
#include "ConstantArguments.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::depends(RcppEigen)]]

// Truncation operator: keeps the k largest coordinates (by absolute value) and normalizes.
// Uses nth_element (O(n) average) to avoid O(n log k) partial_sort and R-side allocation.
inline Eigen::VectorXd truncateVector(const Eigen::VectorXd& v, int k)
{
  const int n = static_cast<int>(v.size());
  if (k >= n) {
    Eigen::VectorXd y = v;
    y.normalize();
    return y;
  }
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::nth_element(idx.begin(), idx.begin() + k, idx.end(),
                   [&](int a, int b){ return std::abs(v[a]) > std::abs(v[b]); });
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < k; ++i) y[idx[i]] = v[idx[i]];
  y.normalize();
  return y;
}

// Compute the value x^T A x -- x is a vector
inline double evaluate(const Eigen::VectorXd& x, const Eigen::MatrixXd& A)
{
  return x.dot(A * x);
}
// Compute the value tr(x^T A x) -- x is a matrix with r vectors (columns)
double evaluate(const Eigen::MatrixXd& x, const Eigen::MatrixXd& A)
{
  return (x.transpose() * A * x).trace();
}

// Evaluate x^T M x using a matvec functor
inline double evaluateFunctor(const Eigen::VectorXd& x,
                               const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& applyM)
{
  return x.dot(applyM(x));
}

// Iterative truncation heuristic (Yuan & Zhang 2013) using a matvec functor.
// Includes early-exit convergence check to avoid unnecessary iterations.
Eigen::VectorXd iterativeTruncHeuristic(int k, const Eigen::VectorXd& beta0,
                                         const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& applyM)
{
  Eigen::VectorXd beta = truncateVector(beta0, k);
  for (int i = 0; i < 100; ++i)
  {
    Eigen::VectorXd prev = beta;
    beta = truncateVector(applyM(beta), k);
    if ((beta - prev).squaredNorm() < 1e-12) break;
  }
  return beta;
}


// Inner routine: sPCA heuristic for a single PC — functor version.
// outerIter gates the random restart budget: full budget on outer iteration 1,
// restartsAfterFirstIter restarts on all subsequent iterations.
void singlePCHeuristic(int k,
                       const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& applyM,
                       int n,
                       const Eigen::VectorXd& beta0,
                       double& lambda_partial,
                       Eigen::VectorXd& x_output,
                       int maxIterTPM = 10,
                       int timeLimitTPM = 20)
{
  Eigen::VectorXd bestBeta = iterativeTruncHeuristic(k, beta0, applyM);
  double bestObj = evaluateFunctor(bestBeta, applyM);

  int countdown = maxIterTPM;
  time_t start = time(0);
  while (countdown > 0 && difftime(time(0), start) < timeLimitTPM)
  {
    Eigen::VectorXd beta(n);
    for (int i = 0; i < n; ++i) beta(i) = R::rnorm(0, 1);
    beta.normalize();

    beta = iterativeTruncHeuristic(k, beta, applyM);
    double obj = evaluateFunctor(beta, applyM);

    if (obj > bestObj)
    {
      bestObj = obj;
      bestBeta = beta;
      countdown = maxIterTPM;
    }
    countdown--;
  }

  lambda_partial = bestObj;
  x_output = bestBeta;
}

// Matrix-based convenience overload for truncatedPowerMethod (single-PC, no deflation).
void singlePCHeuristic(int k, const Eigen::MatrixXd& prob_Sigma, const Eigen::VectorXd& beta0,
                       double& lambda_partial, Eigen::VectorXd& x_output,
                       int maxIter = 10, int timeLimit = 20)
{
  int n = prob_Sigma.rows();
  auto applyM = [&](const Eigen::VectorXd& v) -> Eigen::VectorXd {
    return prob_Sigma * v;
  };
  singlePCHeuristic(k, applyM, n, beta0, lambda_partial, x_output,
                      maxIter, timeLimit);
}

// Computes the orthogonality/uncorrelatedness violation of a family of r vectors x.
double fnviolation(const Eigen::MatrixXd& x,
                    const Eigen::MatrixXd& Sigma,
                    int feasibilityConstraintType)
{
  const int r = x.cols();
  double v = 0;
  Eigen::MatrixXd y = x.transpose() * x;

  if (feasibilityConstraintType == 0) {
    // Orthogonality: sum of |x_i^T x_j - delta_ij|
    for (int i = 0; i < r; ++i) {
      for (int j = i; j < r; ++j) {
        v += std::fabs(y(i, j) - (i == j ? 1.0 : 0.0));
      }
    }
  } else {
    // Uncorrelatedness: sum of |x_i^T Sigma x_j| for i != j, plus ||x_i||^2 - 1 on diagonal
    Eigen::MatrixXd C = x.transpose() * Sigma * x;
    for (int i = 0; i < r; ++i) {
      for (int j = i; j < r; ++j) {
        v += std::fabs(i == j ? y(i, j) - 1.0 : C(i, j));
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
    int feasibilityConstraintType = 0, // 0 = orthogonality constraints, 1 = uncorrelatedness constraints
    double feasibilityTolerance = 1e-4,
    double stallingTolerance = 1e-8,
    int maxIterTPM = 20,
    int timeLimitTPM = 20,
    int restartsAfterFirstIter = 10) // random restart budget for outer iterations >= 2
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
  double scalingLambda = 1; // For memory: Scaling factor used to ensure the matrix passed to TPW is PSD

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

      // Build deflation functor — avoids materialising sigma_current (saves O(n^2) allocation + symmetrisation).
      // applyM(beta) = Sigma*beta - W * diag(d) * W^T * beta, where s-th column of W is w_s and d_s = lambda * weight_s.
      int nOther = r - 1;
      Eigen::MatrixXd W(n, nOther);
      Eigen::VectorXd d(nOther);
      {
        int col = 0;
        for (int s = 0; s < r; s++)
        {
          if (s != t)
          {
            if (feasibilityConstraintType == 0) {
              W.col(col) = x_current.col(s);
            } else {
              W.col(col) = Sigma * x_current.col(s);
            }
            d(col) = theLambda * weights[s];
            col++;
          }
        }
      }

      auto applyM = [&Sigma, &W, &d, &scalingLambda, nOther](const Eigen::VectorXd& beta) -> Eigen::VectorXd {
        Eigen::VectorXd y = Sigma * beta;
        if (nOther > 0) {
          Eigen::VectorXd c = W.transpose() * beta;
          y.noalias() -= W * (d.asDiagonal() * c);
          y.noalias() += scalingLambda * d.sum() * beta;
        }
        return y;
      };

      // Warm start: power iterations on iteration 1 (no eigensolver), previous solution thereafter.
      Eigen::VectorXd beta0;
      if (theIter == 1) {
        beta0 =  Eigen::VectorXd::Ones(n);
        beta0.normalize();
      } else {
        beta0 = x_current.col(t);
      }

      singlePCHeuristic(ks[t], applyM, n, beta0, lambda_partial, x_output,
            (theIter == 1 ? maxIterTPM : restartsAfterFirstIter), timeLimitTPM);

      x_current.col(t) = x_output;

      if (theIter == 1) // Initialize the weights on each PC at the first iteration
      {
        if (t==1 && feasibilityConstraintType == 1) { // For the uncorrelatedness constraints, we need to ensure the matrix passed to TPW is PSD, which requires an appropriate scaling of the penalty parameter lambda (scalingLambda). We set this scaling to lambda_partial^2, where lambda_partial is the variance explained by the first PC, which is an upper bound on the largest sparse eigenvalue of Sigma and thus ensures the matrix passed to TPW is PSD.
          scalingLambda = lambda_partial*lambda_partial; // Initial upper bound on the largest sparse eigenvalue of Sigma for the PSD shift
        }
        if (feasibilityConstraintType == 0) {
              weights[t] = lambda_partial;
        } else {
          weights[t] = 1.0;
        }        
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
    
    bool stopCriterion = (theIter == maxIter) || (std::fabs(ofv_prev - ofv_overall) < stallingTolerance && violation < feasibilityTolerance);
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

    if (violation < feasibilityTolerance || (theIter == maxIter && ofv_best < 0)) //If current solution is feasible (within tolerance) or if we reached the last iteration and no feasible solution was found (ofv_best still <0)
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

    if (std::fabs(ofv_prev - ofv_overall) < stallingTolerance && violation < feasibilityTolerance) //If the algorithm is stalling (in terms of objective value) and the current solution is feasible (within tolerance)
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
  if (violation_best > feasibilityTolerance) // Warning if the best solution found is not feasible (within tolerance)
  {
    warning("Warning: Algorithm terminated without finding a feasible solution (within %i tolerance). Best solution found is %i feasible", feasibilityTolerance, violation_best);
  }
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
