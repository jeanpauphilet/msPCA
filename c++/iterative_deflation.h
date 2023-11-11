#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>

using namespace std;
using namespace Eigen;

// Assuming utilities.h is the equivalent of utilities.jl
#include "utilities.h"

// Function definition
vector<VectorXd> findmultPCs_deflation(MatrixXd Sigma, int r, vector<int> ks, int numIters = 200, bool verbose = true, double violation_tolerance = 1e-4) {
    int n = Sigma.rows();

    double ofv_best = -1e10;
    int violation_best = n;
    MatrixXd x_best = MatrixXd::Zero(n, r);
    MatrixXd x_current = MatrixXd::Zero(n, r);
    double violation_init = max((x_current.transpose() * x_current - MatrixXd::Identity(r, r)).cwiseAbs().sum(), 1e-7);
    double ofv_prev = 0.0;
    double ofv_overall = 0.0;

    // Step size tuning
    VectorXd weights = VectorXd::Zero(r);
    double theLambda = 0.0;

    double stepSize = 0.0;
    int slowPeriod = ceil(0.15 * numIters);
    int fastPeriod = ceil(0.75 * numIters);

    if (verbose) {
        cout << "---- Iterative deflation algorithm for sparse PCA with multiple PCs ---" << endl;
        cout << "Dimension: " << n << endl;
        cout << "Number of PCs: " << r << endl;
        cout << "Sparsity pattern: ";
        for (int k : ks) cout << k << " ";
        cout << endl;
        printf(" %10s | %20s | %25s | %10s \n", "Iteration", "Objective value", "Orthogonality Violation", "Time");
    }

    time_t start_time = time(0);
    for (int theIter = 1; theIter <= numIters; theIter++) {
        theLambda += stepSize;

        for (int t = 0; t < r; t++) {
            MatrixXd sigma_current = Sigma;
            for (int s = 0; s < r; s++) {
                if (s != t) {
                    sigma_current -= theLambda * weights[s] * x_current.col(s) * x_current.col(s).transpose();
                }
            }
            sigma_current = (sigma_current + sigma_current.transpose()) / 2;
            SelfAdjointEigenSolver<MatrixXd> es(sigma_current);
            double lambda0 = -es.eigenvalues().minCoeff() + 1e-4;
            sigma_current += lambda0 * MatrixXd::Identity(n, n);
        }
    }
    // Return the best solution
    vector<VectorXd> result;
    for (int i = 0; i < r; i++) {
        result.push_back(x_best.col(i));
    }
    return result;
}