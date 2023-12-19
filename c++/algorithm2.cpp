// to compile
// g++ -shared -o algorithm2.so algorithm2.h -I/path/to/eigen3/includes
// EXAMPLE : g++ -shared -o algorithm2.so algorithm2.h -I/usr/include/eigen3

#include <cstdio>
#include <ctime>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

struct Problem
{
    VectorXd data;
    MatrixXd Sigma;
};

double evaluate(VectorXd &solution, Problem &prob)
{
    return (solution.transpose() * prob.Sigma)[0];
}

vector<int> selectperm2(VectorXd &x, int k)
{
    vector<int> z(x.size());
    for (int i = 0; i < z.size(); i++)
    {
        z[i] = i;
    }
    partial_sort(
        x.begin(), x.begin() + k, x.end(), [&](int i, int j) -> bool
        { return x[i] > x[j]; });
    z.resize(k);
    return z;
}

VectorXd Hk(VectorXd origlist, int sparsity, vector<int> &support)
{
    VectorXd list = origlist;
    VectorXd kparse = VectorXd::Zero(list.size());
    int nbIndicesToKeep = 0;
    for (int s : support)
    {
        nbIndicesToKeep += s == 1;
    }
    double dummyValue = list[0] - 1;
    for (double x : list)
    {
        dummyValue = min(dummyValue, x - 1);
    }
    for (int i = 0; i < list.size(); i++)
    {
        if (support[i] > -1)
        {
            list[i] = dummyValue;
        }
    }
    vector<int> newIndices = selectperm2(list, sparsity - nbIndicesToKeep);
    int nbIndices = 0;
    for (int i = 0; i < kparse.size(); i++)
    {
        if (support[i] == 1 || newIndices[i] == 1)
        {
            nbIndices++;
            kparse[i] = origlist[i];
        }
    }
    kparse = kparse / nbIndices;
    return kparse;
}

VectorXd eigSubset(Problem &prob, vector<int> &support, int k, VectorXd &beta0)
{
    VectorXd beta = Hk(beta0, k, support);
    for (int i = 1; i < 100; i++)
    {
        beta = Hk(prob.Sigma * beta, k, support);
    }
    return beta;
}

pair<double, VectorXd> subset(Problem prob, int k, int timeLimit = 720, vector<int> support = {0}, int countdown = 100)
{
    int n = prob.Sigma.rows();
    Eigen::SelfAdjointEigenSolver<MatrixXd> solver(prob.Sigma);
    int index;
    solver.eigenvalues().maxCoeff(&index);
    VectorXd beta0 = solver.eigenvectors().col(index);
    if (support.size() == 1)
    {
        support.resize(n);
        fill(support.begin(), support.end(), -1);
    }
    VectorXd bestBeta = eigSubset(prob, support, k, beta0);
    double bestObj = evaluate(bestBeta, prob);
    time_t start = time(0);
    int margin = 1;
    while (countdown > 0 && difftime(time(0), start) < timeLimit)
    {
        VectorXd beta = VectorXd::Random(n).cwiseAbs(); // uniform on (-1, 1)
        beta = beta / beta0.norm();
        beta = eigSubset(prob, support, k, beta);
        double obj = evaluate(beta, prob);
        if (obj > bestObj)
        {
            bestObj = obj;
            bestBeta = bestBeta;
            countdown = 100;
        }
        countdown--;
    }
    return {bestObj, bestBeta};
}

double fnviolation(MatrixXd &x)
{
    double v = 0;
    for (int j = 0; j < x.cols(); j++)
    {
        v += abs(x.col(j).squaredNorm() - 1);
    }
    return v;
}

tuple<double, double, double, MatrixXd> cpp_findmultPCs_deflation(
    MatrixXd &Sigma,
    int r,
    int *ks, // size r
    int numIters = 200,
    bool verbose = true,
    double violation_tolerance = 1e-4)
{
    int n = Sigma.rows();
    double ofv_best = -1e10;
    double violation_best = n;
    MatrixXd x_best = MatrixXd::Zero(n, r);
    MatrixXd x_current = MatrixXd::Zero(n, r);
    double violation_init = max(1e-7, fnviolation(x_current));
    double ofv_prev = 0;
    double ofv_overall = 0;

    VectorXd weights = VectorXd::Zero(r);
    double theLambda = 0;
    double stepSize = 0;
    int slowPeriod = ceil(0.15 * numIters);
    int fastPeriod = ceil(0.75 * numIters);

    if (verbose)
    {
        printf("---- Iterative deflation algorithm for sparse PCA with multiple PCs ---\n");
        printf("Dimension: %d\n", n);
        printf("Number of PCs: %d\n", r);
        printf("Sparsity pattern:");
        for (int t = 0; t < r; t++)
        {
            printf(" %d", ks[t]);
        }
        printf("\n\n %10s | %20s | %25s | %10s\n", "Iteration", "Objective value", "Orthogonality Violation", "Time");
    }

    time_t start_time = time(0);
    for (int theIter = 1; theIter < numIters; theIter++)
    {
        theLambda += stepSize;
        for (int t = 0; t < r; t++)
        {
            MatrixXd sigma_current = Sigma;
            for (int s = 0; s < r; s++)
            {
                if (s != t)
                {
                    sigma_current -= theLambda * weights[s] * x_current.col(s) * x_current.col(s).transpose();
                }
            }
            sigma_current = (sigma_current + sigma_current.transpose()) / 2;
            Eigen::SelfAdjointEigenSolver<MatrixXd> solver(sigma_current);
            double lambda0 = -solver.eigenvalues().minCoeff() + 1e-4;
            for (int i = 0; i < sigma_current.rows(); i++)
            {
                sigma_current(i, i) += lambda0;
            }
            solver = Eigen::SelfAdjointEigenSolver<MatrixXd>(sigma_current);
            auto [lambda_partial, x_output] = subset({solver.operatorSqrt(), sigma_current}, ks[t], 20);
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

        double violation = max(1e-7, fnviolation(x_current));
        stepSize = (theIter < fastPeriod ? 0.01 : 0.05) * (theIter < slowPeriod ? violation : ofv_overall / violation);
        if (verbose)
        {
            if (numIters <= 25 || theIter % 10 == 0)
            {
                double timespan = difftime(time(0), start_time);
                printf(" %10d | %20.3f | %25.2e | %10.3f \n", theIter, ofv_overall / Sigma.trace(), violation, timespan);
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

        if (abs(ofv_prev - ofv_overall) < 1e-8 && violation < violation_tolerance)
        {
            if (ofv_best < 0)
                x_best = x_current;
            ofv_best = (x_current.transpose() * Sigma * x_current).trace();
        }
        break;
    }
    double runtime = difftime(time(0), start_time);
    violation_best = fnviolation(x_best);
    ofv_best = (x_best.transpose() * Sigma * x_best).trace();
    return {ofv_best, violation_best, runtime, x_best};
}

extern "C"
void findmultPCs_deflation(
    double *ofv_best,       // size ()
    double *violation_best, // size ()
    double *runtime,        // size ()
    double *x_best,         // size (n, r)
    int n,
    int r,
    double *Sigma, // size (n)
    int *ks,   // size (r)
    int numIters = 200,
    int verbose = true,
    double violation_tolerance = 1e-4)
{
    MatrixXd S(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            S(i, j) = Sigma[i * n + j];
        }
    }
    auto ans = cpp_findmultPCs_deflation(S, r, ks, numIters, verbose, violation_tolerance);
    *ofv_best = get<0>(ans);
    *violation_best = get<1>(ans);
    *runtime = get<2>(ans);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < r; j++)
        {
            x_best[i * r + j] = get<3>(ans)(i, j);
        }
    }
}
