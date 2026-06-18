#pragma once
// Covariance-operator abstraction for msPCA.
//
// The sparse-PCA core (see msPCA_R_CPP.cpp) accesses the empirical covariance
// matrix Sigma only through a small, fixed set of operations: the matvec
// Sigma * v, the trace tr(Sigma), the quadratic form tr(x^T Sigma x), the small
// Gram matrix x^T Sigma x, and a leading-eigenvector seed. Abstracting these
// behind CovOp lets a single core serve two input modes:
//
//   * DenseOp -- wraps an explicit p x p covariance/correlation matrix. 
//   * GramOp  -- wraps the (centered/scaled) n x p data matrix X and applies
//                Sigma = s * X^T X implicitly as s * X^T(X v), with s = 1/(n-1)
//                or 1/n. It never materialises the p x p matrix, turning each
//                matvec into an O(np) operation instead of O(p^2). 

#include <RcppEigen.h>

struct CovOp {
  // Ambient dimension p (length of each loading vector).
  virtual int dim() const = 0;

  // Sigma * v.
  virtual Eigen::VectorXd apply(const Eigen::VectorXd& v) const = 0;

  // tr(Sigma) -- total variance.
  virtual double trace() const = 0;

  // tr(x^T Sigma x) -- the objective (sum of per-PC variances).
  virtual double quadForm(const Eigen::MatrixXd& x) const = 0;

  // x^T Sigma x -- small r x r Gram matrix used for the uncorrelatedness violation.
  virtual Eigen::MatrixXd gram(const Eigen::MatrixXd& x) const = 0;

  // Leading-eigenvector seed for the Truncated Power Method.
  virtual Eigen::VectorXd seed() const = 0;

  virtual ~CovOp() = default;
};

// Dense covariance/correlation matrix (p x p). Backward-compatible behaviour.
struct DenseOp : public CovOp {
  const Eigen::MatrixXd& S;

  explicit DenseOp(const Eigen::MatrixXd& S_) : S(S_) {}

  int dim() const override { return static_cast<int>(S.rows()); }

  Eigen::VectorXd apply(const Eigen::VectorXd& v) const override { return S * v; }

  double trace() const override { return S.trace(); }

  double quadForm(const Eigen::MatrixXd& x) const override {
    return (x.transpose() * S * x).trace();
  }

  Eigen::MatrixXd gram(const Eigen::MatrixXd& x) const override {
    return x.transpose() * S * x;
  }

  Eigen::VectorXd seed() const override {
    // Exact leading eigenvector via a self-adjoint eigensolver (original behaviour).
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(S);
    int index;
    solver.eigenvalues().maxCoeff(&index);
    return solver.eigenvectors().col(index);
  }
};

// Raw-data operator: Sigma = invDivisor * X^T X, applied implicitly.
// X is the (already centered, and optionally scaled) n x p data matrix.
struct GramOp : public CovOp {
  const Eigen::MatrixXd& X; // nObs x p
  double s;                 // 1/(nObs - 1) or 1/nObs

  GramOp(const Eigen::MatrixXd& X_, double s_) : X(X_), s(s_) {}

  int dim() const override { return static_cast<int>(X.cols()); }

  // Sigma * v = s * X^T (X v): two matrix-vector passes over X, cost O(np).
  //
  // CRITICAL: the parenthesisation X^T (X v) must be preserved. Writing it as
  // `X.transpose() * X * v` would make Eigen evaluate left-to-right and first
  // build the p x p Gram matrix X^T X -- O(np^2) time and O(p^2) memory -- which
  // defeats the entire purpose of this operator. Keep the inner product explicit.
  Eigen::VectorXd apply(const Eigen::VectorXd& v) const override {
    const Eigen::VectorXd Xv = X * v;          // n-vector, O(np)
    Eigen::VectorXd out(X.cols());             // p-vector
    out.noalias() = X.transpose() * Xv;        // O(np), no aliasing temporary
    out *= s;
    return out;
  }

  // tr(Sigma) = s * ||X||_F^2 = s * sum_j ||X_j||^2. Single O(np) pass.
  double trace() const override { return X.squaredNorm() * s; }

  // tr(x^T Sigma x) = s * ||X x||_F^2. Cost O(npr); never forms the p x p matrix.
  double quadForm(const Eigen::MatrixXd& x) const override {
    const Eigen::MatrixXd Xx = X * x;          // n x r, O(npr)
    return Xx.squaredNorm() * s;
  }

  // x^T Sigma x = s * (X x)^T (X x), an r x r matrix.
  // Cost O(npr) for X x plus O(nr^2) for the small Gram; never forms the p x p matrix.
  Eigen::MatrixXd gram(const Eigen::MatrixXd& x) const override {
    const Eigen::MatrixXd Xx = X * x;          // n x r, O(npr)
    Eigen::MatrixXd G(x.cols(), x.cols());     // r x r
    G.noalias() = Xx.transpose() * Xx;
    G *= s;
    return G;
  }

  // Leading-eigenvector seed via matvec power iteration (no p x p, no O(p^3) solver).
  Eigen::VectorXd seed() const override {
    const int p = dim();
    Eigen::VectorXd b = Eigen::VectorXd::Ones(p);
    b.normalize();
    for (int i = 0; i < 100; ++i) {
      Eigen::VectorXd nb = apply(b);
      double nm = nb.norm();
      if (nm < 1e-12) break;
      Eigen::VectorXd next = nb / nm;
      if ((next - b).squaredNorm() < 1e-16) { b = next; break; }
      b = next;
    }
    return b;
  }
};
