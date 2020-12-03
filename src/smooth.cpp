#include "smooth.h"
#include "cotmatrix.h"
#include "massmatrix.h"
#include <igl/edge_lengths.h>
#include <Eigen/SparseCholesky>
void smooth(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & G,
    double lambda,
    Eigen::MatrixXd & U)
{
  // Replace with your code
  Eigen::MatrixXd l;
  igl::edge_lengths(V, F, l);

  Eigen::SparseMatrix<double> L;
  cotmatrix(l, F, L);

  Eigen::DiagonalMatrix<double,Eigen::Dynamic> M;
  massmatrix(l, F, M);

  // To solve Mut = (M-lambda*L)u(t+1), let b=Mut, let A = M-lambda*L

  // A = M-lambda*L
  Eigen::SparseMatrix<double> A = - lambda * L;
  for (int i = 0; i < M.rows(); i++) {
    A.coeffRef(i, i) += M.diagonal()(i);
  }
  
  // b = M*ut in our case M*G
  Eigen::MatrixXd b = M*G;
  // solve Ax=b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  solver.compute(A);
  
  U = solver.solve(b);

}
