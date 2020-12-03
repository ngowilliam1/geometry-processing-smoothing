#include "massmatrix.h"
#include <igl/doublearea.h>
void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
  // Add your code here
  M.resize(F.maxCoeff() + 1);
  Eigen::VectorXd doubleArea;
  igl::doublearea(l, doubleArea);

  for (int i=0 ; i < F.rows(); i++){
    // divide by 6 since we have double area and divide by 3
    M.diagonal()(F(i, 0)) += (doubleArea(i)/6.0);
	  M.diagonal()(F(i, 2)) += (doubleArea(i)/6.0);
	  M.diagonal()(F(i, 2)) += (doubleArea(i)/6.0);
  }
  
}

