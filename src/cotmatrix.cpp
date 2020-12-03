#include "cotmatrix.h"
#include "cmath"
double cot(double a, double b, double c){
  // calculate area: https://en.wikipedia.org/wiki/Heron's_formula
  double s = (a+b+c) / 2.0;
  double area = sqrt(s*(s-a)*(s-b)*(s-c));
  // law of sines
  double sinA = 2.0*area/(b*c);
  // law of cosines
  double cosA = (pow(b, 2.0) + pow(c, 2.0) - pow(a, 2.0))/(2.0 * b * c);
  return cosA/sinA;
}


void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{
  // Add your code here
  typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	tripletList.reserve(F.rows()*12);

  for (int i =0; i< F.rows(); i++){
    double a = l(i,0);
    double b = l(i,1);
    double c = l(i,2);

    double cotA = cot(a,b,c);
    double cotB = cot(b,a,c);
    double cotC = cot(c,a,b);
    // index to insert in matrix L
    int A = F(i,0);
    int B = F(i,1);
    int C = F(i,2);

    tripletList.push_back(T(A, B, cotC / 2.0));
    tripletList.push_back(T(B, A, cotC / 2.0));
    tripletList.push_back(T(A, C, cotB / 2.0));
    tripletList.push_back(T(C, A, cotB / 2.0));
    tripletList.push_back(T(B, C, cotA / 2.0));
    tripletList.push_back(T(C, B, cotA / 2.0));
    
    tripletList.push_back(T(C, C, -cotB / 2.0));
    tripletList.push_back(T(C, C, -cotA / 2.0));
    tripletList.push_back(T(B, B, -cotC / 2.0));
    tripletList.push_back(T(B, B, -cotA / 2.0));
    tripletList.push_back(T(A, A, -cotC / 2.0));
    tripletList.push_back(T(A, A, -cotB / 2.0));
    
  }
  L.resize((F.maxCoeff() + 1), (F.maxCoeff() + 1));
  L.setFromTriplets(tripletList.begin(), tripletList.end());
}

