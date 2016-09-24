/*
 * Compile using : 
 * g++ -std=c++17 -Wall -Wextra -pedantic -pthread -O3 -Wshadow -Wformat
 * -Wdouble-equal -Wcast-qual -Wcast-align -Weffc++ -I/usr/include/python2.7/
 * -lboost_python -lboost_system -lpython2.7 -shared -fPIC -rdynamic hhpy.cpp
 * -o libhouseholder.so
 */


#include <boost/python.hpp>

#include <cmath>
#include <iostream>
#include <vector>

using namespace std;
using namespace boost;
using namespace boost::python;

using vvd = vector<vector<double>>;
using vd = vector<double>;

auto getTrigonalizedMatrix(vvd mat) {
  auto dim = mat.size();

  vvd ans(dim, vd(dim, 0.0));

  size_t k{ 0 };
  bool found{ false };

  while (not found && k < dim - 1) {	
		vvd x(dim, vd(1, 0.0));
		vvd y(1, vd(dim, 0.0));
		vvd tmp(dim, vd(dim, 0.0));

    double sum = 0.0;
    for (size_t i = k + 1; i < dim; ++i)
      sum += (mat[i][k] * mat[i][k]);
    sum = pow(sum, 0.5);

		double alpha = -1.0 * (mat[k + 1][k] / (abs(mat[k + 1][k]))) * sum;
		double r = (0.5 * pow(alpha, 2.0)) - (0.5 * alpha * mat[k + 1][k]);
		r = pow(r, 0.5);

		x[k + 1][0] = (1 / (2 * r)) * (mat[k + 1][k] - alpha);

		for (size_t i = k + 2; i < dim; ++i)
			x[i][0] = (mat[i][k] / (2 * r));

    for (size_t i = 0; i < dim; ++i)
      y[0][i] = x[i][0];

    for (size_t i = 0; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
				ans[i][j] = x[i][0] * y[0][j];
        if (i == j)
          ans[i][j] = 1 - 2 * ans[i][j];
        else
          ans[i][j] = -2 * ans[i][j];
      }
    }

    for (size_t i = 0; i < dim; ++i)
      for (size_t j = 0; j < dim; ++j)
        for (size_t l = 0; l < dim; ++l)
          tmp[i][j] += ans[i][l] * mat[l][j];

		swap(mat, tmp);

    for (size_t i = 0; i < dim; ++i) {
      for (size_t j = 0; j < dim; ++j) {
        tmp[i][j] = 0;
        for (size_t l = 0; l < dim; ++l)
          tmp[i][j] += mat[i][l] * ans[l][j];
      }
    }

		swap(mat, tmp);

    found = true;
    for (size_t i = 0; i < dim && found; ++i)
      for (size_t j = 0; j < dim && found; ++j)
        if (abs(mat[i][j]) > 0.0 && (j < (i - 1) || j > (i + 1)))
          found = false;
    k++;
  }
  return mat;
}

void printMatrix(vvd mat) {
  cout.setf(ios::right | ios::scientific | ios::showpos);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j)
      cout << mat[i][j] << "\t";
    cout << endl;
  }
}

int main() {
  vector<vector<double>> mat = {{	  5,   -2, -0.5,  1.5},
                                {	 -2,    5,  1.5, -0.5},
                                {-0.5,  1.5,    5,   -2},
                                { 1.5, -0.5,   -2,    5}};
  mat = getTrigonalizedMatrix(mat);
  printMatrix(mat);

  return 0;
}

auto pythonWrapper(list &pyMat) {
  vvd mat;
  auto dim = len(pyMat);
  for (int i = 0; i < dim; ++i) {
    mat.emplace_back();
    list pyMatRow = extract<list>(pyMat[i]);
    for (int j = 0; j < dim; ++j)
      mat.back().emplace_back(extract<double>(pyMatRow[j]));
  }
  mat = getTrigonalizedMatrix(mat);
  list ansMat;
  for (int i = 0; i < dim; ++i) {
    list tmpList;
    for (int j = 0; j < dim; ++j)
      tmpList.append(mat[i][j]);
    ansMat.append(tmpList);
  }
  return ansMat;
}

BOOST_PYTHON_MODULE(libhouseholder) {
  PyEval_InitThreads();
  def("GetTrigonalizedMatrix_", pythonWrapper);
}
