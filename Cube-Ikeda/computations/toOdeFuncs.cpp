#include "toOdeFuncs.h"

/* functions needed for calculating approxMatrix */
DVector calculate_nodes(double tau, // time delay parameter
                        int N       // dimention
) {
  DVector s(N + 1);
  const double step = M_PI / double(N);
  for (int i = 0; i <= N; i++)
    s[i] = tau / 2. * (cos(double(i) * step) - 1.);
  return s;
}

DVector calculate_c(const DVector &s, // chebyshev nodes
                    int N             // dimention
) {
  DVector c(N + 1);
  for (int i = 0; i <= N; i++) {
    c[i] = 1.;
    for (int j = 0; j <= N; j++) {
      if (j != i)
        c[i] *= (s[i] - s[j]);
    }
  }
  return c;
}

/* computing approxMatrix */

DMatrix compute_approxMatrix(
    double tau,       // time delay parameter
    int N,            // dimention
    string filename // whether or not to save the matrix in .txt
) {
  DMatrix M(N + 1, N + 1);

  DVector s = calculate_nodes(tau, N);
  DVector c = calculate_c(s, N);

  for (int i = 0; i <= N; i++) {
    double tmp = 0.;
    for (int j = 0; j <= N; j++) {
      if (i != j) {
        double val = c[i] / (c[j] * (s[i] - s[j]));
        M[i][j] = val;
        tmp -= val;
      }
    }
    M[i][i] = tmp;
  }
  if (filename.compare("")) {
    std::ofstream out(filename);
    for (int i = 0; i <= N; i++)
      out << M[i] << "\n";
    out.close();
  }
  return M;
}
