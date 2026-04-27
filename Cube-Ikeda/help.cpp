#include "computations/integrateFuncs.h"
#include "computations/newton.h"
#include "computations/toOdeFuncs.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
using namespace capd;
using namespace std;

// the matrix needed for pseudospectral approximation
DMatrix approxMatrix;

/*
Pseudospectral approximation of Cubic Ikeda
f(x(t), x(t-tau)) = a(x(t-tau) - x(t-tau)^3)
*/

void approx_CubicIkeda(capd ::autodiff ::Node /*t*/, // unused time variable
                       capd ::autodiff ::Node in[],
                       int dimIn, // input variables x1 ,... , xn
                       capd ::autodiff ::Node out[],
                       int dimOut, // output : function values
                       capd ::autodiff ::Node params[],
                       int noParam // parameters
)
{
    int N = dimIn - 1;

    // out[0] = f(x(t), x(t-tau))
    //        = f(in[0], in[N])
    //        = a(in[N]-in[N])^3
    capd ::autodiff ::Node a = params[0];
    capd ::autodiff ::Node xDelayed = in[N];
    out[0] = a * xDelayed * (1 - xDelayed * xDelayed);

    for (int i = 1; i <= N; i++)
    {
        capd ::autodiff ::Node result(0.);
        for (int j = 0; j <= N; j++)
            result += approxMatrix[i][j] * in[j];
        out[i] = result;
    }
}

// MAIN FUNCTION // ---------------------------------------------------------------------

int main()
{
    const int N = 6;
    int dimIn = N + 1, dimOut = N + 1, noParams = 1, highestDerivative = 2;
    double a = 1.53;
    int taylorOrder = 20;
    double returnTime = 0.;
    double chaoticA = 1.546;
    double tau = 1;
    string filename = "";

    approxMatrix = compute_approxMatrix(tau, N, filename);

    DMap CubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
                    highestDerivative);
    CubicIkeda.setParameters({a});
    DOdeSolver solver(CubicIkeda, taylorOrder);
    DCoordinateSection section(N + 1, 0, 0);
    DPoincareMap pm(solver, section, poincare::MinusPlus);

    DVector candidate = getCandidate(CubicIkeda, pm, N, a, chaoticA);
    CubicIkeda.setParameters({chaoticA});

    DMatrix Dphi(N + 1, N + 1);
    DVector P = pm(candidate, Dphi);
    DMatrix DP = pm.computeDP(P, Dphi);
    cout << "P(candidate) - candidate:\n\t" << P - candidate;

    auto v = CubicIkeda(candidate);
    cout << "\nf(candidate):\n\t" << v << endl;
    v.normalize();
    cout << "normalized f(candidate):\n\t" << v << endl;

    cout << "DP:\n"
         << DP << endl;

    cout << "DP * f(candidate):\n\t" << DP * v << endl;

    DVector eigenRealPart(N + 1), eigenImPart(N + 1);
    DMatrix vectorRealPart(N + 1, N + 1), vectorImPart(N + 1, N + 1);
    alglib::computeEigenvaluesAndEigenvectors(DP, eigenRealPart, eigenImPart,
                                              vectorRealPart, vectorImPart);

    cout << "eigenvalues:\n\t" << eigenRealPart << endl;
    cout << "eigenvectors:\n"
         << vectorRealPart << endl;
}
