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
    int dimIn = N + 1, dimOut = N + 1, noParams = 1, highestDerivative = 1;
    double a = 1.53;
    int taylorOrder = 20;
    double returnTime = 0.;
    double chaoticA = 1.542;
    double tau = 1;
    string filename = "";

    approxMatrix = compute_approxMatrix(tau, N, filename);

    DMap CubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
                    highestDerivative);
    CubicIkeda.setParameters({a});

    DOdeSolver solver(CubicIkeda, taylorOrder);

    /* defining Poincare map */
    DCoordinateSection section(N + 1, 0, 0);
    DPoincareMap pm(solver, section, poincare::MinusPlus);

    cout << "here\n";
    DVector periodicPoint = getCandidate(CubicIkeda, pm, N, a, chaoticA);

    // DTimeMap timeMap(solver);
    // getSolutionCurve(timeMap, periodicPoint, 100., 0., filenameSolutionCurve); // looks good
    
    DMatrix Dphi(N + 1, N + 1);
    DVector P = pm(periodicPoint, Dphi);
    DMatrix DP = pm.computeDP(P, Dphi);
    DVector eigenRealPart(N + 1), eigenImPart(N + 1);
    DMatrix vectorRealPart(N + 1, N + 1), vectorImPart(N + 1, N + 1);
    alglib::computeEigenvaluesAndEigenvectors(DP, eigenRealPart, eigenImPart,
                                              vectorRealPart, vectorImPart);

    cout << "eigenvectors:\n" << vectorRealPart << endl;

    /* eigenvectors corresponding to given eigenvalue*/
    std::map<double, DVector> eigenVectors;
    for (int i = 0; i < N + 1; i++)
    {
        eigenVectors[eigenRealPart[i]] = vectorRealPart.column(i);
    }

    std::sort(eigenRealPart.begin(), eigenRealPart.end(), [](double i, double j)
              { return abs(i) > abs(j); });

    cout << endl
         << "real part of eigenvalues:\n"
         << eigenRealPart << endl
         << "imaginary part of eigenvalues:\n"
         << eigenImPart << endl
         << endl;
    /* so we see that we only have real eigenvalues :) */


    /* we define r specifically for N = 6 */
    interval data[] = {interval(0), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1)};
    IVector r(N + 1, data);
    r *= 1e-9;

    IMap ICubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
                     highestDerivative);
    ICubicIkeda.setParameters({chaoticA});

    IVector x0(periodicPoint);

    IMatrix B(N + 1, N + 1);
    auto fx0 = ICubicIkeda(x0);
    fx0.normalize();
    // cout << fx0 << endl;
    B[0][0] = fx0[0];
    for (int i = 1; i < N + 1; i++)
    {
        B[i][0] = fx0[i];
        auto vec = eigenVectors[eigenRealPart[i - 1]];
        vec.normalize();
        for (int j = 0; j < N + 1; j++)
        {
            B[j][i] = vec[j];
        }
    }
    
    IMatrix A = matrixAlgorithms::gaussInverseMatrix(B);
    
    C1Rect2Set s(x0, B, r);
    IOdeSolver ISolver(ICubicIkeda, taylorOrder);
    
    /* defining Poincare map */
    ICoordinateSection ISection(N + 1, 0, 0);
    IPoincareMap IPm(ISolver, ISection, poincare::MinusPlus);
    
    interval IReturnTime = 0;
    x0[0] = 0.;
    IVector P0 = IPm(s, x0, A, IReturnTime);

    // cout << r << endl;
    // cout << B * r << endl;
    // cout << x0 << endl;
    // cout << "\n" << B << endl << endl;

    cout <<  "r: " << r << endl
         << endl;
    cout << "P0: " << P0 << endl
         << endl;

    plotRectangles(r, P0);
}
