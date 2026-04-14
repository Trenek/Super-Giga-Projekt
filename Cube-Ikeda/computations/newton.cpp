#include "newton.h"

DVector getZero(DPoincareMap& pm, DVector& x0, int k){
    DVector xn(x0);
    int N = xn.dimension() - 1;
    DVector ones(N+1);
    for(int i = 0; i<=N; i++)
        ones[i] = 1;
    double rt = 0.;
    DMatrix Dphi(N+1, N+1);
    DMatrix Id = vectalg::Matrix<double, 0, 0>::Identity(N+1);
    DVector P = pm(xn,Dphi);  // but we actually care about P-Id
    DMatrix DF = pm.computeDP(P,Dphi) - Id;
    DMatrix invDF = matrixAlgorithms::gaussInverseMatrix(DF);

    for(int i = 1; i < k; i++){
        xn = xn - invDF * (P-xn);
        P = pm(xn,Dphi);
        DF = pm.computeDP(P,Dphi) - Id;
        invDF = matrixAlgorithms::gaussInverseMatrix(DF);
    }
    return xn;
}