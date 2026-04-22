#include <capd/capdlib.h>

#include "mackeyGlass.hpp"
#include "draw.hpp"

#define N 5

capd::LDVector Newton(capd::LDVector u, capd::LDPoincareMap map) {
    capd::vectalg::MaxNorm<capd::LDVector, capd::LDMatrix> maxNorm;

    capd::LDMatrix Dphi(N + 1, N + 1);
    capd::LDVector P = map(u, Dphi);

    size_t i = 0;

    while (maxNorm(P - u) >= 10e-18) {
        std::print("Iteracja {}: Blad {}\n", i, maxNorm(P - u));
        u -= capd::matrixAlgorithms::gaussInverseMatrix(
            map.computeDP(P, Dphi) - capd::LDMatrix::Identity(N + 1)
        ) * (P - u);
        P = map(u, Dphi);
    }

    return u;
}

void getEigenValues(double n, capd::LDVector u0, capd::LDPoincareMap map) {
    capd::LDMatrix dPhi(N + 1, N + 1);

    capd::LDVector u1 = map(u0, dPhi);
    capd::LDMatrix dP = map.computeDP(u0, dPhi);

    capd::LDVector eigRe(N + 1);
    capd::LDVector eigIm(N + 1);
    capd::LDMatrix vectRe(N + 1, N + 1);
    capd::LDMatrix vectIm(N + 1, N + 1);

    capd::alglib::computeEigenvaluesAndEigenvectors(dP, eigRe, eigIm, vectRe, vectIm);

    std::print("n = {}, Wartosci wlasne {}\n", n, capd::alglib::eigenvaluesToString(eigRe, eigIm));
    std::print("n = {}, Wektory wlasne {}\n", n, capd::alglib::eigenvectorsToString(vectRe, vectIm));
}

int main() {
    class gnuPlotManager manager{{
        {
            .name = "Newton",
            .file = "new.dat",

            .xName = "n",
            .yName = "xN"
        },
    }};

    constexpr uint32_t order = 20;
    capd::LDMap f{mackeyGlass<N>, N + 1, N + 1, 1};
    capd::LDOdeSolver solver{f, order}; {
        solver.setStep(0.01);
    }
    capd::LDCoordinateSection section{N + 1, 0, 0.6};
    capd::LDPoincareMap map{solver, section, capd::poincare::MinusPlus};

    capd::LDVector cyclic;
    capd::LDVector u(N + 1); {
        for (auto &e : u) {
            e = 1.1; 
        }
        u[0] = 0.6; 
    }
    for (double n = 8.7; n < 9; n += 0.01) {
        f.setParameter(0, n);

        cyclic = Newton(u, map);
        getEigenValues(n, cyclic, map);

        std::print("Orbita: {}\n", cyclic);
        manager.print(0, "{} {}\n", n, cyclic[N]);
        manager.fflush();

        manager.initGNUPlot();
    }

    std::cout << map(u) - u << std::endl;

    return 0;
}

// metoda newtona - orbita
// postawić zbiór
// policzyć własności własne
// Relacja Nakrywająca
