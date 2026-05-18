#include <capd/capdlib.h>

#include "mackeyGlass.hpp"
#include "draw.hpp"

#define N 5

static void setGNUPlot(int id, struct thing &drawer) {
    fprintf(drawer.gnuplot, "set term qt %d size 800,600\n", id);
    fprintf(drawer.gnuplot, "set title '%s'\n", drawer.name);
    fprintf(drawer.gnuplot, "set xlabel '%s'\n", drawer.xName);
    fprintf(drawer.gnuplot, "set ylabel '%s'\n", drawer.yName);
    fprintf(drawer.gnuplot, "set pointsize 1.5\n");
    fprintf(drawer.gnuplot, "set grid\n");

    fprintf(drawer.gnuplot, "pause 1\n");

    fprintf(drawer.gnuplot, "plot "
            "\"%s\" index 0 with lines lw 2 lc rgb 'orange' title 'X', "
            "\"%s\" index 1 with lines lw 2 lc rgb 'cyan' title 'f^-1(X) Min', "
            "\"%s\" index 2 with lines lw 2 lc rgb 'red' title 'f^-1(X) Max'\n",
    drawer.file, drawer.file, drawer.file);

    fprintf(drawer.gnuplot, "bind 'q' 'exit'\n");
    fprintf(drawer.gnuplot, "while (1) {\n");
    fprintf(drawer.gnuplot, "    pause 1\n");
    fprintf(drawer.gnuplot, "    replot\n");
    fprintf(drawer.gnuplot, "}\n");
    fflush(drawer.gnuplot);
}

static capd::LDVector Newton(capd::LDVector u, capd::LDPoincareMap map) {
    capd::vectalg::MaxNorm<capd::LDVector, capd::LDMatrix> maxNorm;

    capd::LDMatrix Dphi(N + 1, N + 1);
    capd::LDVector P = map(u, Dphi);

    size_t i = 0;

    while (maxNorm(P - u) >= 10e-18) {
        u -= capd::matrixAlgorithms::gaussInverseMatrix(
            map.computeDP(P, Dphi) - capd::LDMatrix::Identity(N + 1)
        ) * (P - u);
        P = map(u, Dphi);
    }

    return u;
}

static capd::LDMatrix calcEigenVector(capd::LDVector u0, capd::LDPoincareMap map) {
    capd::LDMatrix dPhi(N + 1, N + 1);

    capd::LDVector u1 = map(u0, dPhi);
    capd::LDMatrix dP = map.computeDP(u0, dPhi);

    capd::LDVector eigRe(N + 1);
    capd::LDVector eigIm(N + 1);
    capd::LDMatrix vectRe(N + 1, N + 1);
    capd::LDMatrix vectIm(N + 1, N + 1);

    capd::alglib::computeEigenvaluesAndEigenvectors(dP, eigRe, eigIm, vectRe, vectIm);

    return capd::matrixAlgorithms::gaussInverseMatrix(vectRe) * dP * vectRe;
}

static void drawInitialRectangle(double side, struct gnuPlotManager *manager) {
    manager->print(0, "{} {}\n", -side, -side);
    manager->print(0, "{} {}\n", side, -side);
    manager->print(0, "{} {}\n", side, side);
    manager->print(0, "{} {}\n", -side, side);
    manager->print(0, "{} {}\n", -side, -side);

    manager->print(0, "\n\n"); 
}

void drawEdgeL(double b0, double b1, double e0, double e1, capd::IMatrix A, struct gnuPlotManager *manager) {
    capd::IVector v(N + 1);
    capd::IVector img;

    for (double t = 0; t <= 1.0; t += 0.01) {
        v[0] = b0 + t * (e0 - b0);
        v[1] = b1 + t * (e1 - b1);
        
        img = A * v;
        manager->print(0, "{} {}\n", img[0].leftBound(), img[1].leftBound());
    }
};

void drawEdgeR(double b0, double b1, double e0, double e1, capd::IMatrix A, struct gnuPlotManager *manager) {
    capd::IVector v(N + 1);
    capd::IVector img;

    for (double t = 0; t <= 1.0; t += 0.01) {
        v[0] = b0 + t * (e0 - b0);
        v[1] = b1 + t * (e1 - b1);
        
        img = A * v;
        manager->print(0, "{} {}\n", img[0].rightBound(), img[1].rightBound());
    }
};

static void checkForCoveringRelationship(capd::LDVector u, capd::LDPoincareMap map, capd::IPoincareMap iMap, struct gnuPlotManager *manager, double side) {
    capd::LDVector fixed = Newton(u, map);
    capd::LDMatrix A = calcEigenVector(fixed, map);

    capd::IVector iFixed{fixed};
    capd::IMatrix iA{A};

    capd::IVector r(N + 1); {
        r[0] = capd::Interval(0);
        for (size_t i = 1; i < N + 1; i += 1) {
            r[i] = capd::Interval(-side, side);
        }
    }

    capd::C0Rect2Set x{iFixed, iA, r};
    capd::C0Rect2Set fx = iMap(x);

    capd::IMatrix B = fx.get_B();
    capd::IMatrix pB = capd::matrixAlgorithms::gaussInverseMatrix(iA) * B;

    drawEdgeL(-side, -side,  side, -side, pB, manager);
    drawEdgeL( side, -side,  side,  side, pB, manager);
    drawEdgeL( side,  side, -side,  side, pB, manager);
    drawEdgeL(-side,  side, -side, -side, pB, manager);

    manager->print(0, "\n\n");

    drawEdgeR(-side, -side,  side, -side, pB, manager);
    drawEdgeR( side, -side,  side,  side, pB, manager);
    drawEdgeR( side,  side, -side,  side, pB, manager);
    drawEdgeR(-side,  side, -side, -side, pB, manager);
}

int main() {
    class gnuPlotManager manager{{
        {
            .name = "Relacja Nakrywajaca",
            .file = "cRel.dat",

            .xName = "n",
            .yName = "s",

            .setGNUPlot = setGNUPlot,
        },
    }};

    constexpr double n = 8.7;
    constexpr uint32_t order = 20;
    capd::LDMap f{mackeyGlass<N>, N + 1, N + 1, 1};
    capd::LDOdeSolver solver{f, order}; {
        solver.setStep(0.01);
    }
    capd::LDCoordinateSection section{N + 1, 0, 0.6};
    capd::LDPoincareMap map{solver, section, capd::poincare::MinusPlus};

    capd::IMap iF{mackeyGlass<N>, N + 1, N + 1, 1};
    capd::IOdeSolver iSolver{iF, order}; {
        iSolver.setStep(0.01);
    }
    capd::ICoordinateSection iSection{N + 1, 0, 0.6};
    capd::IPoincareMap iMap{iSolver, iSection, capd::poincare::MinusPlus};

    double side = 10e-10;
    capd::LDVector u(N + 1); {
        for (auto &e : u) {
            e = 1.1; 
        }
        u[0] = 0.6; 
    }

    drawInitialRectangle(side, &manager);

    // for (double n = 8.7; n < 9; n += 0.001) {
        f.setParameter(0, n);
        iF.setParameter(0, n);
        checkForCoveringRelationship(u, map, iMap, &manager, side);

        manager.fflush();
        manager.initGNUPlot();
    // }

    return 0;
}
// drugi Puancare
// Puancare map (set, Macierz, punkt (wektor))
// C0Rect2Set <<<<<------------ lepsze
// C1Rect2Set(x0, macierz współrzędnych, promień) <------- przekomplikowane
// Tripleton??
//
// Interwały
