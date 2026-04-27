#include "newton.h"

DVector getZero(DPoincareMap &pm, DVector &x0, double precision)
{
    DVector xn(x0);
    int N = xn.dimension() - 1;

    DMatrix Dphi(N + 1, N + 1);
    DMatrix Id = vectalg::Matrix<double, 0, 0>::Identity(N + 1);

    DVector P = pm(xn, Dphi); // but we actually care about P-Id
    DMatrix DF = pm.computeDP(P, Dphi) - Id;
    DMatrix invDF = matrixAlgorithms::gaussInverseMatrix(DF);

    int i = 0;
    capd::vectalg::MaxNorm<DVector, DMatrix> maxNorm;
    while (maxNorm(P - xn) > precision)
    {
        i++;
        xn = xn - invDF * (P - xn);
        P = pm(xn, Dphi);
        DF = pm.computeDP(P, Dphi) - Id;
        invDF = matrixAlgorithms::gaussInverseMatrix(DF);
    }
    // cout << "\nNewton's method (" << i << " steps), found:\n"
    //      << xn << "\n\n";
    return xn;
}

DVector getCandidate(DMap &map, DPoincareMap &pm, int N, double aStart, double aEnd)
{
    /* searching for stationary point for lower papameter a */
    map.setParameters({aStart});
    double precision = 1e-15;

    DVector start(N + 1);
    start[0] = 0;
    start[N] = 1.44;
    for (int i = 1; i < N; i++)
        start[i] = 1 / 4.;

    DVector candidate = getZero(pm, start, precision);

    double a = aStart;
    while (a < aEnd)
    {   
        a += .002;
        map.setParameters({a});
        candidate = getZero(pm, candidate, precision);
    }

    map.setParameters({aEnd});
    candidate = getZero(pm, candidate, precision);

    return candidate;
}

void plotRectangles(const IVector &rect1, const IVector &rect2)
{
    int N = rect1.dimension() - 1;
    FILE *gp = popen("gnuplot -persistent", "w");

    const char *colors[] = {
        "#F781BF",
        "#984EA3",
        "#377EB8",
        "#4DAF4A",
        "#685117",
    };
    int nColors = sizeof(colors) / sizeof(colors[0]);

    int nPlots = N;
    int cols = 3;
    int rows = (nPlots + cols - 1) / cols;

    fprintf(gp, "set terminal qt size 1600,900 font 'Sans,8'\n");
    fprintf(gp, "set multiplot layout 2,3 spacing 0.06,0.1\n");
    fprintf(gp, "set lmargin 10\n");
    fprintf(gp, "set rmargin 6\n");
    fprintf(gp, "set tmargin 3\n");
    fprintf(gp, "set bmargin 6\n");

    fprintf(gp, "unset xtics\n"); 
    fprintf(gp, "set format y '%%.1e'\n");

    for (int k = 2; k < N + 1; k++)
    {
        const char *col2 = colors[(k - 2) % nColors];

        double x1 = rect1[1].leftBound(), x2 = rect1[1].rightBound();
        double y1 = rect1[k].leftBound(), y2 = rect1[k].rightBound();
        double rx1 = rect2[1].leftBound(), rx2 = rect2[1].rightBound();
        double ry1 = rect2[k].leftBound(), ry2 = rect2[k].rightBound();

        double xmin = std::min(x1, rx1);
        double xmax = std::max(x2, rx2);
        double ymin = std::min(y1, ry1);
        double ymax = std::max(y2, ry2);

        double xpad = (xmax - xmin) * 0.05 + 1e-30;
        double ypad = (ymax - ymin) * 0.05 + 1e-30;

        fprintf(gp, "set lmargin 14\n");
        
        fprintf(gp, "set xlabel 'coord 1'\n");
        fprintf(gp, "set ylabel 'coord %d' offset 1,0\n", k); 
        fprintf(gp, "set xrange [%e:%e]\n", xmin - xpad, xmax + xpad);
        fprintf(gp, "set yrange [%e:%e]\n", ymin - ypad, ymax + ypad);
        fprintf(gp, "set key top right\n");
        fprintf(gp, "set format x ''\n"); 
        fprintf(gp, "set format y '%%.2e'\n");

        fprintf(gp, "plot '-' with filledcurves closed lc rgb '#E41A1C' fs solid 0.4 title 'rect1', "
                    "'-' with filledcurves closed lc rgb '%s' fs solid 0.4 title 'rect2'\n",
                col2);

        fprintf(gp, "%e %e\n%e %e\n%e %e\n%e %e\n%e %e\ne\n",
                x1, y1, x2, y1, x2, y2, x1, y2, x1, y1);
        fprintf(gp, "%e %e\n%e %e\n%e %e\n%e %e\n%e %e\ne\n",
                rx1, ry1, rx2, ry1, rx2, ry2, rx1, ry2, rx1, ry1);
    }

    fprintf(gp, "unset multiplot\n");
    fflush(gp);
    pclose(gp);
}