#include "integrateFuncs.h"

void getSolutionCurve(DTimeMap timeMap, DVector start,
                                         double time, double skip,
                                         string filename) {
  DTimeMap::SolutionCurve solution(0.);
  timeMap(time, start, solution);

  int N = start.dimension() - 1;

  if (filename.compare("")) {
    ofstream csv(filename);
    csv << "x_present,x_delayed\n";
    for (double t = skip; t <= time; t += 0.05) {
      DVector point = solution(t);
      csv << point[0] << "," << point[N] << "\n";
    }
    csv.close();
  }
}

void getPoincareValues(DPoincareMap &pm, DVector &x, const string &filename) {
  ofstream csv(filename);
  csv << x << "\n";
  DVector returnPoint = pm(x);
  csv << returnPoint << "\n";
  for (int i = 0; i < 10; i++) {
    returnPoint = pm(returnPoint);
    csv << returnPoint << "\n";
  }
  csv.close();
}

// gnuplot code written by Claude
void plotBifurcationDiagram(DMap &CubicIkeda, DPoincareMap &pm, DVector &x,
                            double aStart, double aEnd, double aFrequency,
                            int noSteps, string filename1) {
  FILE *gp = popen("gnuplot", "w"); // no -persistent needed for file output
  string filename = "images/chaos.png";
  fprintf(gp, "set terminal pngcairo size 1200,800\n");
  fprintf(gp, "set output '%s'\n", filename.c_str());
  fprintf(gp, "set title 'Bifurcation diagram - Cubic Ikeda DDE'\n");
  fprintf(gp, "set xlabel 'a'\n");
  fprintf(gp, "set ylabel 'x(t-{/Symbol t}) at section'\n");
  fprintf(gp, "set pointsize 1.5\n");

  // this gave me confirmation that for a = 1.535 there is the 2-periodic orbit
  fprintf(gp, "set arrow from %f, graph 0 to %f, graph 1 nohead lc rgb 'red' lw 1.5\n", 1.5385, 1.5385);

  fprintf(gp, "plot '-' with dots lc rgb '#000000' notitle\n");

  double aIncrease = (aEnd-aStart)/aFrequency;
  DVector returnPoint(x);
  int N = returnPoint.dimension() - 1;
  for (double a = aStart; a <= aEnd; a += aIncrease) {
    cout << "plotting bifurcation diagram for a = " <<  a << endl;
    CubicIkeda.setParameters({a});
    try {
      for (int k = 0; k < 3. / 4. * noSteps; k++) {
        returnPoint = pm(returnPoint);
      }
      for (int k = 3. / 4. * noSteps; k < noSteps; k++) {
        returnPoint = pm(returnPoint);
        fprintf(gp, "%f %f\n", a, returnPoint[N]);
      }
    } catch (exception &e) {
      cout << "Exception at a=" << a << ": " << e.what() << "\n";
      returnPoint = x;
      continue; // we skip this a and try the next one
    }
  }

  fprintf(gp, "e\n"); // end inline data — PNG is written here
  fflush(gp);

  pclose(gp);
}
