#include "integrateFuncs.h"

DTimeMap::SolutionCurve getSolutionCurve(DTimeMap timeMap, DVector start,
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

  return solution;
}

void getPoincareValues(DPoincareMap &pm, DVector &x, const string &filename) {
  ofstream csv(filename);
  csv << x << "\n";
  double returnTime = 0;
  DVector returnPoint = pm(x, returnTime);
  csv << returnPoint << "\n";
  for (int i = 0; i < 10; i++) {
    returnPoint = pm(returnPoint, returnTime);
    csv << returnPoint << "\n";
  }
  csv.close();
}

// gnuplot code written by Claude
void plotBifurcationDiagram(DMap &CubicIkeda, DPoincareMap &pm, DVector &x,
                            double aStart, double aEnd, double aIncrease,
                            int noSteps, string filename) {
  FILE *gp = popen("gnuplot", "w"); // no -persistent needed for file output

  fprintf(gp, "set terminal pngcairo size 1200,800\n");
  fprintf(gp, "set output '%s'\n", filename.c_str());
  fprintf(gp, "set title 'Bifurcation diagram - Cubic Ikeda DDE'\n");
  fprintf(gp, "set xlabel 'a'\n");
  fprintf(gp, "set ylabel 'x(t-{/Symbol t}) at section'\n");
  fprintf(gp, "set pointsize 1.5\n");
  fprintf(gp, "plot '-' with dots lc rgb '#000000' notitle\n");

  DVector returnPoint(x);
  int N = returnPoint.dimension() - 1;
  for (double a = aStart; a <= aEnd; a += aIncrease) {
    cout << "a = " <<  a << endl;
    CubicIkeda.setParameters({a});
    try {
      for (int k = 0; k < 3. / 4. * noSteps; k++) {
        double returnTime = 0.;
        returnPoint = pm(returnPoint, returnTime);
      }
      for (int k = 3. / 4. * noSteps; k < noSteps; k++) {
        double returnTime = 0.;
        returnPoint = pm(returnPoint, returnTime);
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
