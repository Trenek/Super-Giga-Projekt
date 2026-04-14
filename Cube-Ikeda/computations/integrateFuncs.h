#include "capd/capdlib.h"
#include <iomanip>
#include <iostream>
using namespace capd;
using namespace std;

/// @brief  Returns a curve (function of (double) time argument)
/// @param pm       Poincaré map
/// @param x        initial state vector.
/// @param time     value of time for integration
/// @param skip     value of time skipped before saving data into csv
/// @param filename path to the output csv file (must be non-empty)
DTimeMap::SolutionCurve getSolutionCurve(DTimeMap timeMap, DVector x,
                                         double time, double skip,
                                         string filename);

/// @brief Computes first 10 values of a Poincare map od Cubic Ikeda function
/// @param pm       Poincaré map
/// @param x        initial state vector.
/// @param filename path to the output txt file (must be non-empty)
void getPoincareValues(DPoincareMap &pm, DVector &x, const string &filename);

/// @brief Uses gnuplot to plot a bifurcation diagram for the Poincaré map
///        of the cubic Ikeda function.
/// @param dmap     map used for integration.
/// @param pm       Poincaré map defined on the solver of dmap and a section.
/// @param x        initial state vector.
/// @param aStart   starting value of parameter a.
/// @param aEnd     ending value of parameter a.
/// @param noSteps  number of iterations computed for each value of a.
/// @param filename path to the output PNG file (must be non-empty).
void plotBifurcationDiagram(DMap &dmap, DPoincareMap &pm, DVector &x,
                            double aStart, double aEnd, double aFrequency,
                            int noSteps, string filename);
//  @param aIncrease  increment of parameter a between successive evaluations.