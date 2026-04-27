#include "capd/capdlib.h"
#include <iomanip>
#include <iostream>
using namespace capd;
using namespace std;

DVector getCandidate(DMap& map, DPoincareMap& pm, int N, double aStart, double aEnd);
void plotRectangles(const IVector& rect1, const IVector& rect2);