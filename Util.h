//
// Created by Dennis Goldfarb on 1/2/18.
//

#ifndef DEMIX_UTIL_H
#define DEMIX_UTIL_H

#include <OpenMS/MATH/MISC/MathFunctions.h>
#include "OpenMS/CONCEPT/Constants.h"

#include "MassToleranceUnit.h"
#include "PrecursorTargetOptions.h"

class Util {

public:
    static double getTol(double mz, double tol, MassToleranceUnit unit);

    static bool withinTol(double mz1, double mz2, double tol, MassToleranceUnit unit);

    static int compareWithTol(double mz1, double mz2, double tol, MassToleranceUnit unit);

    static void populateOptionsForIsolationWindow(int minCharge, int maxCharge, int maxIsotope,
                                                        double *charge2prob, double isolationWidth,
                                                        double isolationCenter, PrecursorTargetOptions &options);
};


#endif //DEMIX_UTIL_H
