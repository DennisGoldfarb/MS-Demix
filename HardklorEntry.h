//
// Created by Dennis Goldfarb on 3/15/18.
//

#ifndef DEMIX_HARDKLORENTRY_H
#define DEMIX_HARDKLORENTRY_H

#include <iostream>

#include "OpenMS/CONCEPT/Constants.h"

class HardklorEntry {

public:
    HardklorEntry(int scanID, int charge, double monoMass, int maxIso)
    {
        this->scanID = scanID;
        this->charge = charge;
        minMz = (monoMass / charge) + OpenMS::Constants::PROTON_MASS_U;
        maxMz = minMz + (maxIso / (double) charge);
    }

    int scanID;
    int charge;
    double minMz;
    double maxMz;
};


#endif //DEMIX_HARDKLORENTRY_H
