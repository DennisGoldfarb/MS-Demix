//
// Created by Dennis Goldfarb on 3/15/18.
//

#ifndef DEMIX_HARDKLORENTRY_H
#define DEMIX_HARDKLORENTRY_H

#include <iostream>

#include "OpenMS/CONCEPT/Constants.h"

class HardklorEntry {

public:
  HardklorEntry(int scanID, int charge, double monoMass, int maxIso, double intensity)
    {
        this->scanID = scanID;
        this->charge = charge;
        minMz = (monoMass / charge) + OpenMS::Constants::PROTON_MASS_U;
        maxMz = minMz + (maxIso * OpenMS::Constants::C13C12_MASSDIFF_U / (double) charge);
	this->intensity = intensity;
    }

    int scanID;
    int charge;
    double minMz;
    double maxMz;
    double intensity;
};


#endif //DEMIX_HARDKLORENTRY_H
