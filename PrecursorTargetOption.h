//
// Created by Dennis Goldfarb on 12/29/17.
//

#ifndef DEMIX_PRECURSORTARGET_H
#define DEMIX_PRECURSORTARGET_H

#include <set>
#include <OpenMS/CONCEPT/Types.h>
#include "OpenMS/CONCEPT/Constants.h"

class PrecursorTargetOption {

public:

    PrecursorTargetOption() {};

 PrecursorTargetOption(double minMass, double maxMass, int charge, int minIso, int maxIso, double likelihood, double abundance) :
    minMass(minMass), maxMass(maxMass), charge(charge), minIso(minIso), maxIso(maxIso), likelihood(likelihood), abundance(abundance)
    {
        fillPrecursorIsotopes();
    };

    bool operator==(const PrecursorTargetOption&) const;
    bool operator!=(const PrecursorTargetOption&) const;

    double getMinMz();

    double minMass;
    double maxMass;
    int charge;
    int minIso;
    int maxIso;
    double likelihood;
    double abundance;

    std::set<OpenMS::UInt> precursorIsotopes;

private:
    void fillPrecursorIsotopes();

};

#endif //DEMIX_PRECURSORTARGET_H
