//
// Created by Dennis Goldfarb on 12/29/17.
//

#ifndef DEMIX_PRECURSORTARGET_H
#define DEMIX_PRECURSORTARGET_H

#include <set>
#include <OpenMS/CONCEPT/Types.h>

class PrecursorTargetOption {

public:

    PrecursorTargetOption() {};

    PrecursorTargetOption(double minMass, double maxMass, int charge, int minIso, int maxIso, double likelihood) :
            minMass(minMass), maxMass(maxMass), charge(charge), minIso(minIso), maxIso(maxIso), likelihood(likelihood)
    {
        fillPrecursorIsotopes();
    };

    bool operator==(const PrecursorTargetOption&) const;
    bool operator!=(const PrecursorTargetOption&) const;


    double minMass;
    double maxMass;
    int charge;
    int minIso;
    int maxIso;
    double likelihood;

    std::set<OpenMS::UInt> precursorIsotopes;

private:
    void fillPrecursorIsotopes();

};

#endif //DEMIX_PRECURSORTARGET_H
