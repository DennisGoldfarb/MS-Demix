//
// Created by Dennis Goldfarb on 12/29/17.
//

#include "PrecursorTargetOption.h"

void PrecursorTargetOption::fillPrecursorIsotopes() {
    for (OpenMS::UInt i = minIso; i <= maxIso; ++i)
    {
        precursorIsotopes.insert(i);
    }
}

bool PrecursorTargetOption::operator== (const PrecursorTargetOption& c) const
{
    return (charge == c.charge &&
            likelihood == c.likelihood &&
            maxIso == c.maxIso &&
            minIso == c.minIso &&
            minMass == c.minMass &&
            maxMass == c.maxMass);
}

bool PrecursorTargetOption::operator!= (const PrecursorTargetOption& c) const
{
    return !(*this == c);
}
