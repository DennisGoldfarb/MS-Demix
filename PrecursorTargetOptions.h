//
// Created by Dennis Goldfarb on 12/29/17.
//

#ifndef DEMIX_PRECURSORTARGETOPTIONS_H
#define DEMIX_PRECURSORTARGETOPTIONS_H

#include <map>
#include "PrecursorTargetOption.h"

class PrecursorTargetOptionsKey
{
public:
    int z;
    int firstIso;
    int lastIso;

    PrecursorTargetOptionsKey() {};

    PrecursorTargetOptionsKey(int z, int firstIso, int lastIso) : z(z), firstIso(firstIso), lastIso(lastIso) {};

    bool operator<( PrecursorTargetOptionsKey const& o) const
    {
        if (z < o.z) return true;
        if (firstIso < o.firstIso) return true;
        return lastIso < o.lastIso;
    }

};



class PrecursorTargetOptions {

public:

    PrecursorTargetOptions() {};

    void addOption(PrecursorTargetOption o);
    bool hasOption(PrecursorTargetOption o);

    void addAbundance(PrecursorTargetOption o);

    int getMaxIsotope();
    int getMaxCharge();

    std::map<PrecursorTargetOptionsKey, PrecursorTargetOption> key2option;

    PrecursorTargetOptions(const PrecursorTargetOptions &o);

    PrecursorTargetOptions& operator=(const PrecursorTargetOptions &o);
};


#endif //DEMIX_PRECURSORTARGETOPTIONS_H
