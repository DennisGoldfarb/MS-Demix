//
// Created by Dennis Goldfarb on 1/5/18.
//

#ifndef DEMIX_DICTIONARYELEMENT_H
#define DEMIX_DICTIONARYELEMENT_H


#include <vector>
#include "PrecursorTargetOption.h"

class DictionaryElement {

public:

    DictionaryElement(PrecursorTargetOption option, int charge, int offset, double monoMz,
                      std::vector<double> intData, std::vector<double> mzData, int index) :
            precursorDetails(option), charge(charge), offset(offset), monoMz(monoMz),
            intData(intData) , mzData(mzData), index(index) {} ;

    std::vector<double> intData;
    std::vector<double> mzData;

    PrecursorTargetOption precursorDetails;
    int charge;
    int offset;
    double monoMz;
    int index;
};


#endif //DEMIX_DICTIONARYELEMENT_H
