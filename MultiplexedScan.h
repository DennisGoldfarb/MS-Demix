//
// Created by Dennis Goldfarb on 12/28/17.
//

#ifndef DEMIX_SCAN_H
#define DEMIX_SCAN_H

#include <vector>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include "OpenMS/CHEMISTRY/AASequence.h"


class MultiplexedScan {

public :

    MultiplexedScan() {};

    MultiplexedScan(OpenMS::MSSpectrum<OpenMS::Peak1D> scan);

    MultiplexedScan(OpenMS::AASequence seqA, OpenMS::AASequence seqB,
                    int chargeA, int chargeB,
                    int minIsoA, int maxIsoA,
                    int minIsoB, int maxIsoB,
                    double isolationWidth,
                    double isolationCenterA, double isolationCenterB,
                    double scanMzStart, double scanMzEnd,
                    std::vector<double> mzData, std::vector<double> intData) :
            seqA(seqA), seqB(seqB), chargeA(chargeA), chargeB(chargeB),
            minIsoA(minIsoA), maxIsoA(maxIsoA), isolationWidth(isolationWidth),
            isolationCenterA(isolationCenterA), isolationCenterB(isolationCenterB),
            scanMzStart(scanMzStart), scanMzEnd(scanMzEnd),
            minIsoB(minIsoB), maxIsoB(maxIsoB), mzData(mzData), intData(intData)
    {};

    OpenMS::AASequence seqA;
    OpenMS::AASequence seqB;

    int chargeA;
    int chargeB;

    int minIsoA;
    int maxIsoA;
    int minIsoB;
    int maxIsoB;

    double isolationWidth;

    double isolationCenterA;
    double isolationCenterB;

    double scanMzStart;
    double scanMzEnd;

    std::vector<double> mzData;
    std::vector<double> intData;


};


#endif //DEMIX_SCAN_H
