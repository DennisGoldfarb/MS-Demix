//
// Created by Dennis Goldfarb on 12/28/17.
//

#include "MultiplexedScan.h"

MultiplexedScan::MultiplexedScan(OpenMS::MSSpectrum<OpenMS::Peak1D> scan)
{
    for (auto peak : scan)
    {
        mzData.push_back(peak.getMZ());
        intData.push_back(peak.getIntensity());
    }
}
