#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <sstream>

#include "OpenMS/CONCEPT/Constants.h"
#include "MultiplexedScan.h"
#include "PrecursorTargetOption.h"
#include "PrecursorTargetOptions.h"
#include "NNLSModel.h"


MultiplexedScan readSpectrumFile(std::string path)
{
    std::ifstream spectrumFile(path);

    // file starts with the isolation width
    double isolationWidth, centerA, centerB;
    int startIsoA, endIsoA, startIsoB, endIsoB, chargeA, chargeB;
    spectrumFile >> isolationWidth >> startIsoA >> endIsoA >> startIsoB >> endIsoB >> chargeA >> chargeB;

    // ignore header lines
    std::string line;
    for (int i = 0; i < 8; ++i)
    {
        std::getline(spectrumFile, line);
        if (i == 3)
        {
            std::istringstream iss(line);
            std::vector<std::string> results((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

            centerA = std::atof(results[7].substr(0, results[7].find("@")).c_str());
            centerB = std::atof(results[8].substr(0, results[8].find("@")).c_str());
        }
    }

    // parse mz - intensity pairs
    std::vector<double> mzData, intData;
    double mz, intensity;
    while (spectrumFile >> mz >> intensity)
    {
        mzData.push_back(mz);
        intData.push_back(intensity);
    }


    OpenMS::AASequence precursorA = OpenMS::AASequence::fromString("DRVYIHPFHL");
    OpenMS::AASequence precursorB = OpenMS::AASequence::fromString("E[-18.010565]LYENKPRRPYIL");

    MultiplexedScan scan(precursorA, precursorB, chargeA, chargeB, startIsoA, endIsoA, startIsoB, endIsoB,
                         isolationWidth, centerA, centerB, 130, 1300, mzData, intData);

    return scan;
}

int main(int argc, char * argv[])
{
    int mode = std::atoi(argv[3]);
    MultiplexedScan scan = readSpectrumFile(std::string(argv[1]) + std::string(argv[2]) + ".txt");

    int minZ = 2;
    int maxZ = 3;
    int maxIso = 3;

    double charge2prob[7] = {0.0, 0.01, 0.5, 0.5, 0.05, 0.03, 0.02};

    PrecursorTargetOptions optionsA, optionsB;

    PrecursorTargetOption targetA(scan.seqA.getMonoWeight(), scan.seqA.getMonoWeight(), scan.chargeA, scan.minIsoA, scan.maxIsoA, 1.0, 1.0);
    PrecursorTargetOption targetB(scan.seqB.getMonoWeight(), scan.seqB.getMonoWeight(), scan.chargeB, scan.minIsoB, scan.maxIsoB, 1.0, 1.0);


    if (mode == 0) // only use targets
    {
        optionsA.addOption(targetA);
        optionsB.addOption(targetB);
    }
    else if (mode == 1) // use targets and options
    {
        optionsA.addOption(targetA);
        optionsB.addOption(targetB);

        Util::populateOptionsForIsolationWindow(minZ, maxZ, maxIso, charge2prob, scan.isolationWidth, scan.isolationCenterA, optionsA);
        Util::populateOptionsForIsolationWindow(minZ, maxZ, maxIso, charge2prob, scan.isolationWidth, scan.isolationCenterB, optionsB);
    }
    else if (mode == 2) // only use options
    {
        Util::populateOptionsForIsolationWindow(minZ, maxZ, maxIso, charge2prob, scan.isolationWidth, scan.isolationCenterA, optionsA);
        Util::populateOptionsForIsolationWindow(minZ, maxZ, maxIso, charge2prob, scan.isolationWidth, scan.isolationCenterB, optionsB);
    }

    std::vector<PrecursorTargetOptions> allOptions;
    allOptions.push_back(optionsA);
    allOptions.push_back(optionsB);

    NNLSModel model(scan, allOptions, 20, MassToleranceUnit::PPM);

    model.writeModel("/Users/dennisg/Documents/Manuscripts/Deconvolution/Src/", argv[2]);

    return 0;
}
