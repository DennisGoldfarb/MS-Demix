//
// Created by Dennis Goldfarb on 2/25/18.
//
#include <iostream>

#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include "Util.h"
#include "NNLSModel.h"

void usage()
{
    std::cout << "usage: ProcessSingleSpectrum mzML_file scan_ID minZ maxZ maxIso output_directory" << std::endl;
    std::cout << "\tmzML_file: path to input .mzML file " << std::endl;
    std::cout << "\tscan_ID: scan_ID to process." << std::endl;
    std::cout << "\tminZ: min precursor charge state to consider." << std::endl;
    std::cout << "\tmaxZ: max precursor charge state to consider." << std::endl;
    std::cout << "\tmaxIso: max precursor isotope to consider." << std::endl;
    std::cout << "\toutput_directory: path to output files" << std::endl;
}

int main(int argc, char * argv[])
{
    //check for correct number of command line arguments
    if (argc != 7) {
        usage();
        return 0;
    }

    std::string mzMLPath = argv[1];
    int scanID = atoi(argv[2]);
    int minZ = atoi(argv[3]);
    int maxZ = atoi(argv[4]);
    int maxIso = atoi(argv[5]);
    std::string outPath = argv[6];

    double charge2prob[7] = {0.0, 0.01, 0.5, 0.5, 0.05, 0.03, 0.02};

    OpenMS::IndexedMzMLFileLoader mzMLDataFile;
    OpenMS::OnDiscPeakMap map;
    
    mzMLDataFile.load(mzMLPath, map);

    OpenMS::MSSpectrum scan = map.getSpectrum(scanID);

    if (scan.getMSLevel() == 2) {

        const OpenMS::Precursor precursorInfo = scan.getPrecursors()[0];

        double isolationWidth = precursorInfo.getIsolationWindowUpperOffset() + precursorInfo.getIsolationWindowLowerOffset();
        double isolationCenter = precursorInfo.getMZ();
        double monoMass = (isolationCenter - OpenMS::Constants::C13C12_MASSDIFF_U) * precursorInfo.getCharge();
        double isotopeSpacing = 1.0 / precursorInfo.getCharge();
        int maxIsolatedIso = isolationWidth/2 / isotopeSpacing;

        PrecursorTargetOption target(monoMass, monoMass, precursorInfo.getCharge(), 0, maxIsolatedIso, 1.0);

        PrecursorTargetOptions options;
        options.addOption(target);

        Util::populateOptionsForIsolationWindow(minZ, maxZ, maxIso, charge2prob, isolationWidth, isolationCenter, options);

        NNLSModel model(scan, options, 20, MassToleranceUnit::PPM);

        model.writeModel("/Users/dennisg/Documents/Manuscripts/Deconvolution/Src/", argv[2]);
    }

    return 0;
}
