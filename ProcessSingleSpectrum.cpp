//
// Created by Dennis Goldfarb on 2/25/18.
//
#include <iostream>

#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include "Util.h"
#include "NNLSModel.h"

void usage() {
    std::cout << "usage: ProcessSingleSpectrum mzML_file scan_ID numScans minZ maxZ maxIso output_directory"
              << std::endl;
    std::cout << "\tmzML_file: path to input .mzML file " << std::endl;
    std::cout << "\tscan_ID: first scan_ID to process." << std::endl;
    std::cout << "\tnumScans: number of scans to process." << std::endl;
    std::cout << "\tminZ: min precursor charge state to consider." << std::endl;
    std::cout << "\tmaxZ: max precursor charge state to consider." << std::endl;
    std::cout << "\tmaxIso: max precursor isotope to consider." << std::endl;
    std::cout << "\toutput_directory: path to output files" << std::endl;
}

void writeScan(OpenMS::MSSpectrum &scan, std::string outPath, int scanID, double monoMass, int z)
{
    std::ofstream out(outPath + std::to_string(scanID) + ".tab");

    out << scan.getNativeID() << "\t" << scan.getNativeID() << "\t" << scanID << "\t" << monoMass << "\t" << z << "\t" << scan.getRT() << std::endl;

    out.close();
}

int main(int argc, char *argv[]) {
    //check for correct number of command line arguments
    if (argc != 8) {
        usage();
        return -1;
    }

    std::string mzMLPath = argv[1];
    int scanID = atoi(argv[2]);
    int numScans = atoi(argv[3]);
    int minZ = atoi(argv[4]);
    int maxZ = atoi(argv[5]);
    int maxIso = atoi(argv[6]);
    std::string outPath = argv[7];

    double charge2prob[7] = {0.0, 0.01, 0.5, 0.5, 0.05, 0.03, 0.02};

    //OpenMS::MzMLFile mzMLDataFile;
    //OpenMS::PeakMap map;
    OpenMS::IndexedMzMLFileLoader mzMLDataFile;
    OpenMS::OnDiscPeakMap map;

    mzMLDataFile.load(mzMLPath, map);

    for (int i = scanID; i < scanID + numScans; ++i)
    {
        OpenMS::MSSpectrum scan = map.getSpectrum(i);

        if (scan.getMSLevel() == 2)
        {
            std::cout << "MS2 scan: " << i << std::endl;

            const OpenMS::Precursor precursorInfo = scan.getPrecursors()[0];

            double isolationWidth = precursorInfo.getIsolationWindowUpperOffset() + precursorInfo.getIsolationWindowLowerOffset();
            double isolationCenter = precursorInfo.getMZ();
            double monoMass = (isolationCenter - OpenMS::Constants::C13C12_MASSDIFF_U) * precursorInfo.getCharge();
            double isotopeSpacing = 1.0 / precursorInfo.getCharge();
            int maxIsolatedIso = isolationWidth / 2 / isotopeSpacing;

            PrecursorTargetOption target(monoMass, monoMass, precursorInfo.getCharge(), 0, maxIsolatedIso, 1.0);

            PrecursorTargetOptions options;
            options.addOption(target);

            Util::populateOptionsForIsolationWindow(minZ, maxZ, maxIso, charge2prob, isolationWidth, isolationCenter, options);

            NNLSModel model(scan, options, 20, MassToleranceUnit::PPM);

            model.writeModel(outPath, std::to_string(i));
            writeScan(scan, outPath, i, monoMass, precursorInfo.getCharge());
        }
    }

    return 0;
}
