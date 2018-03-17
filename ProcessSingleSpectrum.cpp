//
// Created by Dennis Goldfarb on 2/25/18.
//
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <OpenMS/FORMAT/IndexedMzMLFileLoader.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include "Util.h"
#include "NNLSModel.h"
#include "HardklorEntry.h"

void usage() {
    std::cout << "usage: ProcessSingleSpectrum mzML_file hardklor_file scanStart numJobs minZ maxZ maxIso output_directory"
              << std::endl;
    std::cout << "\tmzML_file: path to input .mzML file " << std::endl;
    std::cout << "\thardklor_file: path to input hardklor results file " << std::endl;
    std::cout << "\tscanStart: first scan_ID to process." << std::endl;
    std::cout << "\tnumJobs: number of current jobs performed." << std::endl;
    std::cout << "\tminZ: min precursor charge state to consider." << std::endl;
    std::cout << "\tmaxZ: max precursor charge state to consider." << std::endl;
    std::cout << "\tmaxIso: max precursor isotope to consider." << std::endl;
    std::cout << "\toutput_directory: path to output files" << std::endl;
}


void writeScan(OpenMS::MSSpectrum &scan, std::string outPath, int scanID, double monoMass, int z)
{
    std::ofstream out(outPath + std::to_string(scanID) + ".tab");

    out << scanID << "\t" << monoMass << "\t" << z << "\t" << scan.getRT() << std::endl;

    out.close();
}

std::vector<HardklorEntry> parseHardklor(std::string hardklorPath, int maxIso)
{
    std::ifstream in(hardklorPath);
    std::vector<HardklorEntry> precursors;

    std::string line;
    char label;
    int scanID, charge;
    double monoMass;

    while (std::getline(in, line))
    {
        std::istringstream iss(line);

        iss >> label;
        if (label == 'P')
        {
            iss >> monoMass >> charge;
	    if (charge > 1)
	    {
                HardklorEntry precursor(scanID, charge, monoMass, maxIso);
                precursors.push_back(precursor);
            }
        } else {
            iss >> scanID;
        }
    }

    return precursors;
}

int getPrecursorOptionsFromHardklor(int minCharge, int maxCharge, int maxIsotope,
                                     double isolationWidth, double isolationCenter,
				    int previousMS1ScanID, PrecursorTargetOptions &options, std::vector<HardklorEntry> &precursors, int MS2scanID, int lastP)
{
    double isoLeft = isolationCenter - (isolationWidth / 2);
    double isoRight = isolationCenter + (isolationWidth / 2);

    //bool targetFound = false;
    // find flanking MS1 scans
    for (int i = lastP; i < precursors.size(); ++i)
    {
        HardklorEntry precursor = precursors[i];

        if (precursor.scanID == previousMS1ScanID)
        {
            // is it in the isolation window
            if ( (isoLeft < precursor.minMz && precursor.minMz < isoRight) ||
                    (isoLeft < precursor.maxMz && precursor.maxMz < isoRight) ||
                    (precursor.minMz < isoLeft && isoRight < precursor.maxMz) )
            {
                double isotopeSpacing = OpenMS::Constants::C13C12_MASSDIFF_U / precursor.charge;
                double monoMass = precursor.minMz * precursor.charge;
                int firstIsotope = std::numeric_limits<int>::max();
                int lastIsotope = std::numeric_limits<int>::min();

                for (int isotopeIndex = 0; isotopeIndex <= maxIsotope; ++isotopeIndex)
                {
                    double mz = precursor.minMz + (isotopeSpacing * isotopeIndex);

                    if (isoLeft < mz && mz < isoRight)
                    {
                        if (isotopeIndex < firstIsotope) firstIsotope = isotopeIndex;
                        if (isotopeIndex > lastIsotope) lastIsotope = isotopeIndex;
                    }
		}
		PrecursorTargetOption option(monoMass, monoMass, precursor.charge, firstIsotope, lastIsotope, 1.0);

		//std::cerr << "DEBUG: " << precursor.scanID << " " << monoMass << " " << precursor.charge << " " << firstIsotope << " " << lastIsotope << " " << isoLeft << " " << isoRight << " " << precursor.minMz << " " << precursor.maxMz << std::endl;

		if (!options.hasOption(option))
		{
		    options.addOption(option);
		} 
		//else 
		//{
		//    targetFound = true;
		//}
            }
        }

	else if (precursor.scanID > previousMS1ScanID)
	{
	  //if (!targetFound) {
	  //  std::cout << "NOT FOUND: " << precursor.scanID << " " << MS2scanID << " " << previousMS1ScanID << std::endl;
	  //}
	  return i-500;
	}
    }
}

int main(int argc, char *argv[]) {
    //check for correct number of command line arguments
    if (argc != 9) {
        usage();
        return -1;
    }

    std::string mzMLPath = argv[1];
    std::string hardklorPath = argv[2];
    int scanStart = atoi(argv[3]);
    int numJobs = atoi(argv[4]);
    int minZ = atoi(argv[5]);
    int maxZ = atoi(argv[6]);
    int maxIso = atoi(argv[7]);
    std::string outPath = argv[8];

    double charge2prob[7] = {0.0, 0.01, 0.5, 0.5, 0.05, 0.03, 0.02};

    //OpenMS::MzMLFile mzMLDataFile;
    //OpenMS::PeakMap map;
    OpenMS::IndexedMzMLFileLoader mzMLDataFile;
    OpenMS::OnDiscPeakMap map;

    mzMLDataFile.load(mzMLPath, map);

    std::vector<HardklorEntry> precursors = parseHardklor(hardklorPath, maxIso);

    int previousMS1ScanID = 1, lastP = 0, numChimera = 0;

    for (int i = std::max(0, scanStart-50); i < scanStart; ++i)
    {
      OpenMS::MSSpectrum scan = map.getSpectrum(i);
      if (scan.getMSLevel() == 1)
	{
	  previousMS1ScanID = i+1;
	}
    } 

    for (int i = scanStart; i < map.size(); i+=numJobs)
    {
        OpenMS::MSSpectrum scan = map.getSpectrum(i);

        if (scan.getMSLevel() == 2)
        {
            std::cout << "MS2 scan: " << i+1 << std::endl;

            const OpenMS::Precursor precursorInfo = scan.getPrecursors()[0];

            double isolationWidth = precursorInfo.getIsolationWindowUpperOffset() + precursorInfo.getIsolationWindowLowerOffset();
            double isolationCenter = precursorInfo.getMZ();
            double monoMass = (isolationCenter - OpenMS::Constants::C13C12_MASSDIFF_U) * precursorInfo.getCharge();
            double isotopeSpacing = 1.0 / precursorInfo.getCharge();
            int maxIsolatedIso = isolationWidth / 2 / isotopeSpacing;

            PrecursorTargetOption target(monoMass, monoMass, precursorInfo.getCharge(), 0, maxIsolatedIso, 1.0);

            PrecursorTargetOptions options;
            //options.addOption(target);

            lastP = getPrecursorOptionsFromHardklor(minZ, maxZ, maxIso, isolationWidth, isolationCenter, previousMS1ScanID, options, precursors, i+1, lastP);
            //Util::populateOptionsForIsolationWindow(minZ, maxZ, maxIso, charge2prob, isolationWidth, isolationCenter, options);

            if (options.key2option.size() > 1)
            {
                std::cout << "Chimeric Scan: " << i << " " << options.key2option.size() << std::endl;

                NNLSModel model(scan, options, 20, MassToleranceUnit::PPM);

                model.writeModel(outPath, std::to_string(i+1));
                writeScan(scan, outPath, i+1, monoMass, precursorInfo.getCharge());
		
		numChimera++;
            }
        }
        else
        {
            previousMS1ScanID = i+1;
        }
    }

    return numChimera;
}
