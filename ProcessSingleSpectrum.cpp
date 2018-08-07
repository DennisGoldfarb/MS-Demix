//
// Created by Dennis Goldfarb on 2/25/18.
//
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
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
    double monoMass, intensity;

    while (std::getline(in, line))
    {
        std::istringstream iss(line);

        iss >> label;
        if (label == 'P')
        {
	  iss >> monoMass >> charge >> intensity;
	    if (charge > 1)
	    {
	      HardklorEntry precursor(scanID, charge, monoMass, maxIso, intensity);
                precursors.push_back(precursor);
            }
        } else {
            iss >> scanID;
        }
    }

    return precursors;
}

double getTotalAbundance(OpenMS::IsotopeDistribution &id, double basePeakIntensity, int firstIsotope, int lastIsotope)
{
  double totalAbundance = 0.0;
  double basePeak = 0.0;
  for (int i = 0; i < id.size(); ++i)
    {
      basePeak = std::max(basePeak, (double) id.getContainer()[i].getIntensity());
    }

  for (int i = firstIsotope; i <= lastIsotope; ++i)
    {
      totalAbundance += basePeakIntensity * id.getContainer()[i].getIntensity() / basePeak;
    }

  return totalAbundance;
}

int getPrecursorOptionsFromHardklor(int minCharge, int maxCharge, int maxIsotope,
                                    double isolationWidth, double isolationCenter,
				    int previousMS1ScanID, int nextMS1ScanID, PrecursorTargetOptions &options, 
				    std::vector<HardklorEntry> &precursors, int MS2scanID, int lastP)
{
    double isoLeft = isolationCenter - (isolationWidth / 2);
    double isoRight = isolationCenter + (isolationWidth / 2);

    int lastScan = 0;
    // find flanking MS1 scans
    for (int i = lastP; i < precursors.size(); ++i)
    {
        HardklorEntry precursor = precursors[i];
	
        if (precursor.scanID == previousMS1ScanID || precursor.scanID == nextMS1ScanID)
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

                if (firstIsotope == 0 && lastIsotope ==0)
		{
		    continue;
		}

		// Determine abundance
		OpenMS::CoarseIsotopePatternGenerator gen(lastIsotope+1);
		OpenMS::IsotopeDistribution id = gen.estimateFromPeptideWeight(monoMass);
		double abundance = getTotalAbundance(id, precursor.intensity, firstIsotope, lastIsotope);

		PrecursorTargetOption option(monoMass, monoMass, precursor.charge, firstIsotope, lastIsotope, 1.0, abundance);
		
		std::cout << "MATCH: " << previousMS1ScanID << " " << isolationCenter << " " << isoLeft << " " << isoRight << std::endl;
		std::cout << precursor.scanID << " " << precursor.minMz << " " << precursor.maxMz << " " << firstIsotope << " " << lastIsotope << " " << precursor.charge << " " << abundance <<std::endl;

		if (!options.hasOption(option))
		{
		    options.addOption(option);
		} 
		else 
		{
		    options.addAbundance(option);
		    std::cout << "ADDING ABUNDANCE" << std::endl; 
		}
            }
        }

	else if (precursor.scanID > nextMS1ScanID)
	{
	  return lastP;
	}

	if (precursor.scanID != lastScan && precursor.scanID < nextMS1ScanID)
          {
            lastScan = precursor.scanID;
            lastP = i;
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

    int previousMS1ScanID = 1, nextMS1ScanID=1, lastP = 0;

    for (int i = scanStart; i < map.size(); i+=numJobs)
    {
        for (int j = i; j >= 0; --j)
	{
	  OpenMS::MSSpectrum scan = map.getSpectrum(j);
	  if (scan.getMSLevel() == 1)
	  {
	      previousMS1ScanID = j+1;
	      break;
	  }
	}

	for (int j = i; j < map.size(); ++j)
	  {
	    OpenMS::MSSpectrum scan = map.getSpectrum(j);
	    if (scan.getMSLevel() == 1)
	      {
		nextMS1ScanID = j+1;
		break;
	      }
          }
	

        OpenMS::MSSpectrum scan = map.getSpectrum(i);

        if (scan.getMSLevel() == 2)
        {
	  std::cout << "MS2 scan: " << i+1 << " " << previousMS1ScanID << " " << nextMS1ScanID << std::endl;

            const OpenMS::Precursor precursorInfo = scan.getPrecursors()[0];

            double isolationWidth = precursorInfo.getIsolationWindowUpperOffset() + precursorInfo.getIsolationWindowLowerOffset();
            double isolationCenter = precursorInfo.getMZ();
            double monoMass = isolationCenter * precursorInfo.getCharge();
            double isotopeSpacing = 1.0 / precursorInfo.getCharge();
            int maxIsolatedIso = isolationWidth / 2 / isotopeSpacing;

            PrecursorTargetOption target(monoMass, monoMass, precursorInfo.getCharge(), 0, maxIsolatedIso, 1.0, 1.0);

            PrecursorTargetOptions options;
            options.addOption(target);

            lastP = getPrecursorOptionsFromHardklor(minZ, maxZ, maxIso, isolationWidth, isolationCenter, previousMS1ScanID, nextMS1ScanID, options, precursors, i+1, lastP);
            //Util::populateOptionsForIsolationWindow(minZ, maxZ, maxIso, charge2prob, isolationWidth, isolationCenter, options);

            if (options.key2option.size() > 1)
            {
                std::cout << "Chimeric Scan: " << i << " " << options.key2option.size() << std::endl;
		//continue;
	    }
	    else 
	      {
		std::cout << "Single Scan: " << i << " " << options.key2option.size() << std::endl;
	      }

	    if (options.key2option.size() == 0)
	      {
		options.addOption(target);
	      }

	    NNLSModel model(scan, options, 20, MassToleranceUnit::PPM);

	    model.writeModel(outPath, std::to_string(i+1));

	    writeScan(scan, outPath, i+1, isolationCenter, precursorInfo.getCharge());

	}
        else
        {
            previousMS1ScanID = i+1;
        }
    }

    return 0;
}
