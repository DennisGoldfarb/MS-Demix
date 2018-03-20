//
// Created by Dennis Goldfarb on 2/25/18.
//
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include "Util.h"
#include "NNLSModel.h"
#include "HardklorEntry.h"

void usage() {
    std::cout << "usage: ProcessSingleSpectrum bullseye_file isolationWidth scanStart numJobs minZ maxZ maxIso output_directory" << std::endl;

    std::cout << "\tbullseye_file: path to input bullseye .ms2 file " << std::endl;
    std::cout << "\tisolationWidth: isolation width of each MS2 scan." << std::endl;
    std::cout << "\tscanStart: first scan_ID to process." << std::endl;
    std::cout << "\tnumJobs: number of current jobs performed." << std::endl;
    std::cout << "\tminZ: min precursor charge state to consider." << std::endl;
    std::cout << "\tmaxZ: max precursor charge state to consider." << std::endl;
    std::cout << "\tmaxIso: max precursor isotope to consider." << std::endl;
    std::cout << "\toutput_directory: path to output files" << std::endl;
}


void writeScan(std::string outPath, int scanID, double monoMass, int z)
{
    std::ofstream out(outPath + std::to_string(scanID) + ".tab");

    out << scanID << "\t" << monoMass << "\t" << z << "\t0.0" << std::endl;

    out.close();
}

double getTotalAbundance(OpenMS::IsotopeDistribution &id, double basePeakIntensity, int firstIsotope, int lastIsotope)
{
    double totalAbundance = 0.0;
    double basePeak = 0.0;
    for (int i = 0; i < id.size(); ++i)
    {
        basePeak = std::max(basePeak, id.getContainer()[i].second);
    }

    for (int i = firstIsotope; i <= lastIsotope; ++i)
    {
        totalAbundance += basePeakIntensity * id.getContainer()[i].second / basePeak;
    }

    return totalAbundance;
}

void addOption(PrecursorTargetOptions &options, double isolationCenter, double isolationWidth, int z,
               double monoNeutralMass, int maxIsotope, double intensity, int scanID)
{
    if ( z == 1) return;
    double isoLeft = isolationCenter - (isolationWidth / 2);
    double isoRight = isolationCenter + (isolationWidth / 2);

    double isotopeSpacing = OpenMS::Constants::C13C12_MASSDIFF_U / z;
    double minMz = ((monoNeutralMass + (OpenMS::Constants::PROTON_MASS_U * (z-1))) / z); //+ OpenMS::Constants::PROTON_MASS_U;
    double monoMass = minMz * z;
    int firstIsotope = std::numeric_limits<int>::max();
    int lastIsotope = std::numeric_limits<int>::min();

    for (int isotopeIndex = 0; isotopeIndex <= maxIsotope; ++isotopeIndex) {
        double mz = minMz + (isotopeSpacing * isotopeIndex);

        if (isoLeft < mz && mz < isoRight) {
            if (isotopeIndex < firstIsotope) firstIsotope = isotopeIndex;
            if (isotopeIndex > lastIsotope) lastIsotope = isotopeIndex;
        }
    }

    if (firstIsotope == std::numeric_limits<int>::max())
    {
      std::cout << "UNMATCHED! " << scanID << " " << monoNeutralMass << " " << minMz << " " << monoMass << " " << firstIsotope << " " << lastIsotope
		<< " " << z  << std::endl;
        return;
    }

    // Determine abundance
    OpenMS::IsotopeDistribution id(maxIsotope + 1);
    id.estimateFromPeptideWeight(monoMass);
    double abundance = getTotalAbundance(id, intensity, firstIsotope, lastIsotope);

    PrecursorTargetOption option(monoMass, monoMass, z, firstIsotope, lastIsotope, 1.0, abundance);

    std::cout << "MATCH: " << isolationCenter << " " << isoLeft << " " << isoRight << std::endl;
    std::cout << scanID << " " << monoNeutralMass << " " << minMz << " " << monoMass << " " << firstIsotope << " " << lastIsotope
              << " " << z << " " << abundance << std::endl;

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


void parseBullseye(std::string bullseyePath, double isolationWidth, int scanStart, int numJobs,
                   int minZ, int maxZ, int maxIso, std::string outPath)
{
    std::ifstream in(bullseyePath);

    std::string line, label;
    char symbol;
    int scanID, charge;
    double isolation_center, monoNeutralMass, rt, abundance, intensity, mz;
    std::vector<double> mzData, intData;
    PrecursorTargetOptions options;

    while (std::getline(in, line))
    {
        std::istringstream iss(line);
        if (!std::isdigit(iss.peek()))
        {
	    iss >> symbol;

            if (symbol == 'S')
            {
                if (mzData.size() > 0 && (scanID % numJobs) == (scanStart - 1) )
                {
                    // create multiplexed scan object
                    if (options.key2option.size() > 1)
                    {
                        std::cout << "Chimeric Scan: " << scanID << " " << options.key2option.size() << std::endl;
                    } 
		    else if (options.key2option.size() == 1)
                    {
                        std::cout << "Single Scan: " << scanID << " " << options.key2option.size() << std::endl;
                    }
		    if (options.key2option.size() > 0)
		    {
                        MultiplexedScan scan;
                        scan.intData = intData;
                        scan.mzData = mzData;
                        scan.isolationWidth = isolationWidth;


                        NNLSModel model(scan, options, 20, MassToleranceUnit::PPM);

                        model.writeModel(outPath, std::to_string(scanID));

                        writeScan(outPath, scanID, isolation_center, charge);
		    } else {
		      std::cout << "No Scan: " << scanID << " " << options.key2option.size() << std::endl;
		    }
                }
                iss >> scanID >> scanID >> isolation_center;
                mzData.clear();
                intData.clear();
                options.key2option.clear();
            }
	    else if ((scanID % numJobs) == (scanStart - 1))
	    {
                if (symbol == 'I')
                {
                    iss >> label;
                    if (label == "EZ")
                    {
                        iss >> charge >> monoNeutralMass >> rt >> abundance;
                        // create precursor option
                        addOption(options, isolation_center, isolationWidth, charge, monoNeutralMass, maxIso, abundance, scanID);
                    }
                }
                else if (symbol == 'Z')
                {
                    if (options.key2option.size() == 0)
                    {
                        iss >> charge >> monoNeutralMass;
                        // create precursor option
                        addOption(options, isolation_center, isolationWidth, charge, monoNeutralMass, maxIso, 1.0, scanID);
                    }
                }
	    }
        }
        else if ((scanID % numJobs) == (scanStart - 1))
        {
            iss >> mz >> intensity;
            mzData.push_back(mz);
            intData.push_back(intensity);
        }
    }
}

int main(int argc, char *argv[]) {
    //check for correct number of command line arguments
    if (argc != 9) {
        usage();
        return -1;
    }

    std::string bullseyePath = argv[1];
    double isolationWidth = atof(argv[2]);
    int scanStart = atoi(argv[3]);
    int numJobs = atoi(argv[4]);
    int minZ = atoi(argv[5]);
    int maxZ = atoi(argv[6]);
    int maxIso = atoi(argv[7]);
    std::string outPath = argv[8];

    parseBullseye(bullseyePath, isolationWidth, scanStart, numJobs, minZ, maxZ, maxIso, outPath);

    return 0;
}
