//
// Created by Dennis Goldfarb on 5/15/18.
//
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iomanip>
#include "Fragment.h"
#include "Util.h"

std::vector<Fragment> getFragments(OpenMS::AASequence pep, std::string name)
{
    // create fragments
    std::vector<Fragment> fragments;

    fragments.push_back(Fragment(3, pep.getMonoWeight(OpenMS::Residue::ResidueType::Full, 3), "Precursor", name));
    fragments.push_back(Fragment(3, pep.getMonoWeight(OpenMS::Residue::ResidueType::Full, 3) - 18.0091422, "Precursor-H2O", name));
    fragments.push_back(Fragment(3, pep.getMonoWeight(OpenMS::Residue::ResidueType::Full, 3) - 17.0086343, "Precursor-NH3", name));

    for (int i = 1; i < pep.size(); ++i)
    {
        std::string zLabel = "";
        for (int z = 1; z <= 3; ++z)
        {
            zLabel += "+";
            fragments.push_back(Fragment(z, pep.getPrefix(i).getMonoWeight(OpenMS::Residue::ResidueType::BIon, z),
                                         "B"+std::to_string(i)+zLabel, name));
            fragments.push_back(Fragment(z, pep.getPrefix(i).getMonoWeight(OpenMS::Residue::ResidueType::BIon, z) - 18.0091422,
                                         "B"+std::to_string(i)+"-H2O"+zLabel, name));
            fragments.push_back(Fragment(z, pep.getPrefix(i).getMonoWeight(OpenMS::Residue::ResidueType::BIon, z) - 17.0086343,
                                         "B"+std::to_string(i)+"-NH3"+zLabel, name));

            fragments.push_back(Fragment(z, pep.getSuffix(i).getMonoWeight(OpenMS::Residue::ResidueType::YIon, z),
                                         "Y"+std::to_string(i)+zLabel, name));
            fragments.push_back(Fragment(z, pep.getSuffix(i).getMonoWeight(OpenMS::Residue::ResidueType::YIon, z) - 18.0091422,
                                         "Y"+std::to_string(i)+"-H2O"+zLabel, name));
            fragments.push_back(Fragment(z, pep.getSuffix(i).getMonoWeight(OpenMS::Residue::ResidueType::YIon, z) - 17.0086343,
                                         "Y"+std::to_string(i)+"-NH3"+zLabel, name));
        }
    }

    // sort by mono mz
    std::sort(fragments.begin(), fragments.end());

    return fragments;
}

void matchAndWrite(std::vector<Fragment> &fragments, std::vector<double> &mzData, std::vector<double> &intData)
{
    double tol = 20;
    int i = 0, j = 0;
    // for each peak, check if it's a fragment
    while (i < mzData.size() && j < fragments.size())
    {
        // do we have a match?
        std::string group = std::to_string(i);
        std::string label = "";
        std::string matched = "unmatched";

        if (Util::withinTol(mzData[i], fragments[j].monoMz, tol, MassToleranceUnit::PPM))
        {
            if (std::abs(mzData[i]-fragments[j].monoMz) < std::abs(mzData[i+1]-fragments[j].monoMz)) {
                label = fragments[j].name;
                matched = fragments[j].pep;
            }
        }

        // write
        std::cout << mzData[i]-0.001 << "\t" << 0 << "\t" << group << "\t" << matched << "\t" << "" << std::endl;
        std::cout << mzData[i] << "\t" << intData[i] << "\t" << group << "\t" << matched << "\t" << label << std::endl;
        std::cout << mzData[i]+0.001 << "\t" << 0 << "\t" << group << "\t" << matched << "\t" << "" << std::endl;

        ++i;
        while (!Util::withinTol(mzData[i], fragments[j].monoMz, tol, MassToleranceUnit::PPM) && mzData[i] >= fragments[j].monoMz)
        {
            ++j;
        }

    }

    while (i < mzData.size())
    {
        // write
        std::string group = std::to_string(i);
        std::string matched = "unmatched";

        std::cout << mzData[i]-0.001 << "\t" << 0 << "\t" << group << "\t" << matched << "\t" << "" << std::endl;
        std::cout << mzData[i] << "\t" << intData[i] << "\t" << group << "\t" << matched << "\t" << "" << std::endl;
        std::cout << mzData[i]+0.001 << "\t" << 0 << "\t" << group << "\t" << matched << "\t" << "" << std::endl;
        ++i;
    }

}

int main(int argc, char * argv[])
{
    std::cout << std::setprecision(8);
    // parse command line arguments
    bool mixed = std::atoi(argv[1]);
    std::string infile = argv[2];
    std::string seq = argv[3];

    // read MS2 file
    std::ifstream spectrumFile(infile);

    // parse mz - intensity pairs
    std::vector<double> mzData, intData;
    double mz, intensity;
    while (spectrumFile >> mz >> intensity)
    {
        //if (intensity > 1000)
        //{
            mzData.push_back(mz);
            intData.push_back(intensity);
        //}
    }

    if (mixed) {
        OpenMS::AASequence pepA = OpenMS::AASequence::fromString("DRVYIHPFHL");
        OpenMS::AASequence pepN = OpenMS::AASequence::fromString("E[-18.010565]LYENKPRRPYIL");

        // create fragments
        std::vector<Fragment> fragments = getFragments(pepA, "A");
        std::vector<Fragment> fragmentsN = getFragments(pepN, "N");

        for (int i = 0; i < fragmentsN.size(); ++i)
        {
            fragments.push_back(fragmentsN[i]);
        }

        std::sort(fragments.begin(), fragments.end());

        matchAndWrite(fragments, mzData, intData);


    } else {

        // get peptide sequence
        OpenMS::AASequence pep;
        if (seq == "A") {
            pep = OpenMS::AASequence::fromString("DRVYIHPFHL");
        } else {
            //pep = OpenMS::AASequence::fromString("E[-18.010565]LYENKPRRPYIL");
            pep = OpenMS::AASequence::fromString(seq);
        }

        std::vector<Fragment> fragments = getFragments(pep, seq);
        std::vector<Fragment> fragments2 = getFragments(OpenMS::AASequence::fromString("VMLMASPSMEDLYHK"), "VMLMASPSMEDLYHK");
        std::vector<Fragment> fragments3 = getFragments(OpenMS::AASequence::fromString("AVFVPDIYSR"), "AVFVPDIYSR");

        for (int i = 0; i < fragments2.size(); ++i)
        {
            fragments.push_back(fragments2[i]);
        }

        for (int i = 0; i < fragments3.size(); ++i)
        {
            fragments.push_back(fragments3[i]);
        }

        std::sort(fragments.begin(), fragments.end());

        matchAndWrite(fragments, mzData, intData);

    }






    return 0;
}