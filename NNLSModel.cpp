//
// Created by Dennis Goldfarb on 12/31/17.
//

#include <OpenMS/CONCEPT/Constants.h>
//#include <OpenMS/CHEMISTRY/IsotopeSplineDB.h>
#include <iostream>
#include "NNLSModel.h"

//static const OpenMS::IsotopeSplineDB* isotopeDB = OpenMS::IsotopeSplineDB::getInstance();

NNLSModel::NNLSModel(MultiplexedScan scan, PrecursorTargetOptions options, double massTolerance,
                     MassToleranceUnit unit) {
    this->scan = scan;
    //this->options.push_back(options);
    initOptions(options);
    this->massTolerance = massTolerance;
    this->unit = unit;
    init_b();
    init_A();
}

NNLSModel::NNLSModel(MultiplexedScan scan, std::vector<PrecursorTargetOptions> options, double massTolerance,
                     MassToleranceUnit unit) {
    this->scan = scan;
    this->options = options;
    //initOptions(options);
    this->massTolerance = massTolerance;
    this->unit = unit;
    init_b();
    init_A();
}

void NNLSModel::initOptions(PrecursorTargetOptions options) {
    double totalAbundance = 0;
    for (auto option : options.key2option) {
        //totalAbundance += option.second.abundance;
        totalAbundance = std::max(option.second.abundance, totalAbundance);
    }

    //std::cout << totalAbundance << std::endl;
    PrecursorTargetOptions acceptedOptions;

    for (auto option : options.key2option) {
        if (option.second.abundance / totalAbundance >= 0.2) {
            //std::cout << "accepted: " << option.second.abundance << " " << option.first.z << " " << option.first.firstIso << " " << option.first.lastIso << std::endl;
            acceptedOptions.addOption(option.second);
        }
    }

    if (acceptedOptions.key2option.size() == 0) {
        this->options.push_back(options);
    } else {
        this->options.push_back(acceptedOptions);
    }
}


NNLSModel::NNLSModel(OpenMS::MSSpectrum scan, PrecursorTargetOptions options, double massTolerance,
                     MassToleranceUnit unit) {
    this->scan = MultiplexedScan(scan);
    //initOptions(options);
    this->options.push_back(options);
    this->massTolerance = massTolerance;
    this->unit = unit;
    init_b();
    std::cout << "b init" << std::endl;
    init_A();
    std::cout << "A init" << std::endl;
}

void NNLSModel::init_b() {
    int maxIsotope = getMaxIsotope();
    int maxCharge = getMaxCharge();
    std::cout << "MaxIso: " << maxIsotope << " MaxCharge: " << maxCharge << std::endl;

    auto itrB = b.begin();
    auto itrMz = mzValues.begin();
    // Create vector b and keep track of mz values of each index
    for (int i = 0; i < scan.mzData.size(); ++i) {
        double currentMz = scan.mzData[i];
        minIntensity = std::min(minIntensity, scan.intData[i]);

        // create all isotope m/z signals we will test against
        std::set<double> isotopeMzs;
        for (int z = 1; z <= maxCharge; z++) {
            for (int j = -maxIsotope; j <= maxIsotope; j++) {
                isotopeMzs.insert(currentMz + ((j * OpenMS::Constants::C13C12_MASSDIFF_U) / z));
            }
        }
        std::vector<double> isotopes;
        for (auto isotope : isotopeMzs) {
            isotopes.push_back(isotope);
        }
        std::sort(isotopes.begin(), isotopes.end());

        int k = 0; //i + 1;
        int l = 0;
        bool foundIsotope = false;
        while (k < scan.mzData.size() && l < isotopes.size()) {
            if (Util::compareWithTol(scan.mzData[k], isotopes[l], massTolerance, unit) == 0) {
                foundIsotope = true;
                break;
            } else if (scan.mzData[k] > isotopes[l]) {
                l++;
            } else {
                k++;
            }
        }

        if (!foundIsotope) {
            continue;
        }

        // check if any of the left isotopes are already in b, if not, add them.
        int centerIsotopeIndex = (isotopes.size() / 2);
        int isotopeIndex = centerIsotopeIndex - 1;
        int indexB = itrB - b.begin();

        while (b.size() > 0 && indexB >= 0
               && Util::compareWithTol(mzValues[indexB], isotopes[0], massTolerance, unit) != -1
               && isotopeIndex >= 0) {
            double tmpMz = mzValues[indexB];
            double isotopeMz = isotopes[isotopeIndex];

            int comparison = Util::compareWithTol(isotopeMz, tmpMz, massTolerance, unit);

            if (comparison == 0) // no need to pad. Already in vector
            {
                isotopeIndex--;
                itrB--;
                itrMz--;
            } else if (comparison == -1) // need to keep looking
            {
                itrB--;
                itrMz--;
            } else // comparison == 1 // need to pad
            {
                // pad
                ++itrB;
                ++itrMz;
                itrMz = mzValues.insert(itrMz, isotopeMz);
                itrB = b.insert(itrB, 0.0);
                isotopeIndex--;
            }

            indexB = itrB - b.begin();
        }

        for (; isotopeIndex >= 0; isotopeIndex--) {
            // pad
            if (itrMz != b.begin() && *itrMz < isotopes[0]) {
                while (*itrMz < isotopes[0] && itrMz != mzValues.end()) {
                    ++itrB;
                    ++itrMz;
                }
            }
            itrMz = mzValues.insert(itrMz, isotopes[isotopeIndex]);
            itrB = b.insert(itrB, 0.0);
        }

        // move back to center
        while (itrMz != mzValues.end() && *itrMz < isotopes[centerIsotopeIndex]) {
            ++itrB;
            ++itrMz;
        }

        // add center
        if (itrMz == mzValues.end() || *itrMz != currentMz) {
            itrMz = mzValues.insert(itrMz, currentMz);
            itrB = b.insert(itrB, scan.intData[i]);
        }

        // check if any of the right isotopes are already in b
        // if not, check if they were observed, if yes add obs, if not pad 0s
        isotopeIndex = centerIsotopeIndex + 1;
        indexB = itrB - b.begin();
        while (itrB != b.end() && isotopeIndex < isotopes.size()) {
            double tmpMz = mzValues[indexB];
            double isotopeMz = isotopes[isotopeIndex];

            int comparison = Util::compareWithTol(isotopeMz, tmpMz, massTolerance, unit);

            if (comparison == 0) // no need to pad. Already in vector
            {
                isotopeIndex++;
                itrB++;
                itrMz++;
            } else if (comparison == 1) // need to keep looking
            {
                itrB++;
                itrMz++;
            } else // comparison == -1 // need to pad
            {
                // pad
                // check if signal was observed
                for (int j = i + 1; j < scan.mzData.size(); j++) {
                    int comparison = Util::compareWithTol(isotopeMz, scan.mzData[j], massTolerance, unit);

                    if (comparison == 0) // observed, add it and stop
                    {
                        itrMz = mzValues.insert(itrMz, scan.mzData[j]);
                        itrB = b.insert(itrB, scan.intData[j]);
                        break;
                    } else if (comparison == -1) // too far. Pad and stop
                    {
                        itrMz = mzValues.insert(itrMz, isotopeMz);
                        itrB = b.insert(itrB, 0.0);
                        break;
                    }
                }
                isotopeIndex++;
            }

            indexB = itrB - b.begin();
        }

        for (; isotopeIndex < isotopes.size(); isotopeIndex++) {
            double isotopeMz = isotopes[isotopeIndex];
            // pad
            // check if signal was observed
            for (int j = i + 1; j < scan.mzData.size(); j++) {
                double obsMz = scan.mzData[j];
                int comparison = Util::compareWithTol(isotopeMz, obsMz, massTolerance, unit);

                if (comparison == 0) // observed, add it
                {
                    itrMz = mzValues.insert(itrMz, obsMz);
                    itrB = b.insert(itrB, scan.intData[j]);
                    break;
                } else if (comparison == -1) // too far. Pad and stop
                {
                    itrMz = mzValues.insert(itrMz, isotopeMz);
                    itrB = b.insert(itrB, 0.0);
                    break;
                }
            }
            ++itrB;
            ++itrMz;
        }

        // move to the end
        itrMz = mzValues.end();
        itrB = b.end();
        --itrMz;
        --itrB;
    }
    std::cout << "minIntensity: " << minIntensity << std::endl;
    std::cout << b.size() << std::endl;
}

int NNLSModel::getMaxIsotope() {
    int maxIsotope = 0;
    for (auto option : options) {
        maxIsotope = std::max(maxIsotope, option.getMaxIsotope());
    }
    return maxIsotope;
}

int NNLSModel::getMaxCharge() {
    int maxCharge = 0;
    for (auto option : options) {
        maxCharge = std::max(maxCharge, option.getMaxCharge());
    }
    return maxCharge;
}

void NNLSModel::init_A() {

    numCol = 0;

    for (int index_options = 0; index_options < options.size(); ++index_options) {
        for (auto option : options[index_options].key2option) {
            double monoPrecursorMass = (option.second.maxMass + option.second.minMass) / 2;

            for (int index_b = 0; index_b < b.size(); ++index_b) {
                if (b[index_b] > 0) // this means we observed a signal
                {
                    double mz = mzValues[index_b];
                    double intensity = b[index_b];

                    for (int fragment_z = 1; fragment_z <= option.first.z; ++fragment_z) {
                        if (fragment_z > 1 && option.first.lastIso ==
                                              0) // no point in trying all fragment charges when only mono-isotope expected.
                        {
                            break;
                        }

                        double isotopeStep = OpenMS::Constants::C13C12_MASSDIFF_U / fragment_z;

                        for (int offset = -option.first.lastIso; offset <= 0; ++offset) {
                            double monoFragMz = mz + (offset * isotopeStep);

                            // find closest index of monoFragMz
                            int tmp_index_b = index_b;
                            while (mzValues[tmp_index_b] > monoFragMz) {
                                tmp_index_b--;
                            }

                            // check if the monoFragMz signal was observed. If yes and offset, skip this offset because it was already done
                            if (offset < 0 &&
                                Util::compareWithTol(mzValues[tmp_index_b], monoFragMz, massTolerance, unit) == 0 &&
                                b[tmp_index_b] > 0) {
                                continue;
                            }

                            double monoFragMass = monoFragMz * fragment_z;

                            if (fragment_z == option.first.z && monoFragMass < monoPrecursorMass * 0.66) {
                                continue;
                            }

                            if (monoFragMass >=
                                monoPrecursorMass + 1 + (option.first.z * OpenMS::Constants::PROTON_MASS_U)) {
                                continue;
                            }
                            if (monoFragMass >= monoPrecursorMass) {
                                monoFragMass = monoPrecursorMass;
                            }
                            // approximate isotopic distribution
                            OpenMS::CoarseIsotopePatternGenerator gen;
                            OpenMS::IsotopeDistribution id;
                            id = gen.estimateForFragmentFromPeptideWeight(monoPrecursorMass, monoFragMass,
                                                                          option.second.precursorIsotopes);
                            //id = isotopeDB->estimateForFragmentFromPeptideWeight(monoPrecursorMass, monoFragMass, option.second.precursorIsotopes);
                            id.renormalize();

                            bool validOffset = true;
                            for (int isoCheck = 0; isoCheck < -offset; ++isoCheck) {
                                if (id.getContainer()[isoCheck].getIntensity() >=
                                    id.getContainer()[-offset].getIntensity()) {
                                    validOffset = false;
                                }
                            }

                            /*if (offset <0 && intensity * id.getContainer()[0].second > minIntensity * 2)
                              {
                            validOffset = false;
                            }*/

                            if (!validOffset) {
                                continue;
                            }

                            /*double basePeak = 0.0;
                            double minPeak = 1.0;
                            for (auto peak : id.getContainer())
                            {
                                basePeak = std::max(basePeak, peak.second);
                                minPeak = std::min(minPeak, peak.second);
                            }

                            if (minPeak == 0.0) continue;*/

                            // create dictionary column
                            std::vector<double> intData(id.size());
                            std::vector<double> mzData(id.size());

                            // find the matching index for each isotope
                            for (int i = 0; i < id.size(); ++i) {
                                double isotopeMz =
                                        monoFragMz + (i * (OpenMS::Constants::C13C12_MASSDIFF_U / fragment_z));

                                if (Util::compareWithTol(mzValues[mzValues.size() - 1], isotopeMz, massTolerance,
                                                         unit) != -1) // make sure it's in the scan range
                                {
                                    while (Util::compareWithTol(mzValues[tmp_index_b], isotopeMz, massTolerance,
                                                                unit) != 0) {
                                        tmp_index_b++;
                                    }

                                    //intData[i] = (id.getContainer()[i].second / basePeak) * (b[index_b] / (id.getContainer()[-offset].second / basePeak));
                                    intData[i] = id.getContainer()[i].getIntensity();
                                    if (i == 0 && intData[i] == 0) {
                                        intData[i] = 0.00000001;
                                    }
                                    mzData[i] = mzValues[tmp_index_b];
                                }
                            }

                            // create dictionary element
                            DictionaryElement column(option.second, fragment_z, offset, monoFragMz, intData, mzData,
                                                     numCol);


                            A.push_back(column);
                            numCol++;
                        }
                    }
                }
            }
        }
    }
    std::cout << numCol << " " << b.size() << std::endl;

}

void NNLSModel::writeModel(std::string path, std::string expName) {
    writeMatrixA(path + "A_" + expName + ".bin");
    writeVectorB(path + "b_" + expName + ".bin");
    writePrecursorOptionIndices(path + "indices_" + expName + ".bin");
    //writePrecursorOptionGroupWeights(path + "groupWeights_" + expName + ".bin");
    //writePrecursorOptionIndividualWeights(path + "individualWeights_" + expName + ".bin");
    writeMZs(path + "mz_" + expName + ".tab");
    writePrecursorOptions(path + "precursorOptions_" + expName + ".tab");
}

void NNLSModel::writeMatrixA(std::string path) {
    std::ofstream out(path);

    double missing_value = 0;
    for (DictionaryElement &e : A) {
        int mzIndex = 0;
        for (double mz : mzValues) {
            if (mzIndex == e.mzData.size() || mz < e.mzData[mzIndex]) {
                out.write(reinterpret_cast<char *>(&missing_value), sizeof(missing_value));
            } else {
                out.write(reinterpret_cast<char *>(&e.intData[mzIndex]), sizeof(e.intData[mzIndex]));
                mzIndex++;
            }
        }

    }
    out.close();
}

void NNLSModel::writeVectorB(std::string path) {
    std::ofstream out(path, std::ofstream::binary);

    for (double v : b) {
        out.write(reinterpret_cast<char *>(&v), sizeof(v));
    }

    out.close();
}

void NNLSModel::writePrecursorOptionIndices(std::string path) {
    std::ofstream out(path);

    if (A.size() > 0) {
        int i = 0;
        out.write(reinterpret_cast<char *>(&i), sizeof(i));

        for (; i < A.size() - 1; ++i) {
            if (A[i].precursorDetails != A[i + 1].precursorDetails) {
                out.write(reinterpret_cast<char *>(&i), sizeof(i));
                int next_index = i + 1;
                out.write(reinterpret_cast<char *>(&next_index), sizeof(next_index));
            }
        }

        out.write(reinterpret_cast<char *>(&i), sizeof(i));
    }

    out.close();
}

void NNLSModel::writePrecursorOptionGroupWeights(std::string path) {
    std::ofstream out(path);

    if (A.size() > 0) {
        double v = A[0].precursorDetails.likelihood * sqrt((A[0].precursorDetails.maxIso + 1));
        out.write(reinterpret_cast<char *>(&v), sizeof(v));

        for (int i = 1; i < A.size(); ++i) {
            if (A[i - 1].precursorDetails != A[i].precursorDetails) {
                v = A[i].precursorDetails.likelihood * sqrt((A[i].precursorDetails.maxIso + 1));
                out.write(reinterpret_cast<char *>(&v), sizeof(v));
            }
        }
    }

    out.close();
}

void NNLSModel::writePrecursorOptionIndividualWeights(std::string path) {
    std::ofstream out(path);

    for (int i = 0; i < A.size(); ++i) {
        int num_isotopes = A[i].precursorDetails.maxIso + 1;
        out.write(reinterpret_cast<char *>(&num_isotopes), sizeof(num_isotopes));
    }

    out.close();
}

void NNLSModel::writeMZs(std::string path) {
    std::ofstream out(path);

    for (int i = 0; i < mzValues.size(); ++i) {
        out << mzValues[i] << std::endl;
    }

    out.close();
}

void NNLSModel::writePrecursorOptions(std::string path) {
    std::ofstream out(path);

    if (A.size() > 0) {
        int i = 0;
        PrecursorTargetOption o = A[i].precursorDetails;
        out << o.minMass << "\t" << o.maxMass << "\t" << o.minIso << "\t" << o.maxIso << "\t" << o.charge << "\t"
            << o.likelihood << "\t" << o.getMinMz() << "\t" << o.abundance << std::endl;

        for (; i < A.size() - 1; ++i) {
            if (A[i].precursorDetails != A[i + 1].precursorDetails) {
                o = A[i + 1].precursorDetails;
                out << o.minMass << "\t" << o.maxMass << "\t" << o.minIso << "\t" << o.maxIso << "\t" << o.charge
                    << "\t" << o.likelihood << "\t" << o.getMinMz() << "\t" << o.abundance << std::endl;
            }
        }
    }

    out.close();
}




