//
// Created by Dennis Goldfarb on 12/31/17.
//

#ifndef DEMIX_NNLSMODEL_H
#define DEMIX_NNLSMODEL_H

#include <fstream>
#include <vector>
#include <algorithm>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include "OpenMS/MATH/MISC/MathFunctions.h"
#include "OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h"
#include "OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h"

#include "PrecursorTargetOptions.h"
#include "MultiplexedScan.h"
#include "MassToleranceUnit.h"
#include "Util.h"
#include "DictionaryElement.h"

class NNLSModel {

public:

    NNLSModel(MultiplexedScan scan, PrecursorTargetOptions options, double massTolerance, MassToleranceUnit unit);

    NNLSModel(OpenMS::MSSpectrum scan, PrecursorTargetOptions options, double massTolerance, MassToleranceUnit unit);

    NNLSModel(MultiplexedScan scan, std::vector<PrecursorTargetOptions> options, double massTolerance, MassToleranceUnit unit);

    int getMaxIsotope();
    int getMaxCharge();

    double massTolerance;
    MassToleranceUnit unit;
    int numCol;
    double minIntensity = 1e9;

    std::vector<double> b;
    std::vector<double> mzValues;
    std::vector<DictionaryElement> A;

    std::vector<double> b2;
    std::vector<double> mzValues2;

    MultiplexedScan scan;
    std::vector<PrecursorTargetOptions> options;

    void writeModel(std::string path, std::string expName);
    void writeMatrixA(std::string path);
    void writeVectorB(std::string path);
    void writePrecursorOptionIndices(std::string path);
    void writePrecursorOptionGroupWeights(std::string path);
    void writePrecursorOptionIndividualWeights(std::string path);
    void writeMZs(std::string path);
    void writePrecursorOptions(std::string path);

private:
    void init_b();
    void init_A();
    void initOptions(PrecursorTargetOptions options);

    void Bcontains(double mz, int indexB);
};



#endif //DEMIX_NNLSMODEL_H
