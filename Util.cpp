//
// Created by Dennis Goldfarb on 1/2/18.
//

#include "Util.h"

double Util::getTol(double mz, double tol, MassToleranceUnit unit) {
    if (unit == MassToleranceUnit::PPM)
    {
        return OpenMS::Math::ppmToMass(tol, mz);
    }
    else
    {
        return tol;
    }
}

bool Util::withinTol(double mz1, double mz2, double tol, MassToleranceUnit unit) {
    double massTol = getTol(mz1, tol, unit);
    return (mz2 >= mz1 - massTol) && (mz2 <= mz1 + massTol);
}

int Util::compareWithTol(double expected, double observed, double tol, MassToleranceUnit unit) {
    double massTol = getTol(expected, tol, unit);
    double low = expected - massTol;
    double high = expected + massTol;
    if ((observed >= low) && (observed <= high))
    {
        return 0;
    }
    else if (observed > high)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

void Util::populateOptionsForIsolationWindow(int minCharge, int maxCharge, int maxIsotope,
                                       double *charge2prob, double isolationWidth,
                                       double isolationCenter, PrecursorTargetOptions &options)
{

    double isolationRight = isolationCenter + isolationWidth/2;
    double isolationLeft = isolationCenter - isolationWidth/2;
    for (int z = minCharge; z <= maxCharge; ++z)
    {
        double isotopeSpacing = OpenMS::Constants::C13C12_MASSDIFF_U / z;

        int firstIsotope = 0, lastIsotope = 0;

        while (firstIsotope <= maxIsotope)
        {
            // where are the first and last isolated isotopes in m/z space relative to the mono-isotope
            double innerLeftBoundary = firstIsotope * isotopeSpacing;
            double innerRightBoundary = lastIsotope * isotopeSpacing;
            double span = innerRightBoundary - innerLeftBoundary;

            // where are the next isotopes in m/z space relative to the monoisotope
            double outerLeftBoundary = firstIsotope == 0 ? std::numeric_limits<double>::lowest() : (firstIsotope - 1) * isotopeSpacing;
            double outerRightBoundary = lastIsotope == maxIsotope ? std::numeric_limits<double>::max() : (lastIsotope + 1) * isotopeSpacing;

            // How far can the isotopes be shifted (take into account the other side's constraints)
            //double maxRight = isolationWidth > outerRightBoundary - innerLeftBoundary ? outerRightBoundary - isolationWidth : innerLeftBoundary;
            //double maxLeft = isolationWidth > innerRightBoundary - outerLeftBoundary ? outerLeftBoundary : innerRightBoundary - isolationWidth;

            // What is the position of innerLeftBoundary when the isotopes are shifted to the right as far as possible?
            // until the outerLeft hits the left isolation boundary, or the innerRight hits the right isolation boundary, whichever comes first;
            double maxRight;
            // The outerLeft hits the left isolation boundary
            if (isolationWidth >= innerRightBoundary - outerLeftBoundary)
            {
                maxRight = isolationLeft + isotopeSpacing; //outerLeftBoundary; // left position of the isotopes
            }
            else // The innerRight hits the right isolation boundary
            {
                maxRight = isolationRight - span;
            }

            // What is the position of innerLeftBoundary when the isotopes are shifted to the left as far as possible?
            // until the InnerLeft hits the left isolation boundary, or the outerRight hits the right isolation boundary, whichever comes first
            double maxLeft;
            // The outerRight hits the right isolation boundary
            if (isolationWidth > outerRightBoundary - innerLeftBoundary)
            {
                maxLeft = isolationRight - isotopeSpacing - span; //outerRightBoundary - span;
            }
            else // The innerLeft hits the left isolation boundary
            {
                maxLeft = isolationLeft; //innerLeftBoundary;
            }

            double possibleIsolationSpan = maxRight - maxLeft;

            // get precursor mass range
            // determine min/max m/z isolated
            double minMonoMz = maxLeft - innerLeftBoundary;
            double maxMonoMz = maxRight - innerLeftBoundary;
            // get min mono neutral mass
            double minMonoMass = (minMonoMz - OpenMS::Constants::PROTON_MASS_U) * z;
            // get max mono neutral mass
            double maxMonoMass = (maxMonoMz - OpenMS::Constants::PROTON_MASS_U) * z;

            // compute likelihood
            double likelihood = possibleIsolationSpan * charge2prob[z];

            PrecursorTargetOption option(minMonoMass, maxMonoMass, z, firstIsotope, lastIsotope, likelihood, 1.0);
            if (!options.hasOption(option))
            {
                options.addOption(option);
            }


            // update sliding window
            // can we add another isotope?
            // if yes, done
            if (lastIsotope < maxIsotope &&
                (lastIsotope + 1 - firstIsotope) * isotopeSpacing <= isolationWidth)
            {
                lastIsotope++;
            }
                // else, can we subtract the lower isotope?
            else
            {
                firstIsotope++;
                if (firstIsotope > lastIsotope)
                {
                    lastIsotope++;
                }
            }
        }
    }
}

