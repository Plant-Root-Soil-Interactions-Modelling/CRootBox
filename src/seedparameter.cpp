// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include <cmath>
#include "seedparameter.h"

namespace CRootBox {

/*
 * class RootSystemParameter
 */

/**
 * Default constructor: No basal roots, not shoot borne roots, planting depth 3 [cm]
 */
SeedSpecificParameter::SeedSpecificParameter() {
    set(3.,1.e9,0.,0, //pd, fB, dB, mB,
        0,1.e9,1.e9,0.,0.,30.);  // nC, fSB, dSB, dRC, nz
}

/**
 * depricated todo no RootSystemTypeParameter yet
 */
void SeedSpecificParameter::read(std::istream & cin) {
    double plantingdepth;
    std::string s; // dummy
    cin  >>  s >> plantingdepth;
    cin >> s >> firstB >> s >> delayB >> s >> maxB >> s >> nC >> s >> firstSB >> s >> delaySB >> s >> delayRC >> s >> nz >> s >> simtime;
    seedPos = Vector3d(0,0,-plantingdepth);
}

/**
 * depricated todo no RootSystemTypeParameter yet
 */
void SeedSpecificParameter::write(std::ostream & cout) const {
    double pd = -seedPos.z;
    cout <<  "plantingdepth\t" << pd << "\n" <<  "firstB\t" << firstB << "\n" <<  "delayB\t" << delayB << "\n"
        <<  "maxB\t" << maxB << "\n" <<  "nC\t" << nC << "\n" <<  "firstSB\t" << firstSB << "\n"
        <<  "delaySB\t" << delaySB << "\n" <<  "delayRC\t" << delayRC << "\n" <<  "nz\t" << nz << "\n" << "simulationTime\t" << simtime << "\n";
}

/**
 * Set all plant parameters
 *
 * @param pd          Planting depth [cm]
 * @param fB          Emergence of first basal root [day]
 * @param dB          Time delay between the basal roots [day]
 * @param mB          Maximal number of basal roots [1]
 * @param nC          Maximal number of roots per root crown [1]
 * @param fSB         First emergence of a shoot borne root [day]
 * @param dSB         Time delay between the shoot borne roots [day]
 * @param dRC         Delay between the root crowns [day]
 * @param nz          Distance between the root crowns along the shoot [cm]
 * @param simtime     Recommended final simulation time (e.g. used in the web interface)
 */
void SeedSpecificParameter::set(double pd, double fB, double dB, int mB, int nC, double fSB, double dSB, double dRC, double nz, double simtime) {
    seedPos=Vector3d(0,0,-pd);
    firstB=fB;
    delayB=dB;
    maxB=mB;
    this->nC=nC;
    firstSB=fSB;
    delaySB=dSB;
    delayRC=dRC;
    this->nz=nz;
    this->simtime=simtime;
}


} // namespace
