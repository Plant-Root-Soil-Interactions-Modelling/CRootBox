// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef SEEDPARAMETER_H_
#define SEEDPARAMETER_H_

#include "mymath.h"
#include "organparameter.h"

/**
 * This file describes the classes SeedSpecificParameter and SeedRandomParameter.
 * SeedSpecificParameter are drawn from the SeedRandomParameter class
 */

namespace CRootBox {

/**
 * SeedSpecificParameter contains all plant specific parameters like planting depth and describing
 * the emergence times of basal and shoot borne roots
 *
 * The model currently rather limited in, and we might replace it, if we come up with something better
 */
class SeedSpecificParameter :public OrganSpecificParameter
{

public:
    SeedSpecificParameter():SeedSpecificParameter(Vector3d(0.,0.,-3),1.e9,0.,0,0.,1.e9,1.e9,1.e9,0.,30.) { }; ///< Default constructor
    SeedSpecificParameter(Vector3d seedPos, double fB, double dB, int mB, int nC, double fSB, double dSB,
    		double dRC, double nz, double simtime): seedPos(seedPos), firstB(fB), delayB(dB),
    	    		maxB(mB), nC(nC), firstSB(fSB), delaySB(dSB), delayRC(dRC), nz(nz), simtime(simtime) { };
    virtual ~SeedSpecificParameter() { };

    /*
     * RootBox and PlantBox plant parameters
     */
    Vector3d seedPos;   ///< Position of the seed [cm]
    //Basal roots (nodal roots)
    double firstB;  ///< Emergence of first basal root [day]
    double delayB;  ///< Time delay between the basal roots [day]
    int maxB;       ///< Maximal number of basal roots [1]
    //Shoot borne roots (crown roots)
    int nC;             ///< Maximal number of roots per root crown [1]
    double firstSB;     ///< First emergence of a shoot borne root [day]
    double delaySB;     ///< Time delay between the shoot borne roots [day]
    double delayRC;     ///< Delay between the root crowns [day]
    double nz;          ///< Distance between the root crowns along the shoot [cm]
    //Simulation parameters
    double simtime;    ///< recommended final simulation time

    std::string toString() const override; ///< for debugging
};



/**
 * Contains a parameter set describing a plant
 */
class SeedRandomParameter :public OrganRandomParameter {
public:

	SeedRandomParameter(Organism* plant); ///< default constructor
    virtual ~SeedRandomParameter() { }; ///< nothing to do

    OrganRandomParameter* copy(Organism* plant_) override;

    OrganSpecificParameter* realize() override; ///< Creates a specific plant from the seed random parameter set

    // DEPRICATED
    void read(std::istream & cin); ///< reads a single root system parameter set
    void write(std::ostream & cout) const; ///< writes a single root system parameter set

    /* Plant parameters */
    Vector3d seedPos;   ///< Mean position of the seed [cm]
    Vector3d seedPoss;  ///< Standard deviation of position  [cm]

    //Basal roots (nodal roots)
    double firstB;  ///< Mean emergence of first basal root [day]
    double firstBs; ///< Standard deviation of emergence of first basal root [day]
    double delayB;  ///< Mean time delay between the basal roots [day]
    double delayBs; ///< Standard deviation of time delay between the basal roots [day]
    double maxB;    ///< Mean maximal number of basal roots [1]
    double maxBs;   ///< Standard deviation of maximal number of basal roots [1]

    //Shoot borne roots (crown roots)
    double nC;          ///< Mean maximal number of roots per root crown [1]
    double nCs;         ///< Standard deviation of maximal number of roots per root crown [1]
    double firstSB;     ///< Mean first emergence of a shoot borne root [day]
    double firstSBs;     ///< Standard deviation of first emergence of a shoot borne root [day]
    double delaySB;     ///< Mean time delay between the shoot borne roots [day]
    double delaySBs;    ///< Standard deviation of time delay between the shoot borne roots [day]
    double delayRC;     ///< Mean delay between the root crowns [day]
    double delayRCs;    ///< Standard deviation of delay between the root crowns [day]
    double nz;          ///< Mean distance between the root crowns along the shoot [cm]
    double nzs;         ///< Standard deviation of distance between the root crowns along the shoot [cm]

    //Simulation parameters
    double simtime;     ///< Mean recommended final simulation time
    double simtimes;    ///< Standard deviation of recommended final simulation time

protected:

    void bindParmaters(); ///<sets up class introspection

};

} // namespace

#endif
