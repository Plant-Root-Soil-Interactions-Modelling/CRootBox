// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOTSYSTEMPARAMETER_H_
#define ROOTSYSTEMPARAMETER_H_

#include "mymath.h"
#include "OrganParameter.h"

namespace CRootBox {

/**
 * RootSystemParameter: contains all plant specific parameters like planting depth and describing the emergence times of basal and shoot borne roots
 *
 * This model is very limited in the moment, and we have to replace it, if we come up with something better
 */
class RootSystemParameter :public OrganParameter
{

public:

    RootSystemParameter(); ///< Default constructor
    virtual ~RootSystemParameter() { };

    virtual void set(double pd, double fB, double dB, int mB, int nC, double fSB, double dSB, double dRC, double nz, double simtime); ///< Sets all the parameters

    virtual void read(std::istream & cin); ///< Read plant parameters
    virtual void write(std::ostream & cout) const; ///< Write plant parameters

    /* Plant parameters */
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
};

} // namespace

#endif
