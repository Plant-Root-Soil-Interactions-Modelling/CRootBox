// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANPARAMETER_H_
#define ORGANPARAMETER_H_

#include <string>
#include <map>
#include <sstream>

namespace CRootBox {

class Organism; // forward declaration

/**
 * Parameter for a specific organ
 */
class OrganParameter {
public:

    OrganParameter() { };

    int subType = -1; // sub type of the organ

};

/**
 * Contains a parameter set describing as single sub type of an organ
 * Specific parameters are created with realize()
 */
class OrganTypeParameter
{
public:

    OrganTypeParameter(Organism* plant); ///< default constructor
    virtual ~OrganTypeParameter() { };

    virtual OrganParameter* realize(); ///< Creates a specific organ from the root parameter set

    virtual std::string toString(); ///< info for debugging

//    virtual void readXML(std::istream & cin); ///< reads a single sub type organ parameter set

//    virtual void writeXML(std::ostream & cout) const; ///< writes a organ root parameter set

    std::string name = "Unnamed organ";
    int organType = 0;
    int subType = 0;

    Organism* plant;

    std::map<std::string, double*> dparam; ///< Parameters with type double that can be read and written
    std::map<std::string, int*> iparam; ///< Parameters with type double that can be read and written

};

} // namespace

#endif
