// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANPARAMETER_H_
#define ORGANPARAMETER_H_

#include <string>
#include <map>

#include "../external/tinyxml2/tinyxml2.h"

namespace CRootBox {

class Organism; // forward declaration

/**
 * Parameters for a specific organ
 */
class OrganParameter {
public:

    virtual ~OrganParameter() { };

    int subType = -1; ///< sub type of the organ

    virtual std::string toString() const; ///< quick info for debugging

};

/**
 * Contains a parameter set describing as single sub type of an organ,
 * specific parameters are then created with realize().
 *
 * Organizes parameters in hash-maps for scalar double and scalar int values.
 * For this reason derived classes getParameter(), toString(), readXML(), and writeXML() should work out of the box.
 * For other parameter types the methods must be overwritten, see RootTypeParameters.
 *
 * The factory function copy() has to be overwritten for each specialization.
 */
class OrganTypeParameter
{
public:

    OrganTypeParameter(Organism* plant); ///< default constructor
    virtual ~OrganTypeParameter() { };

    virtual OrganTypeParameter* copy(Organism* plant); ///< copies the root type parameter into a new plant

    virtual OrganParameter* realize(); ///< creates a specific organ from the root parameter set

    virtual double getParameter(std::string name) const; // get a scalar parameter

    virtual std::string toString(bool verbose = true) const; ///< info for debugging

    virtual void readXML(tinyxml2::XMLElement* element); ///< reads a single sub type organ parameter set
    void readXML(std::string name); ///< reads a single sub type organ parameter set
    virtual tinyxml2::XMLElement* writeXML(tinyxml2::XMLDocument& doc, bool comments = true) const; ///< writes a organ root parameter set
    void writeXML(std::string name) const; ///< writes a organ root parameter set

    std::string name = "organ";
    int organType = 0;
    int subType = 0;

    Organism* plant;

protected:

    std::map<std::string, double*> dparam; ///< Parameters with type double that can be read and written
    std::map<std::string, int*> iparam; ///< Parameters with type double that can be read and written
    std::map<std::string, double*> param_sd; ///< Deviations of parameters
    std::map<std::string, std::string> description; ///< Parameter descriptions

};

} // namespace

#endif
