// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANPARAMETER_H_
#define ORGANPARAMETER_H_

#include <string>
#include <map>

#include "../external/tinyxml2/tinyxml2.h"

/**
 * This file describes the classes OrganSpecificParameter and RootRandomParameter.
 * OrganSpecificParameter are drawn from the RootRandomParameter class
 */

namespace CRootBox {

class Organism; // forward declaration

/**
 * Parameters for a specific organ
 */
class OrganSpecificParameter {
public:

    virtual ~OrganSpecificParameter() { };

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
class OrganRandomParameter
{
public:

    OrganRandomParameter(Organism* plant); ///< default constructor
    virtual ~OrganRandomParameter() { };

    virtual OrganRandomParameter* copy(Organism* plant); ///< copies the root type parameter into a new plant

    virtual OrganSpecificParameter* realize(); ///< creates a specific organ from the root parameter set

    virtual double getParameter(std::string name) const; // get a scalar parameter

    virtual std::string toString(bool verbose = true) const; ///< info for debugging

    virtual void readXML(tinyxml2::XMLElement* element); ///< reads a single sub type organ parameter set
    void readXML(std::string name); ///< reads a single sub type organ parameter set
    virtual tinyxml2::XMLElement* writeXML(tinyxml2::XMLDocument& doc, bool comments = true) const; ///< writes a organ root parameter set
    void writeXML(std::string name) const; ///< writes a organ root parameter set

    void bindParameter(std::string name, int* i, std::string descr = "", double* dev = nullptr); ///< binds integer to parameter name
    void bindParameter(std::string name, double* d, std::string descr = "", double* dev = nullptr); ///< binds double to parameter name

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
