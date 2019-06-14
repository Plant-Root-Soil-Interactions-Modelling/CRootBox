// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "OrganParameter.h"

namespace CRootBox {

/**
 *
 */
OrganTypeParameter::OrganTypeParameter(Organism* plant): plant(plant)
{
    iparam["organType"] = &organType;
    iparam["subType"] = &subType;
}

/**
 *
 */
OrganTypeParameter* OrganTypeParameter::copy() const
{
    OrganTypeParameter* o = new OrganTypeParameter(plant); // dparam and iparam must be handled in the constructor
    o->name == name;
    o->organType = organType;
    o->subType = subType;
    return o;
}

/**
 *
 */
OrganParameter* OrganTypeParameter::realize()
{
    OrganParameter* op = new OrganParameter();
    op->subType = subType;
    return op;
}


/**
 * Quick info about the object for debugging
 */
std::string OrganTypeParameter::toString()
{
       std::stringstream str;
       str << "Name " << name << ", organ type "<< organType << ", sub type " << subType;
       return str.str();
}



} // namespace
