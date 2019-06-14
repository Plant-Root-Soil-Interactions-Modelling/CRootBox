// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Organism.h"

#include "OrganParameter.h"

#include <stdexcept>
#include <iostream>

namespace CRootBox {

Organism::Organism(const Organism& o): organParam(o.organParam), simtime(o.simtime),
    organId(o.organId), nodeId(o.nodeId), gen(o.gen), UD(o.UD), ND(o.ND)
{ }

Organism::~Organism()
{
    // TODO delete all parameter type classes
}


/**
 * Returns the organ type parameter
 */
OrganTypeParameter* Organism::getOrganTypeParameter(int otype, int subtype) const
{
    try {
        std::cout << "reading ogran type " << otype << " sub type " << subtype << "\n";
        return organParam[otype].at(subtype);
    } catch(const std::out_of_range& oor) {
        std::cout << "Organism::getOrganTypeParameter: Organ type parameter of sub type " << subtype << " wast not set \n" << std::flush;
        throw;
    }
}

/**
 *  Sets the parameter type
 *  Deletes the old parameter if there is one, takes ownership of the new one
 */
void Organism::setOrganTypeParameter(OrganTypeParameter* p)
{
    int otype = p->organType;
    int subtype = p->subType;
    try {
        delete organParam[otype][subtype];
    } catch (std::exception& e) {
        // did not exist, nothing to delete
    }
    std::cout << "setting ogran type " << otype << " sub type " << subtype << "\n";
    organParam[otype][subtype] = p;
}

}
