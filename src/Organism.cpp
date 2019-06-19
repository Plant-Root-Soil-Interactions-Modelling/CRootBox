// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Organism.h"

#include "OrganParameter.h"

#include <stdexcept>
#include <iostream>

namespace CRootBox {

/**
 * Copy constructor
 */
Organism::Organism(const Organism& o): organParam(o.organParam), simtime(o.simtime),
    organId(o.organId), nodeId(o.nodeId), gen(o.gen), UD(o.UD), ND(o.ND)
{ }

/*
 * Destructor: deletes all organ type parameters
 */
Organism::~Organism()
{
    for (int ot = 0; ot < numberOfOrganTypes; ot++) {
        std::vector<OrganTypeParameter*> deleteMe = getOrganTypeParameter(ot);
        for (int i=0; i<deleteMe.size(); i++) {
            delete deleteMe[i];
        }
    }
}

/**
 * Copies the organ type parameters of a specific organ type into a vector
 */
std::vector<OrganTypeParameter*> Organism::getOrganTypeParameter(int otype) const
{
    std::vector<OrganTypeParameter*>  otps = std::vector<OrganTypeParameter*>(0);
    for (auto& otp : organParam[otype]) {
        otps.push_back(otp.second);
    }
    return otps;
}

/**
 * Returns an organ type parameter of a specific organ type and sub type
 */
OrganTypeParameter* Organism::getOrganTypeParameter(int otype, int subtype) const
{
    try {
//        std::cout << "reading organ type " << otype << " sub type " << subtype <<": ";
//        for (auto& p : organParam[otype]) {
//            std::cout << p.first;
//        }
//        std::cout << "\n" << std::flush;;
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
