// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "macropores.h"

namespace CRootBox {

/**
 *  @copydoc Root::createLateral
 *
 *  creates a MacroPoreRoot, instead of a normal one
 */
void Root::createLateral(bool verbose)
{
    int lt = getRootTypeParameter()->getLateralType(nodes.back());
    if (lt>0) {
        double ageLN = this->calcAge(length); // age of root when lateral node is created
        double ageLG = this->calcAge(length+param()->la); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
        double delay = ageLG-ageLN; // time the lateral has to wait
        Root* lateral = new MacroPoreRoot(plant, lt,  heading(), delay,  this, length, nodes.size()-1);
        children.push_back(lateral);
        lateral->simulate(age-ageLN,verbose); // pass time overhead (age we want to achieve minus current age)
    }
}

/**
 *  @copydoc Root::getIncrement
 *
 *  adds macro pore model to the increment
 */
Vector3d MacroPoreRoot::getIncrement(const Vector3d& p, double sdx)
{
    Vector3d sv = Root::getIncrement(p, sdx).normalize();
    if (((MacroPoreRootSystem*)plant)->poreGeometry==nullptr) { // no pores defined
        return sv.times(sdx);
    } else {
        if (((MacroPoreRootSystem*)plant)->poreGeometry->getDist(p)<0) { // inside the pore
            auto sv1 = ((MacroPoreRootSystem*)plant)->applyPoreConductivities(sv);
            // std::cout << "Length before " << sv.length() << ", length after " << sv1.length() << "\n";
            sv1.normalize();
            return sv1.times(sdx);
        } else {
            return sv.times(sdx);
        }
    }
}

}
