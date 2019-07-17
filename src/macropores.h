// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef MACROPORES_H_
#define MACROPORES_H_

#include "RootSystem.h"

namespace CRootBox {

/**
 * Include the macro pore model of Landl et al. 2016 into the RootBox model
 *
 * D. Leitner, 2019
 */
class MacroPoreRootSystem : public RootSystem {
public:

    // todo overwrite initialize to create MacroPoreRoot as base organs

    // Macro pores
    void setPoreGeometry(SignedDistanceFunction* geom) { poreGeometry = geom; }
    void setPoreConductivity(Matrix3d m) { poreConductivity = m; }
    void setPoreLocalAxes(Matrix3d m) {
        poreLocalAxes = m;
        invPoreLocalAxes = m.inverse();
    }
    Vector3d applyPoreConductivities(Vector3d v) {
        return poreLocalAxes.times(poreConductivity.times(invPoreLocalAxes.times(v))); // Landl et al. 2016, Eqn (12), K*v = [M*K'*(M^-1)]*v
    }

private:

    // todo currently everything is constant... // (move to specific application)
    Matrix3d poreConductivity = Matrix3d();
    Matrix3d poreLocalAxes= Matrix3d(); ///< Copies a scalar root parameter that is constant per root to a vector
    Matrix3d invPoreLocalAxes = Matrix3d();
    SignedDistanceFunction* poreGeometry = nullptr;

};

/**
 * An extension so that roots grow in a macropore geometry
 */
class MacroPoreRoot :public Root {

    MacroPoreRoot(Organism* rs, int type, Vector3d pheading, double delay, Root* parent, double pbl, int pni) :Root(rs, type, pheading, delay, parent, pbl, pni) { }

protected:

    Vector3d getIncrement(const Vector3d& p, double sdx) override; ///< called by createSegments, to determine growth direction
    void createLateral(bool silence) override; ///< creates a new lateral, called by Root::simulate()

};

} // end namespace CRootBox

#endif
