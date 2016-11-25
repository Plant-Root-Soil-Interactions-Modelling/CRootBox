#ifndef DUMUXROOTSYSTEM_H_
#define DUMUXROOTSYSTEM_H_

//#define HAVE_DUNE // put comment here, if you do not want to use the dune coupling

#ifdef HAVE_DUNE

#include "RootSystem.h"

#include "config.h"
#include "dune/foamgrid/foamgrid.hh"
#include "dune/grid/yaspgrid.hh"
#include "dune/grid/common/gridfactory.hh"
#include "dune/grid/common/gridinfo.hh"
#include "dune/common/fvector.hh"


/**
 * A RootSystem,
 * with additional getters and setters, for the link to FoamGrid<1,3>
 *
 * TODO TODO TODO
 */
class DumuxRootSystem : public RootSystem {

public:

    //typedef Dune::FoamGrid<1,3> Grid;
    typedef Dune::YaspGrid<3> Grid;
    typedef Grid::ctype CType;
    typedef Dune::FieldVector<CType,3> FVector3;

    Grid* getGrid(Grid& grid); ///< returns or updates a foamgrid (TODO update)

    void getInitialNodes(); ///< Call this first to set DuMux indices for the initial nodes
    void findNewNodes(double t0, double t1); ///< Find the nodes created between t0 and t1, and create the corresponding members
    std::vector<Vector3d> getNewNodes() { return newnodes; } ///< the new nodes that were created between t0 and t1
    std::vector<int> getDumuxIds() { return dumuxids; } ///< the Dumux identifiers of the existing nodes to which the new nodes will be attached

    void setDumuxIds(std::vector<int> dids); ///< set the Dumux identifiers after creating the new nodes

private:

    void clearLastTimeStep(); ///< Clears newnodes, dumuxids, rbids, and adjust the size of idmap

    std::vector<Vector3d> newnodes; ///< set by setNewNodes()
    std::vector<int> dumuxids; ///< corresponding Dumux identifiers of predecessing nodes
    std::vector<int> rbids; ///< Rootbox identifiers for the new nodes

    std::vector<int> idmap; ///< maps rootbox ids to dumux ids, i.e. dumuxid = idmap[rbid]
};


#endif /* HAVE_DUMUX */

#endif /* DUMUXROOTSYSTEM_H_ */
