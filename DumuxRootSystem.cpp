#include "DumuxRootSystem.h"

#ifdef HAVE_DUNE


DumuxRootSystem::Grid* DumuxRootSystem::getGrid(Grid& grid)
{
    try {
    Dune::GridFactory<Grid> factory;           ;
    // fetch stuff
    auto roots = getRoots();
    auto nodes = getNodes(RootSystem::ot_segments,roots);
    // copy vertices
    for (const auto& n: nodes) {
        FVector3 pos({n.x,n.y,n.z});
        factory.insertVertex(pos);
    }
    // create cells
    Dune::GeometryType gt(Dune::GeometryType::simplex,1);
    std::vector<unsigned int> line(2);
    for (auto const& r:roots) {
        for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
            line[0] = r->getNodeId(i);
            line[1] = r->getNodeId(i+1);
            factory.insertElement(gt,line);
        }
    }
    // create grid
    Grid* grid_ = factory.createGrid();
    return grid_;
    } catch (Dune::Exception & e) {
        return nullptr;
    }

}



/**
 * Call this first to set DuMux indices for the initial nodes
 * (for tap root system this is exactly one node)
 */
void DumuxRootSystem::getInitialNodes()
{
    clearLastTimeStep();
    for (auto const& b: baseRoots) {
        newnodes.push_back(b->getNode(0));
        rbids.push_back(b->getNodeId(0));
        dumuxids.push_back(-1); // -1 should indicate there is no predecessor
    }
}

/**
 * Find new nodes created between t0 and t1, and create the corresponding members
 * does not find new basal roots (TODO)
 *
 * @param t0    start time [day]
 * @param t1    end time [day]
 */
void DumuxRootSystem::findNewNodes(double t0, double t1)
{
    clearLastTimeStep();

    std::vector<Root*> roots = getRoots();
    for (auto const& r:roots) {
        for (size_t i=1; i<r->getNumberOfNodes(); i++) {

            // cout << r->getNodeId(i) << " " << r->getNodeId(i+1) << " ";
            double t = r->getNodeETime(i);

            if ((t>=t0) && (t<t1)) {

                newnodes.push_back(r->getNode(i));
                rbids.push_back(r->getNodeId(i)); // node with rootbox index i
                int rbi = r->getNodeId(i-1); // root box index of the predecessor
                dumuxids.push_back(idmap.at(rbi));  // look up dumux index

            }

        }
    }

}

/**
 * Stores the DuMux identifier or the newly created nodes
 *
 * @param dids      the node indices from dumux corresponing to newnodes
 */
void DumuxRootSystem::setDumuxIds(std::vector<int> dids)
{
    for(size_t i = 0; i<dids.size();  i++) {
        //cout << "       idmap["<<rbids.at(i)<<"] = "<< dids.at(i) << "\n";
        idmap.at(rbids.at(i)) = dids.at(i);
    }
}


/**
 * Clears newnodes, dumuxids, rbids, and adjust the size of idmap
 */
void DumuxRootSystem::clearLastTimeStep()
{
    int non = getNumberOfNodes();
    newnodes.clear();
    dumuxids.clear();
    rbids.clear();
    while (int(idmap.size())<non) { // append zeros to achieve right length
        idmap.push_back(0);
    }
}

#endif /* HAVE_DUNE */
