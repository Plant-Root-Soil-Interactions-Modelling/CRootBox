#ifndef ORGANISM_H_
#define ORGANISM_H_

#include "mymath.h"

#include <random>
#include <map>
#include <array>

namespace CRootBox {

class Organ;
class OrganTypeParameter;

/**
 * Organism
 *
 * Base class of Plant or RootSystem
 *
 * Manages the OrganTypeParameters
 * Offers an interface for the simulation loop (initialize, simulate, getSimTime)
 * Collects node and line segment geometry from the organ tree
 * Holds global node index and organ index counter
 * Holds random numbers generator for the organ classes
 */
class Organism {

public:

    enum OrganTypes { ot_organ = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, ot_leaf = 4 };

    Organism() { };
    Organism(const Organism& o); ///< copy constructor

    virtual ~Organism();

    /* organ parameter management */
    OrganTypeParameter* getOrganTypeParameter(int otype, int subType) const;
    std::vector<OrganTypeParameter*> getOrganTypeParameter(int otype) const;
    void setOrganTypeParameter(OrganTypeParameter* p);

    /* for the simulation loop */
    virtual void initialize() { }
    virtual void simulate(double dt, bool verbose = false) { }
    double getSimTime() const { return simtime; } ///< returns the current simulation time

    /* id management */
    int getOrganIndex() { organId++; return organId; } ///< returns next unique organ id
    int getNodeIndex() { nodeId++; return nodeId; } ///< returns next unique node idvecto

    /* geometry */
    virtual int getNumberOfNodes() const { return nodeId+1; } ///< Number of nodes of the root system (including nodes for seed, root crowns, and artificial shoot)
    virtual std::vector<Vector3d> getNodes() const { return std::vector<Vector3d>(0); }
    virtual std::vector<Vector2i> getSegments(int otype=-1) const { return std::vector<Vector2i>(0); }
    virtual std::vector<double> getSegmentCTs(int otype=-1) const { return std::vector<double>(0); }
    virtual std::vector<Organ*> getSegmentOrigins(int otype=-1) const { return std::vector<Organ*>(0); } ///< Copies a pointer to the root containing the segment

    /**
     * Sets the seed of the organisms random number generator.
     * In order to obtain two exact same organisms call before Organism::initialize().
     *
     * @param seed      random number generator seed
     */
    void setSeed(unsigned int seed) { this->gen = std::mt19937(seed); } ///< Sets the seed of the organisms random number generator
    double rand() { return UD(gen); } ///< Uniformly distributed random number (0,1)
    double randn() { return ND(gen); } ///< Normally distributed random number (0,1)

protected:

    static const int numberOfOrganTypes = 5;
    std::array<std::map<int, OrganTypeParameter*>, numberOfOrganTypes> organParam;

    double simtime = 0;

    int organId = -1;
    int nodeId = -1;

    std::mt19937 gen;
    std::uniform_real_distribution<double> UD;
    std::normal_distribution<double> ND;
};

} // namespace

#endif
