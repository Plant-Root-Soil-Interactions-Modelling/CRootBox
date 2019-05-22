#ifndef PLANTBASE_H_
#define PLANTBASE_H_

#include "mymath.h"

#include <random>
#include <vector>

namespace CRootBox {

class Organ;
class OrganTypeParameter;

/**
 * Base class of Plant or RootSystem
 */
class PlantBase {

public:

    enum OrganTypes { ot_organ = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, ot_leaf = 4};

    virtual ~PlantBase() { };

    virtual OrganTypeParameter* getOrganTypeParameter(int otype, int subtype) const { return nullptr; }

    /* for the simulation loop */
    virtual void initialize() { };
    virtual void simulate(double dt, bool verbose = false) { };

    /* id management */
    int getOrganIndex() { organId++; return organId; } ///< returns next unique organ id
    int getNodeIndex() { nodeId++; return nodeId; } ///< returns next unique node id

    /* geometry */
    virtual int getNumberOfNodes() const { return nodeId+1; } ///< Number of nodes of the root system (including nodes for seed, root crowns, and artificial shoot)
    virtual std::vector<Vector3d> getNodes() const { return std::vector<Vector3d>(0); };
    virtual std::vector<Vector2i> getSegments(int otype= -1) const { return std::vector<Vector2i>(0); } ;
    virtual std::vector<double> getSegmentCTs() const { return std::vector<double>(0); } ;
    virtual std::vector<Organ*> getSegmentsOrigin(unsigned int otype= -1) const {
        return std::vector<Organ*>(0);
    } ///< Copies a pointer to the root containing the segment

    /* random stuff */
    /**
     * Sets the seed of the root systems random number generator.
     * To obtain two exact same root system call PlantBase::setSeed() before initialize().
     *
     * @param seed      random number generator seed
     */
    void setSeed(unsigned int seed) { this->gen = std::mt19937(seed); }; ///< help fate (sets the seed of all random generators)

    double rand() { return UD(gen); } ///< Uniformly distributed random number (0,1)

    double randn() { return ND(gen); } ///< Normally distributed random number (0,1)

protected:

    int organId = -1;
    int nodeId = -1;

    std::mt19937 gen;
    std::uniform_real_distribution<double> UD;
    std::normal_distribution<double> ND;

};

} // namespace

#endif
