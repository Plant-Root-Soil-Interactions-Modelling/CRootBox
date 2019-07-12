// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGANISM_H_
#define ORGANISM_H_

#include "mymath.h"

#include "../external/tinyxml2/tinyxml2.h"

#include <random>
#include <map>
#include <array>

namespace CRootBox {

class Organ;
class OrganTypeParameter;

/**
 * Organism
 *
 * Base class of Plant (CPlantBox) or RootSystem (CRootBox)
 *
 * Manages the OrganTypeParameters
 * Offers an interface for the simulation loop (initialize, simulate, ...)
 * Collects node and line segment geometry from the organ tree
 * Can collect information about the last time step
 * Supports RSML
 * Holds global node index and organ index counter
 * Holds random numbers generator for the organ classes
 */
class Organism {

public:

    enum OrganTypes { ot_organ = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, ot_leaf = 4 }; ///< coarse organ classification

    Organism() { }; ///< empty constructor
    Organism(const Organism& o); ///< copy constructor
    virtual ~Organism();

    /* organ parameter management */
    OrganTypeParameter* getOrganTypeParameter(int otype, int subType) const; ///< returns the respective the type parameter
    std::vector<OrganTypeParameter*> getOrganTypeParameter(int ot) const; ///< returns all type parameters of an organ type (e.g. root)
    void setOrganTypeParameter(OrganTypeParameter* p); ///< sets an organ type parameter, subType and organType defined within p

    /* initialization and simulation */
    void addOrgan(Organ* o) { baseOrgans.push_back(o); } ///< adds an organ, takes ownership
    virtual void initialize(); ///< overwrite for initialization jobs
    virtual void simulate(double dt, bool verbose = false); ///< calls the base organs simulate methods
    double getSimTime() const { return simtime; } ///< returns the current simulation time

    /* organs as sequential list */
    std::vector<Organ*> getOrgans(int ot=-1) const; ///< sequential list of organs
    virtual std::vector<double> getParameters(std::string name, int ot = -1, std::vector<Organ*> organs = std::vector<Organ*>(0)) const; ///< parameter value per organ
    double getSummed(std::string name, int ot = -1) const;

    /* geometry */
    int getNumberOfNodes() const { return nodeId+1; } ///< number of nodes of the organism
    int getNumberOfOrgans() const { return organId+1; } ///< number of nodes of the organism
    std::vector<std::vector<Vector3d>> getPolylines(int ot=-1) const; ///< nodes per organ
    std::vector<std::vector<double>> getPolylinesIds(int ot=-1) const; ///< global node indices per organ
    std::vector<std::vector<double>> getPolylinesCTs(int ot=-1) const; ///< node creation times per organ
    virtual std::vector<Vector3d> getNodes() const; //
    virtual std::vector<double> getNodeCTs() const;
    virtual std::vector<Vector2i> getSegments(int ot=-1) const;
    virtual std::vector<double> getSegmentCTs(int ot=-1) const;
    virtual std::vector<Organ*> getSegmentOrigins(int ot=-1) const; ///< Copies a pointer to the root containing the segment

    /* last time step */


    /* io */
    virtual std::string toString() const;
//    void readParameter(std::string name);
//    void writeParameter(std::string name);
    virtual void writeRSML(std::string name) const; ///< writes a RSML file
    virtual tinyxml2:: XMLElement* getRSMLMetadata(tinyxml2::XMLDocument& doc) const;
    virtual tinyxml2:: XMLElement* getRSMLScene(tinyxml2::XMLDocument& doc) const;
    int getRSMLSkip() { return rsmlSkip; }
    std::vector<std::string> getRSMLProperties() { return rsmlProperties; }

    /* id management */
    int getOrganIndex() { organId++; return organId; } ///< returns next unique organ id, only organ constructors should call this
    int getNodeIndex() { nodeId++; return nodeId; } ///< returns next unique node id, only organ constructors should call this

    /* random number generator */
    void setSeed(unsigned int seed); ///< Sets the seed of the organisms random number generator
    double rand() { return UD(gen); } ///< Uniformly distributed random number (0,1)
    double randn() { return ND(gen); } ///< Normally distributed random number (0,1)

protected:

    std::vector<Organ*> baseOrgans;  ///< base organs of the root system

    static const int numberOfOrganTypes = 5;
    std::array<std::map<int, OrganTypeParameter*>, numberOfOrganTypes> organParam;

    double simtime = 0;
    int organId = -1;
    int nodeId = -1;

    std::vector<std::string> rsmlProperties = { "organType", "subType","length", "age"  };
    int rsmlSkip = 0; // skips points

    std::mt19937 gen;
    std::uniform_real_distribution<double> UD;
    std::normal_distribution<double> ND;

};

} // namespace

#endif
