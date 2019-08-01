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
 * Collect parameters from the organs
 * Can collect information about the last time step
 * Supports RSML
 * Holds global node index and organ index counter
 * Holds random numbers generator for the organ classes
 */
class Organism {

public:

    enum OrganTypes { ot_organ = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, ot_leaf = 4 }; ///< coarse organ classification
    static std::vector<std::string> organTypeNames; ///< names of the organ types
    static int organTypeNumber(std::string name); ///< organ type number from a string
    static std::string organTypeName(int ot); ///< organ type name from an organ type number

    Organism() { }; ///< empty constructor
    Organism(const Organism& o); ///< copy constructor
    virtual ~Organism(); ///< destructor

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
    virtual std::vector<double> getParameter(std::string name, int ot = -1, std::vector<Organ*> organs = std::vector<Organ*>(0)) const; ///< parameter value per organ
    double getSummed(std::string name, int ot = -1) const; ///< summed up parameters

    /* geometry */
    int getNumberOfOrgans() const { return organId+1; } ///< number of nodes of the organism
    int getNumberOfNodes() const { return nodeId+1; } ///< number of nodes of the organism
    virtual int getNumberOfSegments(int ot=-1) const; ///< number of segments of the organism
    std::vector<std::vector<Vector3d>> getPolylines(int ot=-1) const; ///< nodes per organ
    std::vector<std::vector<double>> getPolylineCTs(int ot=-1) const; ///< node creation times per organ
    virtual std::vector<Vector3d> getNodes() const; ///< nodes of the organ
    virtual std::vector<double> getNodeCTs() const; ///< node creation times, corresponding to Organism::getNodes
    virtual std::vector<Vector2i> getSegments(int ot=-1) const; ///< line segment containing two node indices, corresponding to Organism::getNodes
    virtual std::vector<double> getSegmentCTs(int ot=-1) const; ///< line creation times, corresponding to Organism::getSegments
    virtual std::vector<Organ*> getSegmentOrigins(int ot=-1) const; ///< Points to the organ which contains the segment, corresponding to Organism::getSegments

    /* last time step */
    int getNumberOfNewNodes() const { return getNumberOfNodes()- oldNumberOfNodes; } ///< The number of new nodes created in the previous time step (ame number as new segments)
    int getNumberOfNewOrgans() const { return getNumberOfOrgans() - oldNumberOfOrgans; }  ///< The number of new roots created in the previous time step
    std::vector<int> getUpdatedNodeIndices() const; ///< Indices of nodes that were updated in the previous time step
    std::vector<Vector3d> getUpdatedNodes() const; ///< New coordinates of the updated nodes
    std::vector<Vector3d> getNewNodes() const; ///< Nodes created in the previous time step
    std::vector<Vector2i> getNewSegments(int ot=-1) const; ///< Segments created in the previous time step
    std::vector<Organ*> getNewSegmentOrigins(int ot=-1) const; ///< Copies a pointer to the root containing the new segments
    std::vector<double> getNewSegmentCTs(int ot=-1) const;  ///< Node emergence times created, created in the previous time step

    /* io */
    virtual std::string toString() const; ///< Quick info for debugging
    void readParameters(std::string name, std::string  basetag = "organism"); ///< reads all organ type parameters from a xml file
    void writeParameters(std::string name, std::string basetag = "organism", bool comments = true) const; ///< write all organ type parameters into a xml file
    virtual void writeRSML(std::string name) const; ///< writes a RSML file
    int getRSMLSkip() const { return rsmlSkip; } ///< skips points in the RSML output (default = 0)
    void setRSMLSkip(int skip) { assert(rsmlSkip>=0 && "rsmlSkip must be >= 0" ); rsmlSkip = skip;  } ///< skips points in the RSML output (default = 0)
    std::vector<std::string>& getRSMLProperties() { return rsmlProperties; } ///< reference to the vector<string> of RSML property names, default is { "organType", "subType","length", "age"  }

    /* id management */
    int getOrganIndex() { organId++; return organId; } ///< returns next unique organ id, only organ constructors should call this
    int getNodeIndex() { nodeId++; return nodeId; } ///< returns next unique node id, only organ constructors should call this

    /* random number generator */
    virtual void setSeed(unsigned int seed); ///< Sets the seed of the organisms random number generator
    virtual double rand() { return UD(gen); } ///< Uniformly distributed random number (0,1)
    virtual double randn() { return ND(gen); } ///< Normally distributed random number (0,1)

protected:

    virtual tinyxml2:: XMLElement* getRSMLMetadata(tinyxml2::XMLDocument& doc) const;
    virtual tinyxml2:: XMLElement* getRSMLScene(tinyxml2::XMLDocument& doc) const;

    std::vector<Organ*> baseOrgans;  ///< base organs of the root system

    static const int numberOfOrganTypes = 5;
    std::array<std::map<int, OrganTypeParameter*>, numberOfOrganTypes> organParam;

    double simtime = 0;
    int organId = -1;
    int nodeId = -1;
    int oldNumberOfNodes = 0;
    int oldNumberOfOrgans = 0;

    std::vector<std::string> rsmlProperties = { "organType", "subType","length", "age"  };
    int rsmlSkip = 0; // skips points

    std::mt19937 gen;
    std::uniform_real_distribution<double> UD;
    std::normal_distribution<double> ND;

};

} // namespace

#endif
