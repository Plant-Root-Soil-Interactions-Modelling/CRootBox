// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGAN_H_
#define ORGAN_H_

#include "mymath.h"

#include <limits>
#include <vector>

namespace CRootBox {

class OrganParameter;
class OrganTypeParameter;
class Organism;

/**
 * Organ
 *
 * Describes a plant organ. Acts as base class for seed, root, stem and leaf.
 *
 * Manages:
 * Organ tree (one parent, multiple children, one plant organism)
 * Organ parameters
 * Organ representation (by nodes, and line segments)
 * Supports RSML output (TODO)
 */
class Organ
{
public:

    Organ(Organism* plant, Organ* parent, int subtype, double delay);
    virtual ~Organ();

    virtual int organType() const; ///< returns the organs type, overwrite for each organ

    /* simulation */
    virtual void simulate(double dt, bool verbose = false); ///< grow for a time span of \param dt

    /* parameters */
    int getId() const { return id; } ///< unique organ id
    OrganParameter* getParam() { return param_; } ///< organ parameters
    OrganTypeParameter* getOrganTypeParameter() const;  ///< organ type parameter
    bool isAlive() { return alive; }
    bool isActive() { return active; }
    double getAge() { return age; }
    double getLength() { return length; }
    void setParent(Organ* p) { parent = p; }
    Organ* getParent() { return parent; }
    virtual double getParameter(std::string name) const; ///< returns an organ parameter

    /* geometry */
    size_t getNumberOfNodes() const { return nodes.size(); } ///< number of nodes of the organ
    Vector3d getNode(int i) const { return nodes.at(i); } ///< i-th node of the organ
    int getNodeId(int i) const { return nodeIds.at(i); } ///< global node index of the i-th node
    double getNodeCT(int i) const { return nodeCTs.at(i); } ///< creation time of the i-th node
    void addNode(Vector3d n, double t); //< adds a node to the root
    void addNode(Vector3d n, int id, double t); //< adds a node to the root
    virtual std::vector<Vector2i> getSegments(int otype = -1) const; // eventually move to Organ1D

    /* for post processing */
    std::vector<Organ*> getOrgans(int otype = -1); ///< the organ including children in a sequential vector
    void getOrgans(int otype, std::vector<Organ*>& v); ///< the organ including children in a sequential vector

    /* IO */
    virtual std::string toString() const; ///< info for debugging todo make pretty
    virtual void writeRSML(std::ostream & cout, std::string indent = "") const; ///< writes a RSML tag todo largely extend

protected:

    /* up and down the organ tree */
    Organism* plant; ///< the plant of which this organ is part of
    Organ* parent; ///< pointer to the parent organ (equals nullptr if it has no parent)
    std::vector<Organ*> children; ///< the successive organs

    /* Parameters that are constant over the organ life time*/
    int id; ///< unique organ id
    OrganParameter* param_; ///< the parameter set of this organ

    /* Parameters are changing over time */
    bool alive = true; ///< true: alive, false: dead
    bool active = true; ///< true: active, false: organ stopped growing
    double age = 0; ///< current age [days]
    double length = 0; ///< length of the organ [cm]

    /* node data */
    std::vector<Vector3d> nodes; ///< nodes of the organ [cm]
    std::vector<int> nodeIds; ///< unique node index
    std::vector<double> nodeCTs; ///< node creation times [days]

};

} // namespace CRootBox

#endif
