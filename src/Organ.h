// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ORGAN_H_
#define ORGAN_H_

#include "mymath.h"

#include "../external/tinyxml2/tinyxml2.h"

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
 * Organ development (@see Organ::simulate)
 * The organ tree: one parent, multiple children, one plant organism
 * This organ's parameters
 * This organ's geometry (by nodes, global node indices, node creation times, and line segments)
 * Information about the last time step: nodes can either move, or be created
 * Post processing and RSML output
 */
class Organ
{
public:

    Organ(int id, const OrganParameter* param, bool alive, bool active, double age, double length,
        bool moved= false, int oldNON = 0); ///< creates everything from scratch
    Organ(Organism* plant, Organ* parent, int organtype, int subtype, double delay); ///< used within simulation
    virtual ~Organ();

    virtual Organ* copy(Organism* plant); ///< deep copies the organ tree

    virtual int organType() const; ///< returns the organs type, overwrite for each organ

    /* development */
    virtual void simulate(double dt, bool verbose = false); ///< grow for a time span of \param dt

    /* tree */
    void setParent(Organ* p) { parent = p; } ///< sets parent organ
    Organ* getParent() const { return parent; } ///< return parent organ, equals nullptr if it has no parent
    void setOrganism(Organism* p) { plant = p; } ///< sets the organism of which the organ is part of
    void addChild(Organ* c); ///< adds an subsequent organ

    /* parameters */
    int getId() const { return id; } ///< unique organ id
    const OrganParameter* getParam() const { return param_; } ///< organ parameters
    OrganTypeParameter* getOrganTypeParameter() const;  ///< organ type parameter
    bool isAlive() const { return alive; } ///< checks if alive
    bool isActive() const { return active; } ///< checks if active
    double getAge() const { return age; } ///< return age of the organ
    double getLength() const { return length; } ///< returns length of the organ

    /* geometry */
    int getNumberOfNodes() const { return nodes.size(); } ///< number of nodes of the organ
    Vector3d getNode(int i) const { return nodes.at(i); } ///< i-th node of the organ
    int getNodeId(int i) const { return nodeIds.at(i); } ///< global node index of the i-th node, i is called the local node index
    double getNodeCT(int i) const { return nodeCTs.at(i); } ///< creation time of the i-th node
    void addNode(Vector3d n, double t); //< adds a node to the root
    void addNode(Vector3d n, int id, double t); //< adds a node to the root
    virtual std::vector<Vector2i> getSegments() const; // /< default, the organ is represented by a polyline

    /* last time step */
    bool hasMoved() { return moved; }; ///< have any nodes moved during the last simulate call
    virtual std::vector<int> getMovedNodeIds() const; ///< global node indices of nodes that have been moved during last time step
    virtual std::vector<Vector3d> getMovedNodes() const; ///< the new node positions corresponding to to the indices from getMovedNodeIds
    int getOldNumberOfNodes() { return oldNumberOfNodes; } ///< the number of nodes before the last simulate call

    /* for post processing */
    std::vector<Organ*> getOrgans(int ot=-1); ///< the organ including children in a sequential vector
    void getOrgans(int otype, std::vector<Organ*>& v); ///< the organ including children in a sequential vector
    virtual double getParameter(std::string name) const; ///< returns an organ parameter

    /* IO */
    virtual std::string toString() const; ///< info for debugging
    virtual void writeRSML(tinyxml2::XMLDocument& doc, tinyxml2::XMLElement* parent) const; ///< writes this organs RSML tag

protected:

    /* up and down the organ tree */
    Organism* plant; ///< the plant of which this organ is part of
    Organ* parent; ///< pointer to the parent organ (nullptr if it has no parent)
    std::vector<Organ*> children; ///< the successive organs

    /* Parameters that are constant over the organ life time */
    const int id; ///< unique organ id
    const OrganParameter* param_; ///< the parameter set of this organ

    /* Parameters are changing over time */
    bool alive = true; ///< true: alive, false: dead
    bool active = true; ///< true: active, false: organ stopped growing
    double age = 0; ///< current age [days]
    double length = 0; ///< length of the organ [cm]

    /* node data */
    std::vector<Vector3d> nodes; ///< nodes of the organ [cm]
    std::vector<int> nodeIds; ///< global node indices
    std::vector<double> nodeCTs; ///< node creation times [days]

    /* last time step */
    bool moved = false; ///< nodes moved during last time step
    int oldNumberOfNodes = 0; ///< number of nodes at the end of previous time step

};

} // namespace CRootBox

#endif
