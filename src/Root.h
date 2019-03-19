// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOT_H_
#define ROOT_H_

#include <iostream>
#include <assert.h>

#include "mymath.h"
#include "sdf.h"
#include "tropism.h"
#include "growth.h"
#include "ModelParameter.h"
#include "RootSystem.h"

namespace CRootBox {

class RootSystem;
class RootState;

/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 */
class Root
{

    friend RootSystem;
    friend RootState;

public:

    Root(RootSystem* rs, int type, Vector3d pheading, double delay, Root* parent, double pbl, int pni); ///< typically called by constructor of Root::createLaterals()
    Root(const Root& r, RootSystem& rs); ///< deep copy of the tree
    virtual ~Root();

    /* Grow */
    void simulate(double dt, bool silence = false); ///< root growth for a time span of \param dt

    /* Roots as sequential list */
    std::vector<Root*> getRoots(); ///< return the root including laterals as sequential vector
    void getRoots(std::vector<Root*>& v); ///< return the root system as sequential vector

    /* Nodes of the root */
    Vector3d getNode(int i) const { return nodes.at(i); } ///< i-th node of the root
    double getNodeETime(int i) const { return netimes.at(i); } ///< creation time of i-th node
    int getNodeId(int i) const {return nodeIds.at(i); } ///< unique identifier of i-th node
    size_t getNumberOfNodes() const {return nodes.size(); }  ///< return the number of the nodes of the root
    void addNode(Vector3d n,double t); //< adds a node to the root

    /* From analytical equations */
    double getCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
    double getLength(double age); ///< analytical length of the root
    double getAge(double length); ///< analytical age of the root

    /* Abbreviations */
    RootTypeParameter* getRootTypeParameter() const;  ///< returns the root type parameter of the root
    double dx() const { return getRootTypeParameter()->dx; } ///< returns the axial resolution

    /* IO */
    void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag
    std::string toString() const;

    /* Parameters */
    RootSystem* rootsystem; ///< the root system this root is part of

    /* Parameters that are given per root that are constant*/
    RootParameter param; ///< the parameters of this root
    Vector3d iheading; ///< the initial heading of the root, when it was created
    int id; ///< unique root id, (not used so far)
    double parent_base_length; ///< length [cm]
    int parent_ni; ///< parent node index

    /* Parameters that are given per root that may change with time */
    bool alive = 1; ///< true: alive, false: dead
    bool active = 1; ///< true: active, false: root stopped growing
    double age = 0; ///< current age [days]
    double length = 0; ///< actual length [cm] of the root. might differ from getLength(age) in case of impeded root growth
    int old_non = 1; ///< number of old nodes, the sign is positive if the last node was updated, otherwise its negative

    /* Up and down the tree*/
    Root* parent; ///< pointer to the parent root (equals nullptr if it is a base root)
    std::vector<Root*> laterals; ///< the lateral roots of this root

    const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise root tip direction can become NaN)

protected:

    void createSegments(double l, bool silence); ///< creates segments of length l, called by Root::simulate()
    virtual Vector3d getIncrement(const Vector3d& p, double sdx); ///< called by createSegments, to determine growth direction
    void createLateral(bool silence); ///< creates a new lateral, called by Root::simulate()

    /* parameters that are given per node */
    std::vector<Vector3d> nodes = std::vector<Vector3d>(0); ///< nodes of the root
    std::vector<int> nodeIds = std::vector<int>(0); ///< unique node identifier
    std::vector<double> netimes = std::vector<double>(0); ///< node emergence times [days]

};

/**
 * Stores a state of the root that can be restored at a later point
 * (for RootSystem::push and RootSystem::pop)
 */
class RootState {

public:

    RootState() { };

    RootState(const Root& r);

    void restore(Root& r);

private:

    /* parameters that are given per root that may change with time */
    bool alive = 1; ///< true: alive, false: dead
    bool active = 1; ///< true: active, false: root stopped growing
    double age = 0; ///< current age [days]
    double length = 0; ///< actual length [cm] of the root. might differ from getLength(age) in case of impeded root growth
    int old_non = 1; ///< number of old nodes, the sign is positive if the last node was updated, otherwise its negative

    /* down the root branch*/
    std::vector<RootState> laterals = std::vector<RootState>(0); ///< the lateral roots of this root

    /* last node */
    Vector3d lNode = Vector3d(0.,0.,0.);
    int lNodeId = 0;
    double lneTime = 0.;
    size_t non = 0;

};

} // end namespace CRootBox

#endif /* ROOT_H_ */
