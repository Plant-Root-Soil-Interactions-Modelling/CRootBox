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

class RootSystem;



/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 */
class Root
{

public:

    Root(RootSystem* rs, int type, Vector3d pheading, double delay, Root* parent, double pbl, int pni); ///< typically called by constructor of RootSystem, or Root::createLaterals()
    virtual ~Root();

    void simulate(double dt); ///< root growth for a time span of \param dt

    double getCreationTime(double lenght); ///< creation time of a node at a length
    double getLength(double age); ///< exact length of the root
    double getAge(double length); ///< exact age of the root

    RootTypeParameter* getRootTypeParameter() const;  ///< Returns the root type parameter of the root
    double dx() const { return getRootTypeParameter()->dx; }; // returns the axial resolution

    std::vector<Root*> getRoots(); ///< return the root system as sequential vector
    void getRoots(std::vector<Root*>& v); ///< return the root system as sequential vector

    /* Nodes of the root */
    void addNode(Vector3d n,double t); //< adds a node to the root
    Vector3d getNode(int i) const { return nodes.at(i); } ///< i-th node of the root
    double getNodeETime(int i) const { return netimes.at(i); } ///< creation time of i-th node
    int getNodeId(int i) const {return nodeIds.at(i); } ///< unique identifier of i-th node
    size_t getNumberOfNodes() const {return nodes.size(); }  ///< return the number of the nodes of the root

    /* IO */
    void writeRSML(std::ostream & cout, std::string indent) const; ///< writes a RSML root tag
    std::string toString() const;

    RootSystem* rootsystem; ///< the root system this root is part of

    /* parameters that are given per root */
    RootParameter param; ///< the parameters of this root
    Vector3d iheading; ///< the initial heading of the root, when it was created
    int id; ///< unique root id, (not used so far)

    bool alive = 1; ///< true: alive, false: dead, (not implemented yet)
    bool active = 1; ///< true: active, false: stopped growing
    double length; ///< length [cm]
    double age; ///< current age [days] (redundant with length)
    int newnodes = 0; ///< number of new nodes created in the last simulation step

    /* parent */
    Root* parent; ///< pointer to the parent root (equals nullptr if it is a base root)
    double parent_base_length; ///< length [cm]
    int parent_ni; ///< parent node index

    /* kids */
    std::vector<Root*> laterals; ///< The lateral roots of this root

protected:

    void createSegments(double l); ///< creates segments of length l, called by Root::simulate()
    void createLateral(); ///< creates a new lateral, called by Root::simulate()

    /* parameters that are given per node */
    std::vector<Vector3d> nodes; ///< nodes of the root
    std::vector<int> nodeIds; ///< unique node identifier
    std::vector<double> netimes; ///< node emergence times [days]

};

#endif /* ROOT_H_ */
