// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOT_H_
#define ROOT_H_

#include "mymath.h"
#include "Organ.h"
#include "Organism.h"
#include "RootParameter.h"

namespace CRootBox {

class RootState;

/**
 * Root
 *
 * Describes a single root, by a vector of nodes representing the root.
 * The method simulate() creates new nodes of this root, and lateral roots in the root's branching zone.
 *
 */
class Root :public Organ
{

    friend RootState;  //defined below

public:

    Root(int id, const OrganParameter* param, bool alive, bool active, double age, double length,
        Vector3d iheading, double pbl, int pni, bool moved= false, int oldNON = 0); // ///< creates everything from scratch
    Root(Organism* rs, int type, Vector3d pheading, double delay, Root* parent, double pbl, int pni); ///< used within simulation
    virtual ~Root() { }; ///< no need to do anything, children are deleted in ~Organ()

    Organ* copy(Organism* rs) override;  ///< deep copies the root tree

    int organType() const override { return Organism::ot_root; }; ///< returns the organs type

    /* Grow */
    void simulate(double dt, bool silence = false) override; ///< root growth for a time span of @param dt

    /* Roots as sequential list */
    double getParameter(std::string name) const override; ///< returns an organ parameter

    /* From analytical equations */
    double calcCreationTime(double lenght); ///< analytical creation (=emergence) time of a node at a length
    double calcLength(double age); ///< analytical length of the root
    double calcAge(double length); ///< analytical age of the root

    /* Abbreviations */
    RootTypeParameter* getRootTypeParameter() const;  ///< root type parameter of this root
    const RootParameter* param() const; ///< root parameter
    double dx() const { return getRootTypeParameter()->dx; } ///< returns the axial resolution

    /* IO */
    std::string toString() const override;

    /* Parameters that are given per root that are constant*/
    Vector3d iHeading; ///< the initial heading of the root, when it was created
    double parentBaseLength; ///< length [cm]
    int parentNI; ///< parent node index

protected:

    bool firstCall = true; ///< firstCall of createSegments in simulate

    void createSegments(double l, bool silence); ///< creates segments of length l, called by Root::simulate()
    virtual Vector3d getIncrement(const Vector3d& p, double sdx); ///< called by createSegments, to determine growth direction
    virtual void createLateral(bool silence); ///< creates a new lateral, called by Root::simulate()
    Vector3d heading(); ///< current growth direction of the root

    const double smallDx = 1e-6; ///< threshold value, smaller segments will be skipped (otherwise root tip direction can become NaN)

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
