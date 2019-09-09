// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOTSYSTEM_H_
#define ROOTSYSTEM_H_

#include <stack>
#include <fstream>

#include "soil.h"
#include "tropism.h"
#include "Organism.h"
#include "rootparameter.h"
#include "Root.h"
#include "seedparameter.h"

namespace CRootBox {

class RootState;
class RootSystemState;

/**
 * RootSystem
 *
 * This class manages model parameter, the simulation,
 * stores the base roots of the root system,
 * and offers utility functions for post processing.
 * More post processing functions can be found in the class SegmentAnalyser
 */
class RootSystem :public Organism
{

    friend RootSystemState;

public:

    enum TropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3 };  ///< root tropism types
    enum GrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // root growth function

    RootSystem(); ///< empty root system
    RootSystem(const RootSystem& rs); ///< copy constructor
    virtual ~RootSystem() { };

    /* Parameter input output */
    RootRandomParameter* getRootTypeParameter(int type) const;///< returns the i-th root parameter set (i=1..n)
    std::vector<RootRandomParameter*> getRootTypeParameter() const; ///< all root type parameters as a vector
    void setRootSystemParameter(SeedRandomParameter& rsp); ///< sets the root system parameters
    SeedRandomParameter* getRootSystemParameter(); ///< gets the root system parameters
    void openFile(std::string filename, std::string subdir="modelparameter/"); ///< reads root parameter and plant parameter
    int readParameters(std::istream & cin); ///< reads root parameters from an input stream
    void writeParameters(std::ostream & os) const; ///< writes root parameters

    /* Simulation */
    void setGeometry(SignedDistanceFunction* geom) { geometry = geom; } ///< optionally, sets a confining geometry (call before RootSystem::initialize())
    void setSoil(SoilLookUp* soil_) { soil = soil_; } ///< optionally sets a soil for hydro tropism (call before RootSystem::initialize())
    void reset(); ///< resets the root class, keeps the root type parameters
    void initialize() override { initialize(4,5); }; ///< creates the base roots, call before simulation and after setting the plant and root parameters
    void initialize(int basal, int shootborne); ///< creates the base roots, call before simulation and after setting the plant and root parameters
    void setTropism(Tropism* tf, int rt = -1); ///< sets a tropism function for a single root type or all root types (defaut)
    void simulate(double dt, bool verbose = false) override; ///< simulates root system growth for time span dt
    void simulate(); ///< simulates root system growth for the time defined in the root system parameters
    void simulate(double dt, double maxinc, ProportionalElongation* se, bool silence = false); // simulates the root system with a maximal overall elongation

    /* sequential */
    std::vector<Root*> getRoots() const; ///< represents the root system as sequential vector of roots and buffers the result

    /* call back function creation */
    void initCallbacks(); ///< sets up callback functions for tropisms and growth functions, called by initialize()
    virtual Tropism* createTropismFunction(int tt, double N, double sigma); ///< Creates the tropisms, overwrite or change this method to add more tropisms
    virtual GrowthFunction* createGrowthFunction(int gft); ///< Creates the growth function per root type, overwrite or change this method to add more tropisms

    /* Analysis of simulation results */
    int getNumberOfSegments(int ot = -1) const override { return nodeId-numberOfCrowns-1; } ///< Number of segments of the root system ((nid+1)-1) - numberOfCrowns - 1 (artificial shoot)
    int getNumberOfRoots(bool all = false) const { if (all) return organId+1; else return getRoots().size(); }
    std::vector<Vector3d> getNodes() const override;
    std::vector<Organ*> getBaseRoots() const { return baseOrgans; } ///< Base roots are tap root, basal roots, and shoot borne roots TODO
    std::vector<Vector2i> getShootSegments() const; ///< Copies the segments connecting tap, basal root, shootborne roots
    std::vector<int> getRootTips() const; ///< Node indices of the root tips
    std::vector<int> getRootBases() const; ///< Node indices of the root bases

    /* dynamics */
    void push(); ///< push current state to a stack
    void pop(); ///< retrieve previous state from stack

    /* Output */
    void write(std::string name) const; /// writes simulation results (type is determined from file extension in name)
    void writeVTP(std::ostream & os) const; ///< writes current simulation results as VTP (VTK polydata file)
    void writeGeometry(std::ostream & os) const; ///< writes the current confining geometry (e.g. a plant container) as paraview python script

    std::string toString() const override; ///< infos about current root system state (for debugging)

private:

    SeedSpecificParameter seedParam;
    SignedDistanceFunction* geometry = new SignedDistanceFunction(); ///< Confining geometry (unconfined by default)
    SoilLookUp* soil = nullptr; ///< callback for hydro, or chemo tropism (needs to set before initialize()) TODO should be a part of tf, or rtparam

    mutable std::vector<Root*> roots = std::vector<Root*>(); // buffer for getRoots()
    int numberOfCrowns = 0;

    std::stack<RootSystemState> stateStack;
};



/**
 * Sores a state of the RootSystem,
 * i.e. all data that changes over time (*), i.e. excluding node data that cannot change
 *
 * (*) excluding changes regarding RootSystemParameter, any RootTypeParameter, confining geometry, and soil
 */
class RootSystemState
{

    friend RootSystem;

public:

    RootSystemState(const RootSystem& rs); ///< create root system state from a rootsystem

    void restore(RootSystem& rs); ///< restore evolved rootsystem back to its previous state

private:

    std::vector<RootState> baseRoots;  ///< Base roots of the root system

    double simtime = 0; ///< simulation time
    int rid = -1; ///< unique root id counter
    int nid = -1; ///< unique node id counter
    int old_non=0; ///< old number of nodes
    int old_nor=0; ///< old number of roots
    int numberOfCrowns = 0; ///< old number of root crowns

    mutable std::mt19937 gen; ///< random generator state
    mutable std::uniform_real_distribution<double> UD;  ///< random generator state
    mutable std::normal_distribution<double> ND; ///< random generator state

};



/**
 * Stores a state of the root that can be restored at a later point
 * (for RootSystem::push and RootSystem::pop)
 */
class RootState {

public:

    RootState() { };

    RootState(const Root& r); ///< create the root state from a root

    void restore(Root& r); ///< restore evolved root back to its previous state

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
    Vector3d lNode = Vector3d(0.,0.,0.); ///< last node
    int lNodeId = 0; ///< last node id
    double lneTime = 0.;  ///< last creation time
    size_t non = 0; ///< number of nodes

};

} // end namespace CRootBox

#endif /* ROOTSYSTEM_H_ */
