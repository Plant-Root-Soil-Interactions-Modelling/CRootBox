// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOTSYSTEM_H_
#define ROOTSYSTEM_H_

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <random>
#include <numeric>
#include <cmath>
#include <stack>
#include <map>

#include "ModelParameter.h"
#include "Organism.h"
#include "Root.h"
#include "soil.h"

namespace CRootBox {

class Root;
class RootState;
class Tropism;
class RootSystem;

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

    RootSystemState(const RootSystem& rs);

    void restore(RootSystem& rs);

private:

    std::vector<RootState> baseRoots;  ///< Base roots of the root system

    // copy because of random generator seeds
    std::vector<Tropism*> tf;
    std::vector<GrowthFunction*> gf;
    std::vector<RootTypeParameter*> rtparam;

    double simtime = 0;
    int rid = -1; // unique root id counter
    int nid = -1; // unique node id counter
    int old_non=0;
    int old_nor=0;
    int numberOfCrowns = 0;
    bool manualSeed = false;

    mutable std::mt19937 gen;
    mutable std::uniform_real_distribution<double> UD;
    mutable std::normal_distribution<double> ND;

};

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

    friend Root;  // obviously :-)
    friend Organ;
    friend RootSystemState;

public:

    enum TropismTypes { tt_plagio = 0, tt_gravi = 1, tt_exo = 2, tt_hydro = 3 };  ///< root tropism types

    enum GrowthFunctionTypes { gft_negexp = 1, gft_linear = 2 }; // root growth function

    enum ScalarTypes { st_type = 0, st_radius = 1, st_order = 2, st_time = 3, st_length = 4, st_surface = 5, st_volume = 6, st_one = 7,
        st_userdata1 = 8, st_userdata2 = 9, st_userdata3 = 10, st_parenttype = 11,
        st_lb = 12, st_la = 13, st_nob = 14, st_r = 15, st_theta = 16, st_rlt = 17,
        st_meanln = 18, st_sdln = 19}; ///< @see RootSystem::getScalar //DEPRICATED
    static const std::vector<std::string> scalarTypeNames; ///< the corresponding names //DEPRICATED

    RootSystem(): Organism() { };
    RootSystem(const RootSystem& rs); //< copy constructor
    virtual ~RootSystem();

    // Parameter input output
    RootTypeParameter* getRootTypeParameter(int type) const;///< returns the i-th root parameter set (i=1..n)
    void setRootSystemParameter(const RootSystemParameter& rsp); ///< sets the root system parameters
    RootSystemParameter* getRootSystemParameter(); ///< gets the root system parameters

    void openFile(std::string filename, std::string subdir="modelparameter/"); ///< reads root parameter and plant parameter
    int readParameters(std::istream & cin); ///< reads root parameters from an input stream
    void writeParameters(std::ostream & os) const; ///< writes root parameters

    // Simulation
    void setGeometry(SignedDistanceFunction* geom) { geometry = geom; } ///< optionally, sets a confining geometry (call before RootSystem::initialize())
    void setSoil(SoilLookUp* soil_) { soil = soil_; } ///< optionally sets a soil for hydro tropism (call before RootSystem::initialize())
    void reset(); ///< resets the root class, keeps the root type parameters

    void initialize() override { initialize(4,5); };
    void initialize(int basal, int shootborne); ///< creates the base roots, call before simulation and after setting the plant and root parameters
    void setTropism(Tropism* tf, int rt = -1);
    void simulate(double dt, bool silence = false) override; ///< simulates root system growth for time span dt
    void simulate(); ///< simulates root system growth for the time defined in the root system parameters
    void simulate(double dt, double maxinc, ProportionalElongation* se, bool silence = false); // simulates the root system with a maximal overall elongation

    // call back functions (todo simplify)
    virtual Root* createRoot(int lt, Vector3d  h, double delay, Root* parent, double pbl, int pni);
    ///< Creates a new lateral root, overwrite or change this method to use more specialized root classes
    virtual Tropism* createTropismFunction(int tt, double N, double sigma);
    ///< Creates the tropisms, overwrite or change this method to add more tropisms TODO a vector<tropism*> might be easier to use
    virtual GrowthFunction* createGrowthFunction(int gft);
    ///< Creates the growth function per root type, overwrite or change this method to add more tropisms

    // Analysis of simulation results
    std::vector<Vector3d> getNodes() const override; ///< Copies all root system nodes into a vector
    std::vector<Vector2i> getSegments(int otype=-1) const override; ///< Copies all root system segment indices into a vector
    std::vector<double> getSegmentCTs(int otype=-1) const override; ///< Copies all node emergence times into a vector
    std::vector<Organ*> getSegmentOrigins(int otype=-1) const override; ///< Copies a pointer to the root containing the segment
    //
    int getNumberOfSegments() const { return nodeId-numberOfCrowns-1; } ///< Number of segments of the root system ((nid+1)-1) - numberOfCrowns - 1 (artificial shoot)
    int getNumberOfRoots(bool all = false) const { if (all) return rid+1; else return getRoots().size(); }
    std::vector<Root*> getRoots() const; ///< Represents the root system as sequential vector of roots and buffers the result
    std::vector<Root*> getBaseRoots() const { return baseRoots; } ///< Base roots are tap root, basal roots, and shoot borne roots
    std::vector<std::vector<Vector3d>> getPolylines() const; ///< Copies the nodes of each root into a vector return all resulting vectors
    std::vector<Vector2i> getShootSegments() const; ///< Copies the segments connecting tap, basal root, shootborne roots
    std::vector<std::vector<double>> getPolylinesNET() const; ///< Copies the node emergence times of each root into a vector and returns all resulting vectors
    std::vector<int> getRootTips() const; ///< Node indices of the root tips
    std::vector<int> getRootBases() const; ///< Node indices of the root bases

    std::vector<double> getScalar(int stype=RootSystem::st_length) const; ///< Copies a scalar root parameter that is constant per root to a vector

    // Dynamic information what happened last time step
    int getNumberOfNewNodes() const { return getNumberOfNodes()-old_non; } ///< The number of new nodes created in the previous time step (ame number as new segments)
    int getNumberOfNewRoots() const { return getRoots().size() -old_nor; }  ///< The number of new roots created in the previous time step
    std::vector<int> getUpdatedNodeIndices() const; ///< Indices of nodes that were updated in the previous time step
    std::vector<Vector3d> getUpdatedNodes() const; ///< Values of the updated nodes
    std::vector<Vector3d> getNewNodes() const; ///< Nodes created in the previous time step
    std::vector<int> getNewNodeIndices() const; ///< Node indices that were created in the previous time step
    std::vector<Vector2i> getNewSegments() const; ///< Segments created in the previous time step
    std::vector<Root*> getNewSegmentsOrigin() const; ///< Copies a pointer to the root containing the new segments
    std::vector<double> getNewNETimes() const;  ///< Node emergence times created, created in the previous time step
    std::vector<double> getNewSegmentsTimes() const; ///< Segment emergence times, of segments created in the previous time step
    void push();
    void pop();

    // Output Simulation results
    void write(std::string name) const; /// writes simulation results (type is determined from file extension in name)
    void writeRSML(std::ostream & os) const; ///< writes current simulation results as RSML
    void writeVTP(std::ostream & os) const; ///< writes current simulation results as VTP (VTK polydata file)
    void writeGeometry(std::ostream & os) const; ///< writes the current confining geometry (e.g. a plant container) as paraview python script

    std::string toString() const; ///< infos about current root system state (for debugging)

    // Macro pores
    void setPoreGeometry(SignedDistanceFunction* geom) { poreGeometry = geom; }
    void setPoreConductivity(Matrix3d m) { poreConductivity = m; }
    void setPoreLocalAxes(Matrix3d m) {
        poreLocalAxes = m;
        invPoreLocalAxes = m.inverse();
    }
    Vector3d applyPoreConductivities(Vector3d v) {
        return poreLocalAxes.times(poreConductivity.times(invPoreLocalAxes.times(v))); // Landl et al. 2016, Eqn (12), K*v = [M*K'*(M^-1)]*v
    }

private:

    // todo currently everything is constant... // (move to specific application)
    Matrix3d poreConductivity = Matrix3d();
    Matrix3d poreLocalAxes= Matrix3d();
    Matrix3d invPoreLocalAxes = Matrix3d();

    const int rsmlReduction = 5; ///< only each n-th node is written to the rsml file (to coarsely adjust axial resolution for output)

    RootSystemParameter rsparam; ///< Plant parameter
    std::vector<Root*> baseRoots;  ///< Base roots of the root system

    std::map<int,GrowthFunction*> gf; ///< Growth function per root type
    std::map<int,Tropism*> tf;  ///< Tropism per root type

    SignedDistanceFunction* geometry = new SignedDistanceFunction(); ///< Confining geometry (unconfined by default)
    SoilLookUp* soil = nullptr; ///< callback for hydro, or chemo tropism (needs to set before initialize()) TODO should be a part of tf, or rtparam

    SignedDistanceFunction* poreGeometry = nullptr;

    int rid = -1; // unique root id counter

    int old_non=0;
    int old_nor=0;
    mutable std::vector<Root*> roots = std::vector<Root*>(); // buffer for getRoots()

    int numberOfCrowns = 0;

    bool manualSeed = false;

    void writeRSMLMeta(std::ostream & os) const;
    void writeRSMLPlant(std::ostream & os) const;

    int getRootIndex() { rid++; return rid; } ///< returns next unique root id, called by the constructor of Root

    std::stack<RootSystemState> stateStack;
};

} // end namespace CRootBox

#endif /* ROOTSYSTEM_H_ */
