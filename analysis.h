#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "RootSystem.h"
#include <set>


/**
 * Meshfree analysis of the root system based on signed distance functions
 *
 * for a start this class should contain what Shehan needs for his experimental set-up
 */
class AnalysisSDF
{

public:

    AnalysisSDF() { }; ///< Creates an empty object (use AnalysisSDF::addSegments)
    AnalysisSDF(const RootSystem& rs); ///< Creates an analyser object containing the segments from the root system
    AnalysisSDF(const AnalysisSDF& a) : nodes(a.nodes), segments(a.segments), ctimes(a.ctimes), segO(a.segO) { }; ///< Copy constructor
    virtual ~AnalysisSDF() { };

    // merge segments
    void addSegments(const RootSystem& rs); ///< adds the segments
    void addSegments(const AnalysisSDF& a); ///< adds the segments

    void crop(SignedDistanceFunction* geometry); ///< crops the data to a geometry
    void filter(int st, double min, double max); ///< filters the segments to the data @see AnalysisSDF::getScalar
    void filter(int st, double value); ///< filters the segments to the data @see AnalysisSDF::getScalar
    void pack(); ///< sorts the nodes and deletes unused nodes

    // some things we might want to know
    std::vector<double> getScalar(int st) const; ///< Returns a specific parameter per root segment @see RootSystem::ScalarType
    double getSummed(int st) const; ///< Sums up the parameter
    double getSummed(int st, SignedDistanceFunction* geometry) const; ///< Sums up the parameter within the geometry
    std::vector<double> distribution(int st, double top, double bot, int n, bool exact=false) const; ///< vertical distribution of a parameter
    std::vector<AnalysisSDF> distribution(double top, double bot, int n) const; ///< vertical distribution of a parameter
    std::vector<std::vector<double>> distribution2(int st, double top, double bot, double left, double right, int n, int m, bool exact=false) const; // 2d distribution (x,z) of a parameter
    std::vector<std::vector<AnalysisSDF>> distribution2(double top, double bot, double left, double right, int n, int m) const; // 2d distribution (x,z) of a parameter

    int getNumberOfRoots() const; ///< number of different roots

    AnalysisSDF foto(const Vector3d& pos, const Matrix3d& ons, double height) const; ///< takes a picture
    AnalysisSDF cut(const SDF_HalfPlane& plane) const; ///< cuts with a plane and returns the intersection

    // some exports
    void write(std::string name, std::vector<double> data = std::vector<double>()) const; ///< writes simulation results (type is determined from file extension in name)
    void writeVTP(std::ostream & os, std::vector<int> types = {RootSystem::st_radius}, std::vector<double> data = std::vector<double>()) const; ///< writes a VTP file
    void writeRBSegments(std::ostream & os) const; ///< Writes the segments of the root system, mimics the Matlab script getSegments()

    // auxiliary
    static Vector3d cut(Vector3d in, Vector3d out, SignedDistanceFunction* geometry); ///< intersects a line with  the geometry

    std::vector<Vector3d> nodes; ///< nodes
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<double> ctimes; ///< creation times of the segments
    std::vector<Root*> segO; ///< to look up things

};

#endif
