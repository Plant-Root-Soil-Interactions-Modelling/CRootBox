#ifndef ANALYSIS_H_
#define ANALYSIS_H_

#include "RootSystem.h"
#include <set>

/**
 * Meshfree analysis of the root system based on signed distance functions.
 */
class SegmentAnalyser
{

public:

    SegmentAnalyser() { }; ///< creates an empty object (use AnalysisSDF::addSegments)
    SegmentAnalyser(const RootSystem& rs); ///< creates an analyser object containing the segments from the root system
    SegmentAnalyser(const SegmentAnalyser& a) : nodes(a.nodes), segments(a.segments), ctimes(a.ctimes), segO(a.segO), lengthFromOrigin(a.lengthFromOrigin) { } ///< copy constructor, does not copy user data
    virtual ~SegmentAnalyser() { }; ///< nothing to do here

    // merge segments
    void addSegments(const RootSystem& rs); ///< adds the segments
    void addSegments(const SegmentAnalyser& a); ///< adds the segments

    // reduce number of segments
    void crop(SignedDistanceFunction* geometry); ///< crops the data to a geometry
    void filter(int st, double min, double max); ///< filters the segments to the data @see AnalysisSDF::getScalar
    void filter(int st, double value); ///< filters the segments to the data @see AnalysisSDF::getScalar
    void pack(); ///< sorts the nodes and deletes unused nodes

    // some things we might want to know
    std::vector<double> getScalar(int st) const; ///< Returns a specific parameter per segment @see RootSystem::ScalarType
    double getSegmentLength(int i) const; ///< returns the length of a segment
    double getSummed(int st) const; ///< Sums up the parameter
    double getSummed(int st, SignedDistanceFunction* geometry) const; ///< Sums up the parameter within the geometry
    std::vector<double> distribution(int st, double top, double bot, int n, bool exact=false) const; ///< vertical distribution of a parameter
    std::vector<SegmentAnalyser> distribution(double top, double bot, int n) const; ///< vertical distribution of a parameter
    std::vector<std::vector<double>> distribution2(int st, double top, double bot, double left, double right, int n, int m, bool exact=false) const; ///< 2d distribution (x,z) of a parameter
    std::vector<std::vector<SegmentAnalyser>> distribution2(double top, double bot, double left, double right, int n, int m) const; ///< 2d distribution (x,z) of a parameter
    // todo distribution3

    // rather specialized things we want to know
    std::vector<Root*> getRoots() const; ///< segment origins
    int getNumberOfRoots() const; ///< number of different roots
    SegmentAnalyser foto(const Vector3d& pos, const Matrix3d& ons, double height) const; ///< takes a picture TODO unfinished, untested
    SegmentAnalyser cut(const SDF_HalfPlane& plane) const; ///< returns the segments intersecting with a plane (e.g. for trenches)

    // User data for export or distributions
    void addUserData(std::vector<double> data, std::string name) { assert(data.size()==segments.size()); userData.push_back(data); userDataNames.push_back(name);}
    ///< adds user data that are written int the VTP file, @see SegmentAnalyser::writeVTP
    void clearUserData() { userData.clear(); userDataNames.clear(); } ///< resets the user data

    // some exports
    void write(std::string name); ///< writes simulation results (type is determined from file extension in name)
    void writeVTP(std::ostream & os, std::vector<int> types = {RootSystem::st_radius}) const; ///< writes a VTP file
    void writeRBSegments(std::ostream & os) const; ///< Writes the segments of the root system, mimics the Matlab script getSegments()
    void writeDGF(std::ostream & os) const; ///< Writes the segments of the root system in DGF format used by DuMux

    // auxiliary
    static Vector3d cut(Vector3d in, Vector3d out, SignedDistanceFunction* geometry); ///< intersects a line with  the geometry

    std::vector<Vector3d> nodes; ///< nodes
    std::vector<Vector2i> segments; ///< connectivity of the nodes
    std::vector<Vector3d> nodesDGF; ///< nodes
    std::vector<Vector2i> segmentsDGF; ///< connectivity of the nodes
    std::vector<double> ctimes; ///< creation times of the segments
    std::vector<double> lengthFromOrigin; ///< distance of segment from the branch's origin
    std::vector<Root*> segO; ///< to look up things

protected:

    std::vector<std::vector<double>> userData; ///< user data attached to the segments (for vtp file), e.g. flux, pressure, etc.
    std::vector<std::string> userDataNames; ///< names of the data added, e.g. "Flux", "Pressure", etc.

};

inline bool operator==(const SegmentAnalyser& lhs, const SegmentAnalyser& rhs){ return (&lhs==&rhs); } // only address wise, needed for boost python indexing suite
inline bool operator!=(const SegmentAnalyser& lhs, const SegmentAnalyser& rhs){ return !(lhs == rhs); }

#endif
