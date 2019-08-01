#ifndef GROWTH_H
#define GROWTH_H

#include<vector>

#include<mymath.h>

namespace CRootBox {

class Organism;
class SoilLookUp;

/**
 * Assembles Doussan linear system
 */
class Doussan {

public:

	Doussan(const Organism& o, const SoilLookUp& soil) :plant(o), soil(soil) { }
	virtual ~Doussan() { }

	/**
	 * constant kr value, or overwrite to add sense
	 */
	virtual double radialConductivity(int ot, int subType, double age) const {
		return kr;
	}

	/**
	 * constant kx value, or overwrite to add sense
	 */
	virtual double axialConductivity(int ot, int subType, double age) const {
		return kx;
	}

	void init(int ot = -1); ///< precompute segments, nodes, conductivities
	std::vector<double> getI(); ///< Indices I of sparse matrix A = sparse((I,J), V)
	std::vector<double> getJ(); ///< Indices J of sparse matrix A = sparse((I,J), V)
	std::vector<double> getV(); ///< Values V of sparse matrix A = sparse((I,J), V)
	std::vector<double> getB(); ///< RHS, of Ax = b, vector b

	double kr = 0; // the constant value of radialConductivity, if not overwritten
	double kx = 0; // the constant value of axialConductivity, if not overwritten

	Organism& plant;
	SoilLookUp& soil;

	std::vector<Vector3d> nodes;
	std::vector<Vector2i> segs;

	std::vector<double> radii;
	std::vector<double> length;
	std::vector<double> kr_; // per seg
	std::vector<double> kx_; // per seg

};


} // name space
