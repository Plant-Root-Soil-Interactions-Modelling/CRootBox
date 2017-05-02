#ifndef SOIL_H
#define SOIL_H

#include "mymath.h"
#include "sdf.h"

class Root;



/**
 * Look up method for a scalar soil property
 */
class SoilProperty
{
public:
    virtual ~SoilProperty() {};


    /**
     * Returns a scalar property of the soil scaled from 0..1
     *
     * @param pos       position [cm], (normally, root->getNode(root->getNumberOfNodes()-1))
     * @param root      the root that wants to know the scalar property
     *                  in some situation this might be usefull (e.g. could increase look up speed from a unstructured mesh)
     * \return          scalar soil property
     */
    virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const { return 1.; } ///< Returns a scalar property of the soil, 1. per default

    virtual std::string toString() const { return "SoilProperty base class"; } ///< Quick info about the object for debugging

};



/**
 * A static soil property that is defined by a signed distance function
 */
class SoilPropertySDF : public SoilProperty
{
public:
    SoilPropertySDF(): SoilPropertySDF(nullptr) { } ///< Default constructor

    /**
     * Creaets the soil property from a signed distance function,
     * inside the geometry the value is largest
     *
     * @param sdf_      the signed distance function representing the geometry
     * @param max_      the maximal value of the soil property
     * @param min_      the minimal value of the soil property
     * @param slope_    scales the linear gradient of the sdf (note that |grad(sdf)|= 1)
     */
    SoilPropertySDF(SignedDistanceFunction* sdf_, double max_=1, double min_=0, double slope_=1) {
        this->sdf=sdf_;
        fmax = max_;
        fmin = min_;
        slope = slope_;
    } ///< Creates the soil property from a signed distance function

    virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const override {
        double c = -sdf->getDist(pos)/slope*2.; ///< *(-1), because inside the geometry the value is largest
        c += (fmax-fmin)/2.; // thats the value at the boundary
        return std::max(std::min(c,fmax),fmin);
    } ///< returns fmin outside of the domain and fmax inside, and a linear ascend according slope

    virtual std::string toString() const { return "SoilPropertySDF"; } ///< Quick info about the object for debugging

    SignedDistanceFunction* sdf; ///< signed distance function representing the geometry
    double fmax; ///< maximum is reached within the geometry at the distance slope
    double fmin; ///< minimum is reached outside of the geometry at the distance slope
    double slope; ///< half length of linear interpolation between fmax and fmin
};



/**
 * SoilPropertySDF scaled from 0..1
 */
class ScaledSoilPropertySDF : public SoilPropertySDF
{
public:
	virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const override {
		double v = SoilPropertySDF::getValue(pos, root);
		return (v-fmin)/(fmax-fmin);
	} ///< SoilPropertySDF::getValue() but scaled from 0 to 1

	virtual std::string toString() const { return "ScaledSoilPropertySDF"; } ///< Quick info about the object for debugging

};


/**
 *  1D Look up table (todo)
 */
class SoilProperty1DTable : public SoilProperty
{
public:

	void linspace(double min, double max, double n) {

	} ///< todo sets the mesh

	int map(double z) const {
		return 0;
	} ///< todo corresponding mapping

	virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const override {
		return data[map(pos.z)];
	} ///< todo

    virtual std::string toString() const { return "SoilProperty1DTable"; } ///< Quick info about the object for debugging

	std::vector<double> data; ///< look up data todo we will need a setter

};



#endif
