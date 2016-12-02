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
    virtual double getRelativeValue(const Vector3d& pos, const Root* root = nullptr) const { return 1.; }; ///< Returns a scalar poperty of the soil

    /**
     * Returns an unscaled scalar property of the soil
     *
     * @param pos       position [cm], (normally, root->getNode(root->getNumberOfNodes()-1))
     * @param root      the root that wants to know the scalar property
     *                  in some situation this might be usefull (e.g. could increase look up speed from a unstructured mesh)
     * \return          scalar soil property
     */
    virtual double getAbsoluteValue(const Vector3d& pos, const Root* root = nullptr) const { return 1.; } ///< Returns a scalar poperty of the soil

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
     * @param max_      the maximal value of the property
     * @param min_      the minimal value of the property
     * @param slope_    scales the linear gradient of the sdf (note that |grad(sdf)|= 1)
     */
    SoilPropertySDF(SignedDistanceFunction* sdf_, double max_=1, double min_=0, double slope_=1) {
        this->sdf=sdf_;
        fmax = max_;
        fmin = min_;
        slope = slope_;
    } ///< Creates the soil property from a signed distance function

    virtual double getRelativeValue(const Vector3d& pos, const Root* root = nullptr) const override {
        return (this->getAbsoluteValue(pos,root)-fmin)/(fmax-fmin);
    } ///<@see SoilProperty::getRelativeValue

    virtual double getAbsoluteValue(const Vector3d& pos, const Root* root = nullptr) const override {
        double c = -sdf->getDist(pos)/slope*2.; ///< *(-1), because inside the geometry the value is largest
        c += (fmax-fmin)/2.; // thats the value at the boundary
        return std::max(std::min(c,fmax),fmin);
    } ///<@see SoilProperty::getAbsoluteValue

    SignedDistanceFunction* sdf;
    double fmax;
    double fmin;
    double slope;
};



#endif
