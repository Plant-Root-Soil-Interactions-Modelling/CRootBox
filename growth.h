#ifndef GROWTH_H
#define GROWTH_H

#include "Root.h"

class Root;

/**
 * Abstract base class to all growth functions: currently LinearGrowth and ExponentialGrowth
 *
 * If new classes are created, they have to be added to the vector gf in in RootSystem.
 */
class GrowthFunction
{
public:
    virtual ~GrowthFunction() {};

    /**
     * Returns root length at root age t
     *
     * @param t     root age [day]
     * @param r     initial growth rate [cm/day]
     * @param k     maximal root length [cm]
     * @param root  points to the root in case more information is needed
     *
     * \return      root length [cm]
     */
    virtual double getLength(double t, double r, double k, Root* root) const
    { throw std::runtime_error( "getLength() not implemented" ); return 0; } ///< Returns root length at root age t

    /**
     * Returns the age of a root of length l
     *
     * @param l     root length [cm]
     * @param r     initial growth rate [cm/day]
     * @param k     maximal root length [cm]
     * @param root  points to the root in case more information is needed
     *
     * \return      root age [day]
     */
    virtual double getAge(double l, double r, double k, Root* root) const
    { throw std::runtime_error( "getAge() not implemented" ); return 0; } ///< Returns the age of a root of length l
};



/**
 * LinearGrowth elongates at constant rate until the maximal length k is reached
 */
class LinearGrowth : public GrowthFunction
{
public:
    virtual double getLength(double t, double r, double k, Root* root) const { return std::min(k,r*t); } ///< @see GrowthFunction
    virtual double getAge(double l, double r, double k, Root* root)  const { return l/r; } ///< @see GrowthFunction
};


/**
 * ExponentialGrowth elongates initially at constant rate r and slows down negative exponentially towards the maximum length k is reached
 */
class ExponentialGrowth : public GrowthFunction
{
public:
    virtual double getLength(double t, double r, double k, Root* root) const { return k*(1-exp(-(r/k)*t)); } ///< @see GrowthFunction
    virtual double getAge(double l, double r, double k, Root* root) const {
        if (l>(0.999*k)) { // 0.999*k is reached in finite time
            l=0.999*k;
        }
        return - k/r*log(1-l/k);
    } ///< @see GrowthFunction
};



#endif
