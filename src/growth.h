#ifndef GROWTH_H
#define GROWTH_H

namespace CRootBox {

class Organ;

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
    virtual double getLength(double t, double r, double k, Organ* o) const
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
    virtual double getAge(double l, double r, double k, Organ* o) const
    { throw std::runtime_error( "getAge() not implemented" ); return 0; } ///< Returns the age of a root of length l


    virtual GrowthFunction* copy() { return new GrowthFunction(*this); }
};



/**
 * LinearGrowth elongates at constant rate until the maximal length k is reached
 */
class LinearGrowth : public GrowthFunction
{
public:
    double getLength(double t, double r, double k, Organ* o) const override { return std::min(k,r*t); } ///< @see GrowthFunction
    double getAge(double l, double r, double k, Organ* o)  const override { return l/r; } ///< @see GrowthFunction

    GrowthFunction* copy() override { return new LinearGrowth(*this); }
};


/**
 * ExponentialGrowth elongates initially at constant rate r and slows down negative exponentially towards the maximum length k is reached
 */
class ExponentialGrowth : public GrowthFunction
{
public:
    double getLength(double t, double r, double k, Organ* o) const override { return k*(1-exp(-(r/k)*t)); } ///< @see GrowthFunction
    double getAge(double l, double r, double k, Organ* o) const override {
        double age = - k/r*log(1-l/k);
        if (std::isfinite(age)) { // the age can not be computed when root length approaches max length
            return age;
        } else {
            return 1.e3; // very old
        }
    } ///< @see GrowthFunction

    GrowthFunction* copy() override { return new ExponentialGrowth(*this); }
};

} // end namespace CRootBox

#endif
