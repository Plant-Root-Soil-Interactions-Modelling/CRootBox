#ifndef TROPISM_H
#define TROPISM_H

#include <chrono>
#include <random>

#include "Root.h"
#include "soil.h"

class Root;
class SoilLookUp;



/**
 * Base class for all tropism functions, e.g. Gravitropism, Plagiotropism, Exotropism...
 */
class Tropism
{
public:

    Tropism():Tropism(0,0) { } ///< Default constructor is TropismFunction(0,0)

    /**
     * Tropism with n_ number of trials and standard deviation of sigma_
     *
     * Always call the constructor, when overwriting the class!
     * Otherwise it will not work, and the mistake is hard to find.
     *
     * @param n_            number of tries
     * @param sigma_        standard deviation of angular change [1/cm]
     */
    Tropism(double n_,double sigma_): n(n_), sigma(sigma_), geometry(nullptr) { }

    virtual ~Tropism() {};

    void setGeometry(SignedDistanceFunction* geom) { geometry = geom; }
    void setTropismParameter(double n_,double sigma_) { n=n_; sigma=sigma_; }

    Vector2d getHeading(const Vector3d& pos, Matrix3d old,  double dx, const Root* root = nullptr);
    ///< changes baseTropism->getHeading() in case geometric boundaries are hit

    /**
     * The objective function of the random optimization of getHeading(). Overwrite this function to implement a tropism.
     *
     * @param pos      current root tip position
     * @param old      rotation matrix, old(:,1) is the root tip heading
     * @param a        rotation angle alpha (angular change)
     * @param b        rotation angle beta (radial change)
     * @param dx       small distance to look ahead
     * @param root     points to the root that called getHeading, just in case something else is needed (i.e. iheading for exotropism)
     *
     * \return         the value minimized by getHeading(), it should be in [0,1], in this way combination of various tropisms will be easier
     */
    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Root* root = nullptr) { std::cout << "TropismFunction::tropismObjective() not overwritten\n"; return 0; }
    ///< The objective function of the random optimization of getHeading().

    virtual Tropism* copy() { return new Tropism(*this); } ///< factory method

    static Vector3d getPosition(const Vector3d& pos, Matrix3d old, double a, double b, double dx);
    ///< Auxiliary function: Applies angles a and b and goes dx [cm] into the new direction

    // random numbers
    void setSeed(double seed) const { gen = std::mt19937(seed); } ///< Sets the seed of the random number generator
    double rand() const { return UD(gen); } ///< Uniformly distributed random number (0,1)
    double randn() const { return ND(gen); } ///< Normally distributed random number (0,1)

protected:

    Vector2d getUCHeading(const Vector3d& pos, Matrix3d old, double dx, const Root* root);
    ///< Get unconfined heading (called by getHeading), dices n times and takes the best shot (according to the objective function)

    double n; ///< Number of trials
    double sigma; ///< Standard deviation

    SignedDistanceFunction* geometry;
    const int alphaN = 20;
    const int betaN = 5;

private:

    mutable std::mt19937 gen = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());  // random stuff
    mutable std::normal_distribution<double> ND = std::normal_distribution<double>(0,1);
    mutable std::uniform_real_distribution<double> UD = std::uniform_real_distribution<double>(0,1);
};



/**
 * Gravitropism: the tendency to grow downwards
 */
class Gravitropism : public Tropism
{

public:

    Gravitropism(double n, double sigma) : Tropism(n,sigma) { } ///< @see TropismFunction

    virtual Tropism* copy() override { return new Gravitropism(*this); } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Root* root = nullptr) override {
        old.times(Matrix3d::rotX(b));
        old.times(Matrix3d::rotZ(a));
        return 0.5*(old.column(0).z+1.); // negative values point downwards, tranformed to 0..1
    }
    ///< TropismFunction::getHeading minimizes this function, @see TropismFunction::getHeading and @see TropismFunction::tropismObjective

};



/**
 * Plagiotropism: the tendency to stay in a horicontal layer
 */
class Plagiotropism : public Tropism
{

public:

    Plagiotropism(double n, double sigma) : Tropism(n,sigma) { } ///< @see TropismFunction

    virtual Tropism* copy() override { return new Plagiotropism(*this); } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Root* root = nullptr) override {
        old.times(Matrix3d::rotX(b));
        old.times(Matrix3d::rotZ(a));
        return std::abs(old.column(0).z); // 0..1
    }
    ///< getHeading() minimizes this function, @see TropismFunction

};



/**
 * Exotropism: the tendency to keep the initial heading
 */
class Exotropism : public Tropism
{

public:

    Exotropism(double n, double sigma) : Tropism(n,sigma) { } ///< @see TropismFunction

    virtual Tropism* copy() override { return new Exotropism(*this); } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Root* root = nullptr) override;
    ///< getHeading() minimizes this function, @see TropismFunction

};



/**
 * Hydrotropism (or Chemotropism, ...): the tendency to grow towards a higher saturation (or concentration, ...)
 */
class Hydrotropism : public Tropism
{

public:

    Hydrotropism(double n, double sigma, SoilLookUp* soil) : Tropism(n,sigma), soil(soil) { } ///< @see TropismFunction

    virtual Tropism* copy() override { return new Hydrotropism(*this); } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Root* root = nullptr) override;
    ///< getHeading() minimizes this function, @see TropismFunction

private:
    SoilLookUp* soil;
};



/**
 * Combined tropisms, creates a linear combination of the respective objective functions
 */
class CombinedTropism : public Tropism
{

public:

    CombinedTropism(double n, double sigma, std::vector<Tropism*> tropisms_, std::vector<double> weights_): Tropism(n,sigma), tropisms(tropisms_), weights(weights_) {
        assert(tropisms.size()>0);
        assert(weights.size()>0);
        assert(tropisms.size()==weights.size());
    } ///< linearly comibines the objective functions of multiple tropisms

    CombinedTropism(double n, double sigma, Tropism* t1, double w1, Tropism* t2, double w2);
    ///< linearly comibines the objective functions of two tropism funcitons

    CombinedTropism(CombinedTropism& t): Tropism(t), weights(t.weights) {
    	tropisms = std::vector<Tropism*>(t.tropisms.size());
    	for (size_t i=0; i<tropisms.size(); i++) {
    		tropisms[i] = t.tropisms[i]->copy();
    	}
    }

    virtual Tropism* copy() override { return new CombinedTropism(*this); } ///< copy constructor

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Root* root = nullptr) override;
    ///< getHeading() minimizes this function, @see TropismFunction

private:
    std::vector<Tropism*> tropisms;
    std::vector<double> weights;
};



#endif
