// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef ROOTPARAMETER_H_
#define ROOTPARAMETER_H_

#include "mymath.h"
#include "OrganParameter.h"
#include "soil.h"
#include "growth.h"
#include "tropism.h"

namespace CRootBox {

class Organism;

/**
 * Parameters of a specific root, its created by RootTypeParameter:realize()
 */
class RootParameter :public OrganParameter
{

public:

    RootParameter(): RootParameter(-1,0.,0.,std::vector<double>(0),0,0.,0.,0.,0.) { } ///< Default constructor
    RootParameter(int type, double lb, double la, const std::vector<double>& ln, int nob, double r, double a, double theta, double rlt):
        OrganParameter(),  lb(lb), la(la), nob(nob), r(r), a(a), theta(theta), rlt(rlt), ln(ln) { subType = type; } ///< Constructor setting all parameters

    /*
     * RootBox parameters per single root
     */
    double lb;              ///< Basal zone [cm]
    double la;              ///< Apical zone [cm]
    int nob;                ///< Number of branches [1], ln.size()== nob-1 for nob>0
    double r;               ///< Initial growth rate [cm day-1]
    double a;               ///< Root radius [cm]
    double theta;           ///< Angle between root and parent root [rad]
    double rlt;             ///< Root life time [day]
    std::vector<double> ln; ///< Inter-lateral distances [cm]

    double getK() const; ///< Returns the exact maximal root length of this realization [cm]

    std::string toString() const override; ///< for debugging

};



/**
 * RootTypeParameter: contains a parameter set describing a root type
 */
class RootTypeParameter :public OrganTypeParameter
{

public:

    RootTypeParameter(Organism* plant); ///< default constructor
    virtual ~RootTypeParameter();

    OrganTypeParameter* copy(Organism* plant) override;

    OrganParameter* realize() override; ///< Creates a specific root from the root parameter set
    int getLateralType(const Vector3d& pos); ///< Choose (dice) lateral type based on root parameter set
    double getK() const { return std::max(nob-1,double(0))*ln+la+lb; }  ///< returns the mean maximal root length [cm]

    std::string toString(bool verbose = true) const override; ///< info for debugging

    void readXML(tinyxml2::XMLElement* element) override; ///< reads a single sub type organ parameter set
    tinyxml2::XMLElement* writeXML(tinyxml2::XMLDocument& doc, bool comments = true) const override; ///< writes a organ root parameter set

    // DEPRICATED
    void read(std::istream & cin); ///< reads a single root parameter set
    void write(std::ostream & cout) const; ///< writes a single root parameter set

    /*
     * RootBox parameters per root type
     */
    double lb = 0.; 	    ///< Basal zone [cm]
    double lbs = 0.;        ///< Standard deviation basal zone [cm]
    double la = 10.;	    ///< Apical zone [cm];
    double las = 0.;    	///< Standard deviation apical zone [cm];
    double ln = 1; 		    ///< Inter-lateral distance [cm]
    double lns = 0.;    	///< Standard deviation inter-lateral distance [cm]
    double nob = 0.;    	///< Number of branches [1]
    double nobs = 0.;   	///< Standard deviation number of branches [1]
    double r = 1;		    ///< Initial growth rate [cm day-1]
    double rs = 0.;	    	///< Standard deviation initial growth rate [cm day-1]
    double a = 0.1; 		///< Root radius [cm]
    double as = 0.; 		///< Standard deviation root radius [cm]
    double colorR = 0.6;	///< Root color (red)
    double colorG = 0.2;	///< Root color (green)
    double colorB = 0.2;	///< Root color (blue)
    int tropismT = 1;	    ///< Root tropism parameter (Type)
    double tropismN = 1.;   ///< Root tropism parameter (number of trials)
    double tropismS = 0.2;  ///< Root tropism parameter (mean value of expected changeg) [1/cm]
    double dx = 0.25; 		///< Maximal segment size [cm]
    double theta = 1.22; 	///< Angle between root and parent root (rad)
    double thetas= 0.; 	    ///< Standard deviation angle between root and parent root (rad)
    double rlt = 1e9;		///< Root life time (days)
    double rlts = 0.;	    ///< Standard deviation root life time (days)
    int gf = 1;			    ///< Growth function (1=negative exponential, 2=linear)
    std::vector<int> successor = std::vector<int>(0);			///< Lateral types [1]
    std::vector<double> successorP = std::vector<double>(0);  	///< Probabilities of lateral type to emerge (sum of values == 1) [1]

    /*
     * Callback functions for the Root (set up by the class RootSystem)
     */
    Tropism* f_tf;  ///< tropism function ( = new Tropism(plant) )
    GrowthFunction* f_gf = new ExponentialGrowth();; ///< growth function
    SoilLookUp* f_se = new SoilLookUp(); ///< scale elongation function
    SoilLookUp* f_sa = new SoilLookUp(); ///< scale angle function
    SoilLookUp* f_sbp = new SoilLookUp(); ///< scale branching probability function

};

} // end namespace CRootBox

#endif
