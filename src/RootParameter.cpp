// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "RootParameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <assert.h>

namespace CRootBox {

/**
 * @return Mean maximal root length of this root type
 */
double RootParameter::getK() const {
    double l = std::accumulate(ln.begin(), ln.end(), 0.);
    return l+la+lb;
}

/**
 * @copydoc OrganParameter::toString()
 */
std::string RootParameter::toString() const
{
    std::stringstream str;
    str << "subType\t" << subType << std::endl;
    str << "lb\t" << lb << std::endl << "la\t" << la << std::endl;
    str << "nob\t" << nob << std::endl << "r\t" << r << std::endl << "a\t" << a << std::endl;
    str << "theta\t" << theta << std::endl << "rlt\t" << rlt << std::endl;
    str << "ln\t";
    for (int i=0; i<ln.size(); i++) {
        str << ln[i] << " ";
    }
    str << std::endl;
    return str.str();
}



/**
 * Default constructor sets up hashmaps for class introspection
 */
RootTypeParameter::RootTypeParameter(Organism* p) :OrganTypeParameter(p)
{
    // base class default values
    name = "undefined";
    organType = Organism::ot_root;
    subType = -1;
    f_tf = new Tropism(p);
    bindParmaters();
}

/**
 * Destructor: delete all callback functions
 */
RootTypeParameter::~RootTypeParameter()
{
    delete f_tf;
    delete f_gf;
    delete f_se;
    delete f_sa;
    delete f_sbp;
}

/**
 * @copydoc OrganTypeParameter::copy()
 */
OrganTypeParameter* RootTypeParameter::copy(Organism* p)
{
    RootTypeParameter* r = new RootTypeParameter(*this); // copy constructor breaks class introspection
    r->plant = p;
    r->bindParmaters(); // fix class introspection
    r->f_tf = f_tf->copy(p); // copy call back classes
    r->f_gf = f_gf->copy();
    r->f_se = f_se->copy();
    r->f_sa = f_sa->copy();
    r->f_sbp = f_sbp->copy();
    return r;
}

/**
 * @copydoc OrganTypeParameter::realize()
 *
 * Creates a specific root from the root type parameters.
 * @return Specific root parameters derived from the root type parameters
 */
OrganParameter* RootTypeParameter::realize()
{
    //& std::cout << "RootTypeParameter::realize(): subType " << subType << "\n" << std::flush;
    double lb_ = std::max(lb + plant->randn()*lbs, 0.); // length of basal zone
    double la_ = std::max(la + plant->randn()*las, 0.); // length of apical zone
    std::vector<double> ln_; // stores the inter-distances
    int nob_ = std::max(round(nob + plant->randn()*nobs), 0.); // maximal number of branches
    for (int i = 0; i<nob_-1; i++) { // create inter-root distances
        double d = std::max(ln + plant->randn()*lns, 1.e-5); // miminum is 1.e-5
        ln_.push_back(d);
    }
    double r_ = std::max(r + plant->randn()*rs, 0.); // initial elongation
    double a_ = std::max(a + plant->randn()*as, 0.); // radius
    double theta_ = std::max(theta + plant->randn()*thetas, 0.); // initial elongation
    double rlt_ = std::max(rlt + plant->randn()*rlts, 0.); // root life time
    OrganParameter* p = new RootParameter(subType,lb_,la_,ln_,nob_,r_,a_,theta_,rlt_);
    return p;
}

/**
 * Choose (dice) lateral type based on root parameters successor and successorP
 *
 * @param pos       spatial position (for coupling to a soil model)
 * @return          root sub type of the lateral root
 */
int RootTypeParameter::getLateralType(const Vector3d& pos)
{
    assert(successor.size()==successorP.size()
        && "RootTypeParameter::getLateralType: Successor sub type and probability vector does not have the same size");
    if (successorP.size()>0) { // at least 1 successor type
        double d = plant->rand(); // in [0,1]
        int i=0;
        double p=successorP.at(i);
        i++;
        while ((p<d) && (i<successorP.size())) {
            p+=successorP.at(i);
            i++;
        }
        if (p>=d) { // success
            // std::cout << "lateral type " << successor.at(i-1) << "\n" << std::flush;
            return successor.at(i-1);
        } else { // no successors
            // std::cout << "no lateral type " << std::flush;
            return -1;
        }
    } else {
        return -1; // no successors
    }
}

/**
 * @copydoc OrganTypeParameter::toString()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
std::string RootTypeParameter::toString(bool verbose) const {

    if (verbose) {
        std::string s = OrganTypeParameter::toString(true);
        std::stringstream str;
        str << "successor\t";
        for (int i=0; i<successor.size(); i++) {
            str << successor[i] << " ";
        }
        str << "\t" << description.at("successor") << std::endl;
        str << "successorP\t";
        for (int i=0; i<successorP.size(); i++) {
            str << successorP[i] << " ";
        }
        str << "\t" << description.at("successorP") << std::endl;
        return s.insert(s.length()-4, str.str());
    } else {
        return OrganTypeParameter::toString(false);
    }

}

/**
 * @copydoc OrganTypeParameter::readXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 *
 * If the parameter successor or successorP are not in the element, they are set to zero size.
 */
void RootTypeParameter::readXML(tinyxml2::XMLElement* element)
{
    OrganTypeParameter::readXML(element);
    tinyxml2::XMLElement* p = element->FirstChildElement("parameter");
    successor.resize(0);
    successorP.resize(0);
    while(p) {
        std::string key = p->Attribute("name");
        if (key.compare("successor")==0)  {
            std::string s = std::string(p->GetText());
            std::stringstream iss(s);
            int i;
            while (iss >> i) {
                successor.push_back(i);
            }
        }
        if (key.compare("successorP")==0)  {
            std::string s = std::string(p->GetText());
            std::stringstream iss(s);
            double d;
            while (iss >> d) {
                successorP.push_back(d);
            }
        }
        p = p->NextSiblingElement("parameter");
    }
    assert(successor.size()==successorP.size() &&
        "RootTypeParameter::readXML: Successor sub type and probability vector does not have the same size" );
}

/**
 * @copydoc OrganTypeParameter::writeXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
tinyxml2::XMLElement* RootTypeParameter::writeXML(tinyxml2::XMLDocument& doc, bool comments) const
{
    tinyxml2::XMLElement* element = OrganTypeParameter::writeXML(doc, comments);
    tinyxml2::XMLElement* p = doc.NewElement("parameter");
    p->SetAttribute("name", "successor");
    std::stringstream str;
    for (auto s: successor) {
        str << s << " ";
    }
    p->SetText(str.str().c_str());
    element->InsertEndChild(p);
    if (comments) {
        std::string str = description.at("successor");
        tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
        element->InsertEndChild(c);
    }
    tinyxml2::XMLElement* p2 = doc.NewElement("parameter");
    p2->SetAttribute("name", "successorP");
    std::stringstream str2;
    for (auto s: successorP) {
        str2 << s << " ";
    }
    p2->SetText(str2.str().c_str());
    element->InsertEndChild(p2);
    if (comments) {
        std::string str = description.at("successorP");
        tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
        element->InsertEndChild(c);
    }
    return element;
}

/**
 * CRootBox parameter write
 * todo Depricated: use writeXML instead
 */
void RootTypeParameter::write(std::ostream & cout) const {
    cout << "# Root type parameter for " << name << "\n";
    cout << "type\t" << subType << "\n" << "name\t" << name << "\n" << "lb\t"<< lb <<"\t"<< lbs << "\n" << "la\t"<< la <<"\t"<< las << "\n"
        << "ln\t" << ln << "\t" << lns << "\n" << "nob\t"<< nob <<"\t"<< nobs << "\n" << "r\t"<< r <<"\t"<< rs << "\n" <<
        "a\t" << a << "\t" << as << "\n" << "color\t"<< colorR <<"\t"<< colorG << "\t" << colorB << "\n"
        << "tropism\t"<< tropismT <<"\t"<< tropismN << "\t" << tropismS << "\n" << "dx\t" << dx << "\n" << "successor\t" << successor.size() << "\t";
    for (size_t i=0; i<successor.size(); i++) {
        cout << successor[i] << "\t";
    }
    cout << "\n" << "successorP\t" << successorP.size() <<"\t";
    for (size_t i=0; i<successorP.size(); i++) {
        cout << successorP[i] << "\t";
    }
    cout << "\n" << "theta\t" << theta << "\t" << thetas << "\n" << "rlt\t" << rlt << "\t" << rlts << "\n" << "gf\t" << gf << "\n";
}

/**
 * @return Mean maximal root length of this root type
 *
 * CRootBox parameter reader
 * todo Depricated: use readXML instead
 */
void RootTypeParameter::read(std::istream & cin) {
    char ch[256]; // dummy
    cin.getline(ch,256);
    std::string s; // dummy
    double k;
    double ks;
    cin >> s >> subType >> s >> name >> s >> lb >> lbs >> s >> la >> las >> s >> ln >> lns >> s >> k >> ks;
    cin >> s >> r >> rs >> s >> a >> as >> s >> colorR >> colorG >> colorB >> s >> tropismT >> tropismN >> tropismS >> s >> dx;
    if (ln > 0) {
        nob=  (k-la-lb)/ln+1;   //conversion, because the input file delivers the lmax value and not the nob value
        nob = std::max(nob,0.);
        nobs = (ks/k - lns/ln)*k/ln; // error propagation
        if (la>0) {
            nobs -= (las/la - lns/ln)*la/ln;
        }
        if (lb>0) {
            nobs -= (lbs/lb - lns/ln)*lb/ln;
        }
        nobs = std::max(nobs,0.);
        if (std::isnan(nobs)) {
            std::cout << "RootTypeParameter::read() nobs is nan \n";
            nobs =0;
        }
    } else {
        nob=0;
        nobs = 0;
    }
    int n;
    cin  >> s >> n;
    successor.clear();
    int is;
    for (int i=0; i<n; i++) {
        cin >> is;
        successor.push_back(is);
    }
    cin >> s >> n;
    successorP.clear();
    double ds;
    for (int i=0; i<n; i++) {
        cin >> ds;
        successorP.push_back(ds);
    }
    cin >> s >> theta >> thetas >> s >> rlt >> rlts >> s >> gf >> s;
}

/**
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void RootTypeParameter::bindParmaters()
{
    iparam["organType"] = &organType; //parameters from OrganTypeparameters
    iparam["subType"] = &subType;
    description["name"]  = "Name of the sub type of the organ, e.g. small lateral";
    description["organType"]  = "Organ type (unspecified organ = 0, seed = 1, root = 2, stem = 3, leaf = 4)";
    description["subType"]  = "Unique identifier of this sub type";
    dparam["lb"] = &lb; // parameters from RootTypeParameter
    param_sd["lb"] = &lbs;
    description["lb"] = "Basal zone [cm]";
    dparam["la"] = &la;
    param_sd["la"] = &las;
    description["la"] = "Apical zone [cm]";
    dparam["ln"] = &ln;
    param_sd["ln"] = &lns;
    description["ln"] = "Inter-lateral distance [cm]";
    dparam["nob"] = &nob;
    param_sd["nob"] = &nobs;
    description["nob"] = "Maximal number of laterals [1]";
    dparam["r"] = &r;
    param_sd["r"] = &rs;
    description["r"] = "Initial growth rate [cm day-1]";
    dparam["a"] = &a;
    param_sd["a"] = &as;
    description["a"] = "Root radius [cm]";
    dparam["colorR"] = &colorR;
    description["colorR"] = "Root color, red component [0.-1.]";
    dparam["colorG"] = &colorG;
    description["colorG"] = "Root color, green component [0.-1.]";
    dparam["colorB"] = &colorB;
    description["colorB"] = "Root color, blue component [0.-1.]";
    iparam["tropismT"] = &tropismT;
    description["tropismT"] = "Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)";
    dparam["tropismN"] = &tropismN;
    description["tropismN"] = "Number of trials of root tropism";
    dparam["tropismS"] = &tropismS;
    description["tropismS"] = "Mean value of expected change of root tropism [1/cm]";
    dparam["dx"] = &dx;
    description["dx"] = "Axial resolution [cm] (maximal segment size)";
    dparam["theta"] = &theta;
    param_sd["theta"] = &thetas;
    description["theta"] = "Angle between root and parent root [rad]";
    dparam["rlt"] = &rlt;
    param_sd["rlt"] = &rlts;
    description["rlt"] = "Root life time [day]";
    iparam["gf"] = &gf;
    description["gf"] = "Growth function (negative exponential = 1, linear = 2)";
    description["successor"] = "Sub type of lateral roots";
    description["successorP"] = "Probability of each sub type to occur";
}

} // end namespace CRootBox