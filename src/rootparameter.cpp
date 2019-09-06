// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "rootparameter.h"

#include "Organism.h"

#include <cmath>
#include <iostream>
#include <chrono>
#include <assert.h>

namespace CRootBox {

/**
 * @return Mean maximal root length of this root type
 */
double RootSpecificParameter::getK() const {
    double l = std::accumulate(ln.begin(), ln.end(), 0.);
    return l+la+lb;
}

/**
 * @copydoc OrganParameter::toString()
 */
std::string RootSpecificParameter::toString() const
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
RootRandomParameter::RootRandomParameter(Organism* p) :OrganRandomParameter(p)
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
RootRandomParameter::~RootRandomParameter()
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
OrganRandomParameter* RootRandomParameter::copy(Organism* p)
{
    RootRandomParameter* r = new RootRandomParameter(*this); // copy constructor breaks class introspection
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
OrganSpecificParameter* RootRandomParameter::realize()
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
    OrganSpecificParameter* p = new RootSpecificParameter(subType,lb_,la_,ln_,nob_,r_,a_,theta_,rlt_);
    return p;
}

/**
 * Choose (dice) lateral type based on root parameters successor and successorP
 *
 * @param pos       spatial position (for coupling to a soil model)
 * @return          root sub type of the lateral root
 */
int RootRandomParameter::getLateralType(const Vector3d& pos)
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
std::string RootRandomParameter::toString(bool verbose) const {

    if (verbose) {
        std::string s = OrganRandomParameter::toString(true);
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
        return OrganRandomParameter::toString(false);
    }

}

/**
 * @copydoc OrganTypeParameter::readXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 *
 * If the parameter successor or successorP are not in the element, they are set to zero size.
 */
void RootRandomParameter::readXML(tinyxml2::XMLElement* element)
{
    OrganRandomParameter::readXML(element);
    tinyxml2::XMLElement* p = element->FirstChildElement("parameter");
    successor.resize(0);
    successorP.resize(0);
    while(p) {
        std::string key = p->Attribute("name");
        if (key.compare("successor")==0)  {
            successor.push_back(p->IntAttribute("type"));
            successorP.push_back(p->DoubleAttribute("percentage"));
        }
        p = p->NextSiblingElement("parameter");
    }
    double p_ = std::accumulate(successorP.begin(), successorP.end(), 0.);
    if (p_<1) {
        std::cout << "RootRandomParameter::readXML: Warning! percentages to not add up to 1. \n";
    }
    assert(successor.size()==successorP.size() &&
        "RootTypeParameter::readXML: Successor sub type and probability vector does not have the same size" );
}

/**
 * @copydoc OrganTypeParameter::writeXML()
 *
 * We need to add the parameters that are not in the hashmaps (i.e. successor, and successorP)
 */
tinyxml2::XMLElement* RootRandomParameter::writeXML(tinyxml2::XMLDocument& doc, bool comments) const
{
    assert(successor.size()==successorP.size() &&
        "RootTypeParameter::writeXML: Successor sub type and probability vector does not have the same size" );
    tinyxml2::XMLElement* element = OrganRandomParameter::writeXML(doc, comments);
    for (int i = 0; i<successor.size(); i++) {
        tinyxml2::XMLElement* p = doc.NewElement("parameter");
        p->SetAttribute("name", "successor");
        p->SetAttribute("number", i);
        p->SetAttribute("type", successor[i]);
        p->SetAttribute("percentage", float(successorP[i]));
        element->InsertEndChild(p);
        if (comments) {
            std::string str = description.at("successor");
            tinyxml2::XMLComment* c = doc.NewComment(str.c_str());
            element->InsertEndChild(c);
        }

    }
    double p_ = std::accumulate(successorP.begin(), successorP.end(), 0.);
    if (p_<1) {
        std::cout << "RootRandomParameter::writeXML: Warning! percentages to not add up to 1. \n";
    }
    return element;
}

/**
 * @return Mean maximal root length of this roo    description["name"]  = "Name of the sub type of the organ, e.g. small lateral";
 * t type
 *
 * CRootBox parameter reader
 * todo Depricated: use readXML instead
 */
void RootRandomParameter::read(std::istream & cin)
{
    std::cout << "RootRandomParameter::read is deprecated, use readXML instead \n";
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
 * CRootBox parameter write
 * todo Depricated: use writeXML instead
 */
void RootRandomParameter::write(std::ostream & cout) const {
    std::cout << "RootRandomParameter::write is deprecated, use writeXML instead \n";
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
 * Sets up class introspection by linking parameter names to their class members,
 * additionally adds a description for each parameter, for toString and writeXML
 */
void RootRandomParameter::bindParmaters()
{
    bindParameter("organType", &organType, "Organ type (unspecified organ = 0, seed = 1, root = 2, stem = 3, leaf = 4)");
    bindParameter("subType", &subType, "Unique identifier of this sub type");
    bindParameter("lb", &lb, "Basal zone [cm]", &lbs);
    bindParameter("la", &la, "Apical zone [cm]", &las);
    bindParameter("ln", &ln, "Inter-lateral distance [cm]", &lns);
    bindParameter("nob", &nob, "Maximal number of laterals [1]", &nobs);
    bindParameter("r", &r, "Initial growth rate [cm day-1]", &rs);
    bindParameter("a", &r, "Root radius [cm]", &as);
    bindParameter("colorR", &colorR, "Root color, red component [0.-1.]");
    bindParameter("colorG", &colorG, "Root color, green component [0.-1.]");
    bindParameter("colorB", &colorB, "Root color, blue component [0.-1.]");
    bindParameter("tropismT", &tropismT, "Type of root tropism (plagio = 0, gravi = 1, exo = 2, hydro, chemo = 3)");
    bindParameter("tropismN", &tropismN, "Number of trials of root tropism");
    bindParameter("tropismS", &tropismS, "Mean value of expected change of root tropism [1/cm]");
    bindParameter("dx", &dx, "Axial resolution [cm] (maximal segment size)");
    bindParameter("theta", &theta, "Angle between root and parent root [rad]", &thetas);
    bindParameter("rlt", &rlt, "Root life time [day]", &rlts);
    bindParameter("gf", &gf, "Growth function number [1]", &rlts);
    // other parameters (descriptions only)
    description["name"]  = "Name of the sub type of the organ, e.g. small lateral";
    description["successor"] = "Sub type of lateral roots";
    description["successorP"] = "Probability of each sub type to occur";
}

} // end namespace CRootBox
