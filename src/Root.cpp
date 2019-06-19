// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Root.h"

#include "RootSystem.h"

namespace CRootBox {

/**
 * Should be only called by Root::createNewRoot().
 * For base roots the initial node and node emergence time (netime) must be set from outside
 *
 * @param rs 			points to RootSystem
 * @param type 		    type of root that is created
 * @param pheading		heading of parent root at emergence
 * @param delay 		to give apical zone of parent time to develop
 * @param parent		parent root
 * @param pbl			parent base length
 * @param pni			parent node index
 */
Root::Root(Organism* rs, int type, Vector3d heading, double delay,  Root* parent, double pbl, int pni) :Organ(rs, parent, Organism::ot_root, type, delay)
{
    double beta = 2*M_PI*plant->rand(); // initial rotation
    Matrix3d ons = Matrix3d::ons(heading);
    double theta = param()->theta;
    if (parent!=nullptr) { // scale if not a baseRoot
        double scale = getRootTypeParameter()->sa->getValue(parent->getNode(pni),this);
        theta*=scale;
    }
    this->iheading = ons.times(Vector3d::rotAB(theta,beta));  // new initial heading
    parent_base_length=pbl;
    parent_ni=pni;
    length = 0;
    // initial node
    if (parent!=nullptr) { // the first node of the base roots must be created in RootSystem::initialize()
        // otherwise, don't use addNode for the first node of the root,
        // since this node exists already and does not need a new identifier
        addNode(parent->getNode(pni), parent->getNodeId(pni), parent->getNodeCT(pni)+delay);
    }
}

/**
 * Copies the root tree
 */
Root::Root(const Root& r, Organism* rs) :Organ(r, rs),
    iheading(r.iheading), parent_base_length(r.parent_base_length), parent_ni(r.parent_ni), old_non(r.old_non), smallDx(r.smallDx)
{ }

/**
 * Simulates growth of this root for a time span dt
 *
 * @param dt       time step [day]
 * @param silence  indicates if status messages are written to the console (cout) (default = false)
 */
void Root::simulate(double dt, bool silence)
{
    old_non = 0; // is set in Root:createSegments, (the zero indicates the first call to createSegments)
    const RootParameter& p = *param(); // rename

    // increase age
    if (age+dt>p.rlt) { // root life time
        dt=p.rlt-age; // remaining life span
        alive = false; // this root is dead
    }
    age+=dt;

    if (alive) { // dead roots wont grow

        // probabilistic branching model
        if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
            double P = getRootTypeParameter()->sbp->getValue(nodes.back(),this);
            if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
                double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
                //std::cout <<P<<", "<<p<< "\n";
                if (plant->rand()>p) { // not rand()<p
                    age -= dt; // the root does not emerge in this time step
                }
            }
        }
        if (age>0) {

            // children first (lateral roots grow even if base root is inactive)
            for (auto l:children) {
                l->simulate(dt,silence);
            }

            if (active) {

                // length increment
                double age_ = calcAge(length); // root age
                double dt_; // time step
                if (age<dt) {
                    dt_= age; // i dont remember what i did here TODO
                } else {
                    dt_=dt;
                }
                double targetlength = calcLength(age_+dt_); // wouldnt legnth+dt_ make more sence  TODO

                double e = targetlength-length; //elongation in time step dt
                double scale = getRootTypeParameter()->se->getValue(nodes.back(),this);
                double dl = std::max(scale*e, 0.); // length increment

                // create geometry
                if (p.nob>0) { // root has children
                    // basal zone
                    if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
                        if (length+dl<=p.lb) {
                            createSegments(dl,silence);
                            length+=dl;
                            dl=0;
                        } else {
                            double ddx = p.lb-length;
                            createSegments(ddx,silence);
                            dl-=ddx; // ddx already has been created
                            length=p.lb;
                        }
                    }
                    // branching zone
                    if ((dl>0)&&(length>=p.lb)) {
                        double s = p.lb; // summed length
                        for (size_t i=0; ((i<p.ln.size()) && (dl>0)); i++) {
                            s+=p.ln.at(i);
                            if (length<s) {
                                if (i==children.size()) { // new lateral
                                    createLateral(silence);
                                }
                                if (length+dl<=s) { // finish within inter-lateral distance i
                                    createSegments(dl,silence);
                                    length+=dl;
                                    dl=0;
                                } else { // grow over inter-lateral distance i
                                    double ddx = s-length;
                                    createSegments(ddx,silence);
                                    dl-=ddx;
                                    length=s;
                                }
                            }
                        }
                        if (p.ln.size()==children.size()) { // new lateral (the last one)
                            createLateral(silence);
                        }
                    }
                    // apical zone
                    if (dl>0) {
                        createSegments(dl,silence);
                        length+=dl;
                    }
                } else { // no laterals
                    if (dl>0) {
                        createSegments(dl,silence);
                        length+=dl;
                    }
                } // if lateralgetLengths
            } // if active
            active = length<(p.getK()-dx()/10); // become inactive, if final length is nearly reached
        }
    } // if alive

    if (old_non==0) { // if createSegments was not called
        old_non = -nodes.size();
    }

}

/**
 * Analytical creation (=emergence) time of a node at a length along the root
 *
 * @param length   length of the root [cm]
 */
double Root::calcCreationTime(double length)
{
    assert(length >= 0 && "Root::getCreationTime() negative length");
    double rootage = calcAge(length);
    assert(rootage >= 0 && "Root::getCreationTime() negative root age");
    if (parent!=nullptr) {
        double pl = parent_base_length+((RootParameter*)parent->getParam())->la; // parent length, when this root was created TODO
        double page=((Root*)parent)->calcCreationTime(pl); // TODO
        assert(page>=0);
        return rootage+page;
    } else {
        return rootage+nodeCTs[0];
    }
}

/**
 * Analytical length of the root at a given age
 *
 * @param age          age of the root [day]
 */
double Root::calcLength(double age)
{
    assert(age >= 0 && "Root::getLength() negative root age");
    return ((RootSystem*)plant)->gf.at(param()->subType)->getLength(age,param()->r,param()->getK(),this);
}

/**
 * Analytical age of the root at a given length
 *
 * @param length   length of the root [cm]
 */
double Root::calcAge(double length)
{
    assert(length >= 0 && "Root::getAge() negative root length");
    return ((RootSystem*)plant)->gf.at(param()->subType)->getAge(length,param()->r,param()->getK(),this);
}

/**
 * Returns the RootTypeParameter from the plant
 */
RootTypeParameter* Root::getRootTypeParameter() const
{
    return (RootTypeParameter*)plant->getOrganTypeParameter(Organism::ot_root, param_->subType);
}

/**
 * Parameters of the specific root
 */
const RootParameter* Root::param() const
{
    return (const RootParameter*)param_;
}

/**
 * Creates a new lateral by calling RootSystem::createNewRoot().
 *
 * Overwrite this method to implement more spezialized root classes.
 */
void Root::createLateral(bool silence)
{
    // std::cout << "createLaterals\n" << std::flush;
    const RootParameter& p = *param(); // rename
    int lt = getRootTypeParameter()->getLateralType(nodes.back());

    if (lt>0) {

        Vector3d h; // old heading
        if (nodes.size()>1) {
            h = nodes.back().minus(nodes.at(nodes.size()-2)); // getHeading(b-a)
        } else {
            h= iheading;
        }

        double ageLN = this->calcAge(length); // age of root when lateral node is created
        double ageLG = this->calcAge(length+p.la); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
        double delay = ageLG-ageLN; // time the lateral has to wait

        Root* lateral = ((RootSystem*)plant)->createRoot(lt,  h, delay,  this, length, nodes.size()-1);
        children.push_back(lateral);
        lateral->simulate(age-ageLN,silence); // pass time overhead (age we want to achieve minus current age)
    }
    // std::cout << "createLaterals end\n" << std::flush;
}

/**
 *  Creates nodes and node emergence times for length l,
 *  and updates the root heading
 *
 *  Cecks that each new segments length is <= dx but >= ddx
 *
 *  @param l       length the root growth [cm]
 */
void Root::createSegments(double l, bool silence)
{
    // std::cout << "createSegments("<< l << ")\n";
    assert(l>0);
    double sl=0; // summed length of created segment

    // shift first node to axial resolution
    int nn = nodes.size();
    if (old_non==0) { // first call of createSegments (in Root::simulate)
        if (nn>1) {
            auto n2 = nodes.at(nn-2);
            auto n1 = nodes.at(nn-1);
            double olddx = n1.minus(n2).length(); // length of last segment
            if (olddx<dx()*0.99) { // shift node instead of creating a new node
                double newdx = std::min(dx()-olddx, l);
                double sdx = olddx + newdx; // length of new segment
                Vector3d newdxv = getIncrement(n2, sdx);
                nodes[nn - 1] = Vector3d(n2.plus(newdxv));
                double et = this->calcCreationTime(length+newdx);
                nodeCTs[nn-1] = std::max(et,plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
                old_non = nn;
                l -= newdx;
                if (l<=0) { // ==0 should be enough
                    return;
                }
            } else {
                old_non = -nn;
            }
        } else {
            old_non = -nn;
        }
    }

    if (l<smallDx) {
        if (!silence) {
            std::cout << "skipped small segment (<"<< smallDx << ") \n";
        }
        return;
    }

    // create n+1 new nodes
    int n = floor(l/dx());
    for (int i = 0; i < n + 1; i++) {

        double sdx; // segment length (<=dx)
        if (i<n) {  // normal case
            sdx = dx();
        } else { // last segment
            sdx = l-n*dx();
            if (sdx<smallDx) {
                if (!silence) {
                    std::cout << "skipped small segment (<"<< smallDx << ") \n";
                }
                return;
            }
        }
        sl += sdx; //

        Vector3d newdx = getIncrement(nodes.back(), sdx);
        Vector3d newnode = Vector3d(nodes.back().plus(newdx));
        double et = this->calcCreationTime(length+sl);
        et = std::max(et,plant->getSimTime()); // in case of impeded growth the node emergence time is not exact anymore, but might break down to temporal resolution
        addNode(newnode,et);

    } // for

}

/**
 * Returns the increment of the next segments
 *
 *  @param p       position of previous node
 *  @param sdx     length of next segment
 */
Vector3d Root::getIncrement(const Vector3d& p, double sdx) {
    Vector3d h; // current heading
    if (nodes.size() > 1) {
        h = nodes.back().minus(nodes.at(nodes.size() - 2));
        h.normalize();
    } else {
        h = iheading;
    }
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = ((RootSystem*)plant)->tf.at(param()->subType)->getHeading(p, ons, sdx, this);
    Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
    if (((RootSystem*)plant)->poreGeometry==nullptr) { // no pores defined
        return sv.times(sdx);
    } else {
        if (((RootSystem*)plant)->poreGeometry->getDist(p)<0) { // inside the pore
            auto sv1 = ((RootSystem*)plant)->applyPoreConductivities(sv);
            // std::cout << "Length before " << sv.length() << ", length after " << sv1.length() << "\n";
            sv1.normalize();
            return sv1.times(sdx);
        } else {
            return sv.times(sdx);
        }
    }
}

// Root::getPoreIncrement

/**
 * Returns the root system as sequential list,
 * copies only roots with more than 1 node.
 *
 * \return sequential list of roots
 */
std::vector<Root*> Root::getRoots()
{
    std::vector<Root*> v = std::vector<Root*>();
    getRoots(v);
    return v;
}

/**
 * Returns the root system as sequential list,
 * copies only roots with more than 1 node.
 *
 * @param v     adds the subrootsystem to this vector
 */
void Root::getRoots(std::vector<Root*>& v)
{
    if (this->nodes.size()>1) {
        v.push_back(this);
    }
    for (Organ* l : this->children) {
        ((Root*)l)->getRoots(v);
    }
}

double Root::getParameter(std::string name) const {
    double d = Organ::getParameter(name);
    if (std::isnan(d)) { // name not found
        // paramter
        if (name=="lb") { return param()->lb; } // basal zone [cm]
        if (name=="la") { return param()->la; } // apical zone [cm]
        if (name=="nob") { return param()->nob; } // number of branches
        if (name=="r"){ return param()->r; }  // initial growth rate [cm day-1]
        if (name=="radius") { return param()->a; } // root radius [cm]
        if (name=="a") { return param()->a; } // root radius [cm]
        if (name=="theta") { return param()->theta; } // angle between root and parent root [rad]
        if (name=="rlt") { return param()->rlt; } // root life time [day]
        if (name=="k") { return param()->getK(); }; // maximal root length
        // todo meanLn, sdLn
        if (name=="volume") { return param()->a*param()->a*M_PI*getLength(); }
        if (name=="surface") { return 2*param()->a*M_PI*getLength(); }
        else {
            return d; // nan
        }
        if (name=="type") { d = this->param_->subType; }  // in CRootBox the subType is often called just type
    } else {
        return d;
    }
}




/**
 * writes RSML root tag
 *
 * @param cout      typically a file out stream
 * @param indent    we care for looks
 */
void Root::writeRSML(std::ostream & cout, std::string indent) const
{
    if (this->nodes.size()>1) {
        cout << indent << "<root id=\"" <<  id << "\">\n";  // open root

        /* geometry tag */
        cout << indent << "\t<geometry>\n"; // open geometry
        cout << indent << "\t\t<polyline>\n"; // open polyline
        // polyline nodes
        cout << indent << "\t\t\t" << "<point ";
        Vector3d v = nodes.at(0);
        cout << "x=\"" << v.x << "\" y=\"" << v.z << "\" z=\"" << v.y << "\"/>\n";
        int n = ((RootSystem*)plant)->rsmlReduction;
        for (size_t i = 1; i<nodes.size()-1; i+=n) {
            cout << indent << "\t\t\t" << "<point ";
            Vector3d v = nodes.at(i);
            cout << "x=\"" << v.x << "\" y=\"" << v.z << "\" z=\"" << v.y << "\"/>\n";
        }
        cout << indent << "\t\t\t" << "<point ";
        v = nodes.at(nodes.size()-1);
        cout << "x=\"" << v.x << "\" y=\"" << v.z << "\" z=\"" << v.y << "\"/>\n";
        cout << indent << "\t\t</polyline>\n"; // close polyline
        cout << indent << "\t</geometry>\n"; // close geometry

        /* properties */
        cout << indent <<"\t<properties>\n"; // open properties
        // TODO
        cout << indent << "\t</properties>\n"; // close properties

        cout << indent << "\t<functions>\n"; // open functions
        cout << indent << "\t\t<function name='emergence_time' domain='polyline'>\n"; // open functions
        cout << indent << "\t\t\t" << "<sample>" << nodeCTs.at(0) << "</sample>\n";
        for (size_t i = 1; i<nodeCTs.size()-1; i+=n) {
            cout << indent << "\t\t\t" << "<sample>" << nodeCTs.at(i) << "</sample>\n";

        }
        cout << indent << "\t\t\t" << "<sample>" << nodeCTs.at(nodeCTs.size()-1) << "</sample>\n";

        cout << indent << "\t\t</function>\n"; // close functions
        cout << indent << "\t</functions>\n"; // close functions

        /* laterals roots */
        for (size_t i = 0; i<children.size(); i++) {
            children[i]->writeRSML(cout,indent+"\t");
        }

        cout << indent << "</root>\n"; // close root
    }
}

/**
 * Quick info about the object for debugging
 */
std::string Root::toString() const
{
    std::stringstream str;
    str << "Root #"<< id <<": type "<<param()->subType << ", length: "<< length << ", age: " <<age<<" with "<< children.size() << " laterals\n";
    return str.str();
}




RootState::RootState(const Root& r): alive(r.alive), active(r.active), age(r.age), length(r.length), old_non(r.old_non)
{
    lNode = r.nodes.back();
    lNodeId = r.nodeIds.back();
    lneTime = r.nodeCTs.back();
    non = r.nodes.size();
    laterals = std::vector<RootState>(r.children.size());
    for (size_t i=0; i<laterals.size(); i++) {
        laterals[i] = RootState(*((Root*)r.children[i]));
    }
}

void RootState::restore(Root& r)
{
    r.alive = alive; // copy things that changed
    r.active = active;
    r.age = age;
    r.length = length;
    r.old_non = old_non;
    r.nodes.resize(non); // shrink vectors
    r.nodeIds.resize(non);
    r.nodeCTs.resize(non);
    r.nodes.back() = lNode; // restore last value
    r.nodeIds.back() = lNodeId;
    r.nodeCTs.back() = lneTime;
    for (size_t i = laterals.size(); i<r.children.size(); i++) { // delete roots that have not been created
        delete r.children[i];
    }
    r.children.resize(laterals.size()); // shrink and restore laterals
    for (size_t i=0; i<laterals.size(); i++) {
        laterals[i].restore(*((Root*)r.children[i]));
    }
}

} // end namespace CRootBox
