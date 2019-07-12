// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Root.h"

namespace CRootBox {

/**
 * Constructs a root from given data.
 * The organ tree must be created, @see Organ::setPlant, Organ::setParent, Organ::addChild
 * Organ geometry must be created, @see Organ::addNode, ensure that this->getNodeId(0) == parent->getNodeId(pni)
 *
 * @param id        the organ's unique id (@see Organ::getId)
 * @param param     the organs parameters set, ownership transfers to the organ
 * @param alive     indicates if the organ is alive (@see Organ::isAlive)
 * @param active    indicates if the organ is active (@see Organ::isActive)
 * @param age       the current age of the organ (@see Organ::getAge)
 * @param length    the current length of the organ (@see Organ::getLength)
 * @param iheading  the initial heading of this root
 * @param pbl       base length of the parent root, where this root emerges
 * @param pni       local node index, where this root emerges
 * @param moved     indicates if nodes were moved in the previous time step (default = false)
 * @param oldNON    the number of nodes of the previous time step (default = 0)
 */
Root::Root(int id, const OrganParameter* param, bool alive, bool active, double age, double length,
    Vector3d iheading, double pbl, int pni, bool moved, int oldNON)
:Organ(id, param, alive, active, age, length, moved,  oldNON ), parentBaseLength(pbl), parentNI(pni)
{ }

/**
 * Constructor: Should be only called during simulation by Root::createLateral().
 * For base roots the initial node and node creation time must be set from outside
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
        double scale = getRootTypeParameter()->f_sa->getValue(parent->getNode(pni),this);
        theta*=scale;
    }
    iHeading = ons.times(Vector3d::rotAB(theta,beta)); // new initial heading
    parentBaseLength = pbl;
    parentNI = pni;
    length = 0;
    // initial node
    if (parent!=nullptr) { // the first node of the base roots must be created in RootSystem::initialize()
        addNode(parent->getNode(pni), parent->getNodeId(pni), parent->getNodeCT(pni)+delay);
    }
}

/**
 * Deep copies the organ into the new plant @param plant.
 * All laterals are deep copied, plant and parent pointers are updated.
 *
 * @param plant     the plant the copied organ will be part of
 */
Organ* Root::copy(Organism* rs)
{
    Root* r = new Root(*this); // shallow copy
    r->param_ = new RootParameter(*param()); // copy parameters
    r->parent = nullptr;
    r->plant = plant;
    for (size_t i=0; i< children.size(); i++) {
        r->children[i] = children[i]->copy(plant); // copy laterals
        r->children[i]->setParent(this);
    }
    return r;
}

/**
 * Simulates the development of the organ in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void Root::simulate(double dt, bool verbose)
{
    firstCall = true;
    moved = false;
    oldNumberOfNodes = nodes.size();

    const RootParameter& p = *param(); // rename

    if (alive) { // dead roots wont grow

        // increase age
        if (age+dt>p.rlt) { // root life time
            dt=p.rlt-age; // remaining life span
            alive = false; // this root is dead todo remaining is not simulated
        }
        age+=dt;

        // probabilistic branching model
        if ((age>0) && (age-dt<=0)) { // the root emerges in this time step
            double P = getRootTypeParameter()->f_sbp->getValue(nodes.back(),this);
            if (P<1.) { // P==1 means the lateral emerges with probability 1 (default case)
                double p = 1.-std::pow((1.-P), dt); //probability of emergence in this time step
                if (plant->rand()>p) { // not rand()<p
                    age -= dt; // the root does not emerge in this time step
                }
            }
        }

        if (age>0) { // unborn  roots have no children

            // children first (lateral roots grow even if base root is inactive)
            for (auto l:children) {
                l->simulate(dt,verbose);
            }

            if (active) {

                // length increment
                double age_ = calcAge(length); // root age
                double dt_; // time step
                if (age<dt) { // the root emerged in this time step, adjust time step
                    dt_= age;
                } else {
                    dt_=dt;
                }

                double targetlength = calcLength(age_+dt_);
                double e = targetlength-length; // unimpeded elongation in time step dt
                double scale = getRootTypeParameter()->f_se->getValue(nodes.back(),this);
                double dl = std::max(scale*e, 0.); // length increment

                // create geometry
                if (p.nob>0) { // root has children
                    // basal zone
                    if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
                        if (length+dl<=p.lb) {
                            createSegments(dl,verbose);
                            length+=dl;
                            dl=0;
                        } else {
                            double ddx = p.lb-length;
                            createSegments(ddx,verbose);
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
                                    createLateral(verbose);
                                }
                                if (length+dl<=s) { // finish within inter-lateral distance i
                                    createSegments(dl,verbose);
                                    length+=dl;
                                    dl=0;
                                } else { // grow over inter-lateral distance i
                                    double ddx = s-length;
                                    createSegments(ddx,verbose);
                                    dl-=ddx;
                                    length=s;
                                }
                            }
                        }
                        if (p.ln.size()==children.size()) { // new lateral (the last one)
                            createLateral(verbose);
                        }
                    }
                    // apical zone
                    if (dl>0) {
                        createSegments(dl,verbose);
                        length+=dl;
                    }
                } else { // no laterals
                    if (dl>0) {
                        createSegments(dl,verbose);
                        length+=dl;
                    }
                } // if lateralgetLengths
            } // if active
            active = length<(p.getK()-dx()/10); // become inactive, if final length is nearly reached
        }
    } // if alive

}

/**
 * Analytical creation (=emergence) time of a point along the root
 *
 * @param length   length along the root, where the point is located [cm]
 * @return         the analytic time when this point was reached by the growing root [day]
 */
double Root::calcCreationTime(double length)
{
    assert(length >= 0 && "Root::getCreationTime() negative length");
    double rootage = calcAge(length);
    assert(rootage >= 0 && "Root::getCreationTime() negative root age");
    if (parent!=nullptr) {
        double pl = parentBaseLength+((RootParameter*)parent->getParam())->la;
        double page=((Root*)parent)->calcCreationTime(pl);
        assert(page>=0 && "Root::getCreationTime() parent root age is negative");
        return rootage+page;
    } else {
        return rootage+nodeCTs[0];
    }
}

/**
 * Analytical length of the root at a given age
 *
 * @param age          age of the root [day]
 * @return             root length [cm]
 */
double Root::calcLength(double age)
{
    assert(age >= 0 && "Root::getLength() negative root age");
    return getRootTypeParameter()->f_gf->getLength(age,param()->r,param()->getK(),this);
}

/**
 * Analytical age of the root at a given length
 *
 * @param length   length of the root [cm]
 */
double Root::calcAge(double length)
{
    assert(length >= 0 && "Root::getAge() negative root length");
    return getRootTypeParameter()->f_gf->getAge(length,param()->r,param()->getK(),this);
}

/**
 * @return The RootTypeParameter from the plant
 */
RootTypeParameter* Root::getRootTypeParameter() const
{
    return (RootTypeParameter*)plant->getOrganTypeParameter(Organism::ot_root, param_->subType);
}

/**
 * @return Parameters of the specific root
 */
const RootParameter* Root::param() const
{
    return (const RootParameter*)param_;
}

/**
 * Creates a new lateral root and passes time overhead
 *
 * Overwrite this method to implement more spezialized root classes.
 *
 * @param verbose   turns console output on or off
 */
void Root::createLateral(bool verbose)
{
    int lt = getRootTypeParameter()->getLateralType(nodes.back());
    if (lt>0) {
        double ageLN = this->calcAge(length); // age of root when lateral node is created
        double ageLG = this->calcAge(length+param()->la); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
        double delay = ageLG-ageLN; // time the lateral has to wait
        Root* lateral = new Root(plant, lt,  heading(), delay,  this, length, nodes.size()-1);
        children.push_back(lateral);
        lateral->simulate(age-ageLN,verbose); // pass time overhead (age we want to achieve minus current age)
    }
}

/**
 * @return Current root heading
 */
Vector3d Root::heading()
{
    if (nodes.size()>1) {
        return nodes.back().minus(nodes.at(nodes.size()-2)); // getHeading(b-a)
    } else {
        return iHeading;
    }
}

/**
 *  Creates nodes and node emergence times for a length l
 *
 *  Checks that each new segments length is <= dx but >= smallDx
 *
 *  @param l        total length of the segments that are created [cm]
 *  @param verbose  turns console output on or off
 */
void Root::createSegments(double l, bool verbose)
{
    assert(l>=0);
    // shift first node to axial resolution
    int nn = nodes.size();
    if (firstCall) { // first call of createSegments (in Root::simulate)
        firstCall = false;
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
                moved = true;
                l -= newdx;
                if (l<=0) { // ==0 should be enough
                    return;
                }
            } else {
                moved = false;
            }
        } else {
            moved = false;
        }
    }
    // create n+1 new nodes
    double sl = 0; // summed length of created segment
    int n = floor(l/dx());
    for (int i = 0; i < n + 1; i++) {

        double sdx; // segment length (<=dx)
        if (i<n) {  // normal case
            sdx = dx();
        } else { // last segment
            sdx = l-n*dx();
            if (sdx<smallDx) { // quit if l is too small
                if (verbose) {
                    std::cout << "skipped small segment (<"<< smallDx << ") \n";
                }
                return;
            }
        }
        sl += sdx;

        Vector3d newdx = getIncrement(nodes.back(), sdx);
        Vector3d newnode = Vector3d(nodes.back().plus(newdx));
        double et = this->calcCreationTime(length+sl);
        et = std::max(et,plant->getSimTime());
        // in case of impeded growth the node emergence time is not exact anymore,
        // but might break down to temporal resolution
        addNode(newnode,et);

    }
}

/**
 * Returns the increment of the next segments
 *
 *  @param p       coordinates of previous node
 *  @param sdx     length of next segment [cm]
 *  @return        the vector representing the increment
 */
Vector3d Root::getIncrement(const Vector3d& p, double sdx) {
    Vector3d h; // current heading
    if (nodes.size() > 1) {
        h = nodes.back().minus(nodes.at(nodes.size() - 2));
        h.normalize();
    } else {
        h = iHeading;
    }
    Matrix3d ons = Matrix3d::ons(h);
    Vector2d ab = getRootTypeParameter()->f_tf->getHeading(p, ons, sdx, this);
    Vector3d sv = ons.times(Vector3d::rotAB(ab.x,ab.y));
    return sv.times(sdx);
    // TODO how to we solve that?! OVERWRITE in specialisation
    //    if (((RootSystem*)plant)->poreGeometry==nullptr) { // no pores defined
    //        return sv.times(sdx);
    //    } else {
    //        if (((RootSystem*)plant)->poreGeometry->getDist(p)<0) { // inside the pore
    //            auto sv1 = ((RootSystem*)plant)->applyPoreConductivities(sv);
    //            // std::cout << "Length before " << sv.length() << ", length after " << sv1.length() << "\n";
    //            sv1.normalize();
    //            return sv1.times(sdx);
    //        } else {
    //            return sv.times(sdx);
    //        }
    //    }

}

/**
 * @copydoc Organ::getParameter
 *
 * Note:
 * lnMean, and lnDev denotes the mean and standard deviation of the inter-lateral distance of this organ
 * ln_mean, and ln_dev is the mean and standard deviation from the root type parameters
 */
double Root::getParameter(std::string name) const {
    if (name=="lb") { return param()->lb; } // basal zone [cm]
    if (name=="la") { return param()->la; } // apical zone [cm]
    if (name=="nob") { return param()->nob; } // number of branches
    if (name=="r"){ return param()->r; }  // initial growth rate [cm day-1]
    if (name=="radius") { return param()->a; } // root radius [cm]
    if (name=="a") { return param()->a; } // root radius [cm]
    if (name=="theta") { return param()->theta; } // angle between root and parent root [rad]
    if (name=="rlt") { return param()->rlt; } // root life time [day]
    if (name=="k") { return param()->getK(); }; // maximal root length [cm]
    if (name=="lnMean") { // mean lateral distance [cm]
        auto& v =param()->ln;
        return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    }
    if (name=="lnDev") { // standard deviation of lateral distance [cm]
        auto& v =param()->ln;
        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
        return std::sqrt(sq_sum / v.size() - mean * mean);
    }
    if (name=="volume") { return param()->a*param()->a*M_PI*getLength(); } // // root volume [cm^3]
    if (name=="surface") { return 2*param()->a*M_PI*getLength(); }
    if (name=="type") { return this->param_->subType; }  // in CRootBox the subType is often called just type
    if (name=="iHeadingX") { return iHeading.x; } // root initial heading x - coordinate [cm]
    if (name=="iHeadingY") { return iHeading.y; } // root initial heading y - coordinate [cm]
    if (name=="iHeadingZ") { return iHeading.z; } // root initial heading z - coordinate [cm]
    if (name=="parentBaseLength") { return parentBaseLength; } // length of parent root where the lateral emerges [cm]
    if (name=="parentNI") { return parentNI; } // local parent node index where the lateral emerges
    return Organ::getParameter(name);
}

/**
 * @return Quick info about the object for debugging
 */
std::string Root::toString() const
{
    std::stringstream str;
    str << "Root #"<< id <<": type "<<param()->subType << ", length: "<< length << ", age: " <<age<<" with "
        << children.size() << " laterals" << std::endl;
    return str.str();
}



/**
 * TODO docme
 */
RootState::RootState(const Root& r): alive(r.alive), active(r.active), age(r.age), length(r.length), old_non(r.oldNumberOfNodes)
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

/**
 * TODO docme
 */
void RootState::restore(Root& r)
{
    r.alive = alive; // copy things that changed
    r.active = active;
    r.age = age;
    r.length = length;
    r.oldNumberOfNodes = old_non;
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
