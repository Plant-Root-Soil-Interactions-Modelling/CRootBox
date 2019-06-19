// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "RootSystem.h"

#include "OrganParameter.h"
#include "Organism.h"



namespace CRootBox {

/**
 * Copy Constructor
 *
 * deep copies the root system
 * does not deep copy geometry, elongation functions, and soil (all not owned by rootsystem)
 * empties buffer
 */
RootSystem::RootSystem(const RootSystem& rs): Organism(rs), rsmlReduction(rs.rsmlReduction), rsparam(rs.rsparam),
    gf(rs.gf), tf(rs.tf), geometry(rs.geometry), soil(rs.soil), rid(rs.rid), old_non(rs.old_non), old_nor(rs.old_nor)
{
    std::cout << "Copying root system ("<<rs.baseRoots.size()<< " base roots) \n";

    // copy base Roots
    baseRoots = std::vector<Root*>(rs.baseRoots.size());
    for (size_t i=0; i<rs.baseRoots.size(); i++) {
        baseRoots[i] = new Root(*rs.baseRoots[i], this); // deep copy root tree
    }
    roots = std::vector<Root*>(0); // new empty buffer

    // deep copy tropisms
    auto it = rs.tf.begin();
    while(it != rs.tf.end()) {
        tf[it->first] = (it->second)->copy();
        ++it;
    }

    // deep copy growth
    auto it2 = rs.gf.begin();
    while(it2 != rs.gf.end()) {
        gf[it2->first] = (it2->second)->copy();
        ++it2;
    }

}

RootTypeParameter* RootSystem::getRootTypeParameter(int type) const
{
    return (RootTypeParameter*) getOrganTypeParameter(Organism::ot_root, type);
} ///< returns the i-th root parameter


std::vector<RootTypeParameter*> RootSystem::getRootTypeParameter() const
{
    std::vector<RootTypeParameter*>  otps = std::vector<RootTypeParameter*>(0);
    for (auto& otp : organParam[Organism::ot_root]) {
        otps.push_back((RootTypeParameter*)otp.second);
    }
    return otps;
} ///< returns all root parameters



void RootSystem::setRootSystemParameter(const RootSystemParameter& rsp)
{
    rsparam = rsp;
} ///< sets the root system parameters

RootSystemParameter* RootSystem::getRootSystemParameter() {
    return &rsparam;
} ///< gets the root system parameters

/**
 * Destructor
 */
RootSystem::~RootSystem()
{
    reset();
}

/**
 * Resets the root system: deletes all roots, deletes the growth functions, deletes the tropisms, sets simulation time to 0
 */
void RootSystem::reset()
{
    for(auto b :baseRoots) {
        delete b;
    }
    baseRoots.clear();
    for(auto f:gf) {
        delete f.second;
    }
    gf.clear();
    for(auto f:tf) {
        delete f.second;
    }
    tf.clear();
    simtime=0;
    rid = -1;
    nodeId = -1;
}

/**
 * Reads the root parameter from a file. Opens plant parameters with the same filename if available,
 * otherwise assumes a tap root system at position (0,0,-3).
 *
 * @param name          filename without file extension
 * @param subdir        directory ("modelparameter/" by default)
 */
void RootSystem::openFile(std::string name, std::string subdir)
{
    std::ifstream fis;
    // open root parameter
    std::string rp_name = subdir;
    rp_name.append(name);
    rp_name.append(".rparam");
    fis.open(rp_name.c_str());
    int c = 0;
    if (fis.good()) { // did it work?
        c = readParameters(fis);
        fis.close();
    } else {
        std::string s = "RootSystem::openFile() could not open root parameter file ";
        throw std::invalid_argument(s.append(rp_name));
    }
    std::cout << "Read " << c << " root type parameters \n"; // debug

    // open plant parameter
    std::string pp_name = subdir;
    pp_name.append(name);
    pp_name.append(".pparam");
    fis.open(pp_name.c_str());
    if (fis.good()) { // did it work?
        rsparam.read(fis);
        fis.close();
    } else { // create a tap root system
        std::cout << "No root system parameters found, using default tap root system \n";
        rsparam = RootSystemParameter();
    }
}

/**
 * Reads parameter from input stream (there is a Matlab script exporting these, @see writeParams.m)
 *
 * @param cin  in stream
 */
int RootSystem::readParameters(std::istream& cin)
{
    int c = 0;
    while (cin.good()) {
        RootTypeParameter* p = new RootTypeParameter(this);
        p->read(cin);
        p->organType = Organism::ot_root;
        setOrganTypeParameter(p);
        c++;
    }
    return c;
}

/**
 * Writes root parameters (for debugging)
 *
 * @param os  out stream
 */
void RootSystem::writeParameters(std::ostream& os) const
{
    for (auto& otp :organParam[Organism::ot_root]) {
        ((RootTypeParameter*)otp.second)->write(os);
    }
}

/**
 * Sets up the base roots according to the plant parameters,
 * a confining geometry, the tropism functions, and the growth functions.
 *
 * Call this method before simulation and after setting geometry, plant and root parameters
 */
void RootSystem::initialize(int basaltype, int shootbornetype)
{
    reset(); // just in case

    // introduce an extra node at nodes[0]
    getNodeIndex(); // increase node index

    // Create root system from the root system parameter
    const double maxT = 365.; // maximal simulation time
    const RootSystemParameter& rs = rsparam; // rename
    Vector3d iheading(0,0,-1);

    // Taproot
    Root* taproot = new Root(this, 1, iheading ,0, nullptr, 0, 0); // tap root has root type 1
    taproot->addNode(rs.seedPos,0);
    baseRoots.push_back(taproot);

    auto& pmap = organParam[Organism::ot_root];
    // Basal roots
    if (rs.maxB>0) {
        if (pmap.find(basaltype)==pmap.end()) { // if the type is not defined, copy tap root
            std::cout << "Basal root type #" << basaltype << " was not defined, using tap root parameters instead\n" << std::flush;
            RootTypeParameter* brtp = (RootTypeParameter*)getRootTypeParameter(1)->copy();
            brtp->subType = basaltype;
            setOrganTypeParameter(brtp);
        }
        int maxB = rs.maxB;
        if (rs.delayB>0) { // limit if possible
            maxB = std::min(maxB,int(ceil((maxT-rs.firstB)/rs.delayB))); // maximal for simtime maxT
        }
        double delay = rs.firstB;
        for (int i=0; i<maxB; i++) {
            Root* basalroot = new Root(this, basaltype, iheading ,delay, nullptr, 0, 0);
            basalroot->addNode(taproot->getNode(0), taproot->getNodeId(0),delay);
            baseRoots.push_back(basalroot);
            delay += rs.delayB;
        }
    }
    // Shoot borne roots
    if ((rs.nC>0) && (rs.delaySB<maxT)) { // if the type is not defined, copy basal root
        if (pmap.find(shootbornetype)==pmap.end()) {
            std::cout << "Shootborne root type #" << shootbornetype << " was not defined, using tap root parameters instead\n";
            RootTypeParameter* srtp = (RootTypeParameter*)getRootTypeParameter(1)->copy();
            srtp->subType = shootbornetype;
            setOrganTypeParameter(srtp);
        }
        Vector3d sbpos = rs.seedPos;
        sbpos.z=sbpos.z/2.; // half way up the mesocotyl
        numberOfCrowns = ceil((maxT-rs.firstSB)/rs.delayRC); // maximal number of root crowns
        double delay = rs.firstSB;
        for (int i=0; i<numberOfCrowns; i++) {
            Root* shootborne0 = new Root(this, shootbornetype, iheading ,delay, nullptr, 0, 0);
            // TODO fix the initial radial heading
            shootborne0->addNode(sbpos,delay);
            baseRoots.push_back(shootborne0);
            delay += rs.delaySB;
            for (int j=1; j<rs.nC; j++) {
                Root* shootborne = new Root(this, shootbornetype, iheading ,delay, nullptr, 0, 0);
                // TODO fix the initial radial heading
                shootborne->addNode(shootborne0->getNode(0), shootborne0->getNodeId(0),delay);
                baseRoots.push_back(shootborne);
                delay += rs.delaySB;
            }
            sbpos.z+=rs.nz;  // move up, for next root crown
            delay = rs.firstSB + i*rs.delayRC; // reset age
        }
    } else {
        numberOfCrowns=0;
    }

    // Create tropisms and growth functions per root type
    for (auto& p_otp :organParam[Organism::ot_root]) {
        RootTypeParameter* rtp = (RootTypeParameter*)p_otp.second;
        int type = rtp->tropismT;
        double N = rtp->tropismN;
        double sigma = rtp->tropismS;
        Tropism* tropism = this->createTropismFunction(type,N,sigma);
        tropism->setGeometry(geometry);
        tf[rtp->subType] = tropism;
        GrowthFunction* gf_ = this->createGrowthFunction(rtp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        gf[rtp->subType] = gf_;
    }

    old_non = baseRoots.size();
}

/**
 * Manually sets a tropism function for a specific or for all root types.
 * Must be called after RootSystem::initialize()
 */
void RootSystem::setTropism(Tropism* tf_, int rt)
{
    if (rt>-1) { // set for a specific root type
        tf.at(rt-1) = tf_;
    } else { // set for all root types (default)
        for (size_t i=0; i<tf.size(); i++) {
            tf.at(i) = tf_;
        }
    }
}

/**
 * Simulates root system growth for time span dt
 *
 * @param dt    	time step [days]
 * @param silence 	indicates if status is written to the console (cout) (default = false)
 */
void RootSystem::simulate(double dt, bool silence)
{
    if (!silence) {
        std::cout << "RootSystem.simulate(dt) from "<< simtime << " to " << simtime+dt << " days \n";
    }
    old_non = getNumberOfNodes();
    old_nor = getRoots().size();
    for (const auto& r : baseRoots) {
        r->simulate(dt, silence);
    }
    simtime+=dt;
    roots.clear(); // empty buffer
}

/**
 * Simulates root system growth for the time span defined in the parameter set
 */
void RootSystem::simulate()
{
    this->simulate(rsparam.simtime);
}

/**
 * Simulates root system growth for the time span dt [days],
 * elongates a maximum of maxinc total length [cm/day]
 * using the proportional elongation se to impede overall growth
 */
void RootSystem::simulate(double dt, double maxinc_, ProportionalElongation* se, bool silence)
{
    const double accuracy = 1.e-3;
    const int maxiter = 20;

    double maxinc = dt*maxinc_;
    auto v_ = this->getScalar(RootSystem::st_length);
    double ol = std::accumulate(v_.begin(), v_.end(), 0.0);

    int i = 0;

    push();
    se->setScale(1.);
    simulate(dt, silence);
    v_ = getScalar(RootSystem::st_length);
    double l = std::accumulate(v_.begin(), v_.end(), 0.0);
    double inc_ = l - ol;
    if (!silence) {
        std::cout << "expected increase is " << inc_ << " maximum is " << maxinc
            << "\n";
    }
    pop();

    if ((inc_>maxinc) && (std::abs(inc_-maxinc)>accuracy)) { // check if we have to perform a binary search

        double sl = 0.; // left
        double sr = 1.; // right

        while ( ((std::abs(inc_-maxinc)) > accuracy) && (i<maxiter) )  { // binary search

            double m = (sl+sr)/2.; // mid
            push();
            se->setScale(m);
            simulate(dt, silence);
            v_ = getScalar(RootSystem::st_length);
            l = std::accumulate(v_.begin(), v_.end(), 0.0);
            inc_ = l - ol;
            pop();
            if (!silence) {/**
             * Sets the seed of the root systems random number generator,
             * and all subclasses using random number generators:
             * @see TropismFunction, @see RootParameter
             * To obtain two exact same root system call before initialize().
             *
             * @param seed      random number generator seed
             */
                std::cout << "\t(sl, mid, sr) = (" << sl << ", " <<  m << ", " <<  sr << "), inc " <<  inc_ << ", err: " << std::abs(inc_-maxinc) << " > " << accuracy << "\n";
            }

            if (inc_>maxinc) { // concatenate
                sr = m;
            } else {
                sl = m;
            }
            i++;

        }
    }
    this->simulate(dt, silence);
}

/**
 * Creates a new lateral root (called by RootSystem::initialize, Root:createLateral)
 *
 * @param lt       lateral root type
 * @param h        initial heading of the new root
 * @param delay    time until the root starts to grow
 * @param parent   parent root
 * @param pbl      parent base length
 * @param pni      parent node index
 *
 */
Root* RootSystem::createRoot(int lt, Vector3d  h, double delay, Root* parent, double pbl, int pni) {
    return new Root(this,lt,h,delay,parent,pbl,pni);
}

/**
 * Creates a specific tropsim,
 * the function must be extended or overwritten to add more tropisms
 */
Tropism* RootSystem::createTropismFunction(int tt, double N, double sigma) {
    switch (tt) {
    case tt_plagio: return new Plagiotropism(N,sigma);
    case tt_gravi: return new Gravitropism(N,sigma);
    case tt_exo: return new Exotropism(N,sigma);
    case tt_hydro: {
        Tropism* gt =  new Gravitropism(N,sigma);
        Tropism* ht= new Hydrotropism(N,sigma,soil);
        Tropism* cht = new CombinedTropism(N,sigma,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
        return cht;
    }
    default: throw std::invalid_argument( "RootSystem::createTropismFunction() tropism type not implemented" );
    }
}

/**
 * Creates the possible growth functions
 * the function must bee extended or overwritten to add more growth function
 */
GrowthFunction* RootSystem::createGrowthFunction(int gft) {
    switch (gft) {
    case gft_negexp: return new ExponentialGrowth();
    case gft_linear: return new LinearGrowth();
    default: throw std::invalid_argument( "RootSystem::createGrowthFunction() growth function type not implemented" );
    }
}

/**
 * Represents the root system as sequential vector of roots,
 * copies the root only, if it has more than 1 node.
 * buffers the result, until next call of simulate(dt)
 *
 * \return sequential vector of roots with more than 1 node
 */
std::vector<Root*> RootSystem::getRoots() const
{
    if (roots.empty()) { // create buffer
        for (const auto& br : this->baseRoots) {
            br->getRoots(roots);
        }
        return roots;
    } else { // return buffer
        return roots;
    }
}

/**
 * Returns the node indices of the root tips
 */
std::vector<int> RootSystem::getRootTips() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<int> tips;
    for (auto& r : roots) {
        tips.push_back(r->getNodeId(r->getNumberOfNodes()-1));
    }
    return tips;
}

/**
 * Returns the positions of the root bases
 */
std::vector<int> RootSystem::getRootBases() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<int> bases;
    for (auto& r : roots) {
        bases.push_back(r->getNodeId(0));
    }
    return bases;
}

/**
 * Copies the nodes of the root systems into a sequential vector,
 * nodes are unique (default). See also RootSystem::getSegments
 */
std::vector<Vector3d> RootSystem::getNodes() const
{
    this->getRoots(); // update roots (if necessary)
    int non = getNumberOfNodes();
    std::vector<Vector3d> nv = std::vector<Vector3d>(non); // reserve big enough vector
    // copy initial nodes (roots might not have developed)
    for (const auto& r : baseRoots) {
        nv.at(r->getNodeId(0)) = r->getNode(0);
    }
    // copy root nodes
    for (const auto& r : roots) {
        for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
            nv.at(r->getNodeId(i)) = r->getNode(i); // pray that ids are correct
        }
    }
    nv.at(0) = Vector3d(0.,0.,0.); // add artificial shoot
    return nv;
}

/**
 * Returns the root system as polylines, i.e. each root is represented by its nodes
 */
std::vector<std::vector<Vector3d>> RootSystem::getPolylines() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<std::vector<Vector3d>> nodes = std::vector<std::vector<Vector3d>>(roots.size()); // reserve big enough vector
    for (size_t j=0; j<roots.size(); j++) {
        std::vector<Vector3d>  rn = std::vector<Vector3d>(roots[j]->getNumberOfNodes());
        for (size_t i=0; i<roots[j]->getNumberOfNodes(); i++) { // loop over all nodes of all roots
            rn.at(i) = roots[j]->getNode(i);
        }
        nodes[j] = rn;
    }
    return nodes;
}

/**
 * Return the segments of the root system at the current simulation time
 */
std::vector<Vector2i> RootSystem::getSegments(int otype) const
{
    this->getRoots(); // update roots (if necessary)
    int nos=getNumberOfSegments();
    std::vector<Vector2i> s(nos);
    int c=0;
    for (const auto& r : roots) {
        for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
            Vector2i v(r->getNodeId(i),r->getNodeId(i+1));
            s.at(c) = v;
            c++;
        }
    }
    return s;
}

/**
 * Return the segments connecting tap root, basal roots, and shoot borne roots.
 *
 * The upper node represents the oldest emerged shoot-borne root, or if none, the node where
 * all basal roots emerge.
 */
std::vector<Vector2i> RootSystem::getShootSegments() const
{
    std::vector<Vector2i> seg = std::vector<Vector2i>(0);

    // connect  basal roots node (seed) to artificial shoot
    seg.push_back(Vector2i(0,baseRoots.at(0)->getNodeId(0)));

    int n1=0, n2=0;
    for (int i=0; i<numberOfCrowns-1; i++) { // connecting root crowns
        int brn = baseRoots.size()-1;
        n1 = baseRoots.at(brn-i*rsparam.nC)->getNodeId(0);
        n2 = baseRoots.at(brn-(i+1)*rsparam.nC)->getNodeId(0);
        seg.push_back(Vector2i(n1,n2));
    }
    if (numberOfCrowns>0) { // connect to basal roots node (seed)
        int ti = baseRoots.at(0)->getNodeId(0);
        seg.push_back(Vector2i(n2,ti));
    }

    return seg;
}

/**
 * Returns pointers of the roots corresponding to each segment
 */
std::vector<Organ*> RootSystem::getSegmentOrigins(int otype) const
{
    this->getRoots(); // update roots (if necessary)
    int nos=getNumberOfSegments();
    std::vector<Organ*> s(nos);
    int c=0;
    for (const auto& r : roots) {
        for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
            s.at(c) = r;
            c++;
        }
    }
    return s;
}

/**
 * Copies the node emergence times of the root system per segment or per node into a sequential vector,
 * see RootSystem::getNodes(), RootSystem::getSegments()
 */
std::vector<double> RootSystem::getSegmentCTs(int otype) const
{
    this->getRoots(); // update roots (if necessary)

    int nos = getNumberOfSegments();
    std::vector<double> netv = std::vector<double>(nos); // reserve big enough vector
    int c = 0;
    for (const auto& r : roots) {
        for (size_t i = 1; i < r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
            netv.at(c) = r->getNodeCT(i); // pray that ids are correct
            c++;
        }
    }
    return netv;

    // CT per NODE ::
    //        int non = getNumberOfNodes();
    //        std::vector<double> nv = std::vector<double>(non); // reserve big enough vector
    //        // copy initial nodes (roots might not have developed)
    //        for (const auto& r : baseRoots) {
    //            nv.at(r->getNodeId(0)) = r->getNodeCT(0);
    //        }
    //        // copy root nodes
    //        for (const auto& r : roots) {
    //            for (size_t i = 0; i < r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
    //                nv.at(r->getNodeId(i)) = r->getNodeCT(i); // pray that ids are correct
    //            }
    //        }
    //        nv.at(0) = 0; // add artificial shoot
    //        return nv;

}

/**
 *  Returns the node emergence times to the corresponding polylines, see also RootSystem::getPolylines
 */
std::vector<std::vector<double>> RootSystem::getPolylinesNET() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<std::vector<double>> times = std::vector<std::vector<double>>(roots.size()); // reserve big enough vector
    for (size_t j=0; j<roots.size(); j++) {
        std::vector<double>  rt = std::vector<double>(roots[j]->getNumberOfNodes());
        for (size_t i=0; i<roots[j]->getNumberOfNodes(); i++) {
            rt[i] = roots[j]->getNodeCT(i);
        }
        times[j] = rt;
    }
    return times;
}


/**
 * Copies a scalar that is constant per root to a sequential vector (one scalar per root).
 *
 * @param stype     a scalar type (@see RootSystem::ScalarTypes). st_time is the emergence time of the root
 */
std::vector<double> RootSystem::getScalar(int stype) const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<double> scalars(roots.size());
    for (size_t i=0; i<roots.size(); i++) {
        double value = 0;
        switch(stype) {
        case st_type:  // type
            value = roots[i]->param()->subType;
            break;
        case st_radius: // root radius
            value = roots[i]->param()->a;
            break;
        case st_order: { // root order (calculate)
            value = 0;
            Organ* r_ = roots[i];
            while (r_->getParent()!=nullptr) {
                value++;
                r_=r_->getParent();
            }
            break;
        }
        case st_time: // emergence time of the root
            value = roots[i]->getNodeCT(0);
            break;
        case st_length:
            value = roots[i]->getLength();
            break;
        case st_surface:
            value =  roots[i]->getLength()*2.*M_PI*roots[i]->param()->a;
            break;
        case st_volume:
            value =  roots[i]->getLength()*M_PI*(roots[i]->param()->a)*(roots[i]->param()->a);
            break;
        case st_one:
            value =  1;
            break;
        case st_parenttype: {
            Root* r_ = roots[i];
            if (r_->getParent()!=nullptr) {
                value = r_->getParent()->getParam()->subType;
            } else {
                value = 0;
            }
            break;
        }
        case st_lb:
            value = roots[i]->param()->lb;
            break;
        case st_la:
            value = roots[i]->param()->la;
            break;
        case st_nob:
            value = roots[i]->param()->nob;
            break;
        case st_r:
            value = roots[i]->param()->r;
            break;
        case st_theta:
            value = roots[i]->param()->theta;
            break;
        case st_rlt:
            value = roots[i]->param()->rlt;
            break;
        case st_meanln: {
            const std::vector<double>& v_ = roots[i]->param()->ln;
            value = std::accumulate(v_.begin(), v_.end(), 0.0) / v_.size();
            break;
        }
        case st_sdln: {
            const std::vector<double>& v_ = roots[i]->param()->ln;
            double mean = std::accumulate(v_.begin(), v_.end(), 0.0) / v_.size();
            double sq_sum = std::inner_product(v_.begin(), v_.end(), v_.begin(), 0.0);
            value = std::sqrt(sq_sum / v_.size() - mean * mean);
            break;
        }
        default:
            throw std::invalid_argument( "RootSystem::getScalar type not implemented" );
        }
        scalars[i]=value;
    }
    return scalars;
}

/**
 * The indices of the nodes that were updated in the last time step
 */
std::vector<int> RootSystem::getUpdatedNodeIndices() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<int> ni = std::vector<int>(0);
    for (const auto& r : roots) {
        if (r->old_non>0){
            ni.push_back(r->getNodeId(r->old_non-1));
        }
    }
    return ni;
}

/**
 * The values of the nodes that were updated in the last time step
 */
std::vector<Vector3d> RootSystem::getUpdatedNodes() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<Vector3d> nv = std::vector<Vector3d>(0);
    for (const auto& r : roots) {
        if (r->old_non>0){
            nv.push_back(r->getNode(r->old_non-1));
        }
    }
    return nv;
}

/**
 * Returns a vector of newly created nodes since the last call of RootSystem::simulate(dt),
 * to dynamically add to the old node vector, see also RootSystem::getNodes
 */
std::vector<Vector3d> RootSystem::getNewNodes() const
{
    roots.clear();
    this->getRoots(); // update roots (if necessary)
    std::vector<Vector3d> nv(this->getNumberOfNewNodes());
    for (const auto& r : roots) {
        int onon = std::abs(r->old_non);
        for (size_t i=onon; i<r->getNumberOfNodes(); i++) { // loop over all new nodes
            nv.at(r->getNodeId(i)-this->old_non) = r->getNode(i); // pray that ids are correct
        }
    }
    return nv;
}

/**
 * Returns a vector of newly created node indices since the last call of RootSystem::simulate(dt),
 * to dynamically add to the old node vector, see also RootSystem::getNodes
 */
std::vector<int> RootSystem::getNewNodeIndices() const
{
    roots.clear();
    this->getRoots(); // update roots (if necessary)
    std::vector<int> nv(this->getNumberOfNewNodes());
    for (const auto& r : roots) {
        int onon = std::abs(r->old_non);
        for (size_t i=onon; i<r->getNumberOfNodes(); i++) { // loop over all new nodes
            nv.at(r->getNodeId(i)-this->old_non) = r->getNodeId(i); // pray that ids are correct
        }
    }
    return nv;
}

/**
 * Returns a vector of newly created segments since the last call of RootSystem::simulate(dt),
 * to dynamically add to the old segment index vector, see also RootSystem::getSegments
 */
std::vector<Vector2i> RootSystem::getNewSegments() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<Vector2i> si(this->getNumberOfNewNodes());
    int c=0;
    for (const auto& r : roots) {
        int onon = std::abs(r->old_non);
        for (size_t i=onon-1; i<r->getNumberOfNodes()-1; i++) {
            Vector2i v(r->getNodeId(i),r->getNodeId(i+1));
            si.at(c) = v;
            c++;
        }
    }
    return si;
}

/**
 * Returns a vector of pointers to the root class containing the segments, for each newly created segment
 */
std::vector<Root*> RootSystem::getNewSegmentsOrigin() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<Root*> si(this->getNumberOfNewNodes());
    int c=0;
    for (auto& r:roots) {
        int onon = std::abs(r->old_non);
        for (size_t i=onon-1; i<r->getNumberOfNodes()-1; i++) {
            si.at(c) = r;
            c++;
        }
    }
    return si;
}

/**
 * Returns a vector of newly created node ermegence times since the last call of RootSystem::simulate(dt),
 * to dynamically add to the old node vector, see also RootSystem::getNodes
 */
std::vector<double> RootSystem::getNewNETimes() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<double> nodeCTs(this->getNumberOfNewNodes());
    for (const auto& r : roots) {
        int onon = std::abs(r->old_non);
        for (size_t i=onon; i<r->getNumberOfNodes(); i++) { // loop over all new nodes
            nodeCTs.at(r->getNodeId(i)-this->old_non) = r->getNodeCT(i); // pray that ids are correct

        }
    }
    return nodeCTs;
}


/**
 * Returns a vector of emergence times of newly created segmeents since the last call of RootSystem::simulate(dt)
 */
std::vector<double> RootSystem::getNewSegmentsTimes() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<double> setimes(this->getNumberOfNewNodes());
    int c=0;
    for (const auto& r : roots) {
        int onon = std::abs(r->old_non);
        for (size_t i=onon-1; i<r->getNumberOfNodes()-1; i++) {
            setimes.at(c) = r->getNodeCT(i);
            c++;
        }
    }
    return setimes;
}

void RootSystem::push()
{
    stateStack.push(RootSystemState(*this));
}

void RootSystem::pop()
{
    RootSystemState& rss = stateStack.top();
    rss.restore(*this);
    stateStack.pop();
}


/**
 * todo
 */
std::string RootSystem::toString() const
{
    return "todo";
}


/**
 * Exports the simulation results with the type from the extension in name
 * (that must be lower case)
 *
 * @param name      file name e.g. output.vtp
 */
void RootSystem::write(std::string name) const
{
    std::ofstream fos;
    fos.open(name.c_str());
    std::string ext = name.substr(name.size()-3,name.size()); // pick the right writer
    if (ext.compare("sml")==0) {
        std::cout << "writing RSML... "<< name.c_str() <<"\n";
        writeRSML(fos);
    } else if (ext.compare("vtp")==0) {
        std::cout << "writing VTP... "<< name.c_str() <<"\n";
        writeVTP(fos);
    } else if (ext.compare(".py")==0)  {
        std::cout << "writing Geometry ... "<< name.c_str() <<"\n";
        writeGeometry(fos);
    } else {
        throw std::invalid_argument("RootSystem::write(): Unkwown file type");
    }
    fos.close();
}

/**
 * Creates an RSML file
 *
 * @param os      typically a file out stream
 */
void RootSystem::writeRSML(std::ostream & os) const
{
    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // i am not using utf-8, but not sure if ISO-8859-1 is correct
    os << "<rsml>\n";
    writeRSMLMeta(os);
    os<< "<scene>\n";
    writeRSMLPlant(os);
    os << "</scene>\n";
    os << "</rsml>\n";
}

/**
 * Writes RML meta data tag
 *
 * @param os      typically a file out stream
 */
void RootSystem::writeRSMLMeta(std::ostream & os) const
{
    os << "<metadata>\n";
    os << "\t<version>" << 1 << "</version>\n";
    os << "\t<unit>" << "cm" << "</unit>\n";
    os << "\t<resolution>" << 1 << "</resolution>\n";
    // fetch time
    //    os << "<last-modified>";
    //    auto t = std::time(nullptr);
    //    auto tm = *std::localtime(&t);
    //    os << std::put_time(&tm, "%d-%m-%Y"); // %H-%M-%S" would do the job for gcc 5.0
    //    os << "</last-modified>\n";
    os << "\t<software>CRootBox</software>\n";
    os << "</metadata>\n";
}

/**
 * Writes RSML plant tag
 *
 * @param os      typically a file out stream
 */
void RootSystem::writeRSMLPlant(std::ostream & os) const
{
    os << "<plant>\n";
    for (const auto& root : baseRoots) {
        root->writeRSML(os,"");
    }
    os << "</plant>\n";
}

/**
 * Writes current simulation results as VTP (VTK polydata file),
 * where each root is represented by a polyline.
 *
 * Use SegmentAnalyser::writeVTP() for a representation based on segments,
 * e.g. for creating a movie (and run the animate.py script), or mapping values to segments
 *
 * @param os      typically a file out stream
 */
void RootSystem::writeVTP(std::ostream & os) const
{
    this->getRoots(); // update roots (if necessary)
    const auto& nodes = getPolylines();
    const auto& times = getPolylinesNET();

    os << "<?xml version=\"1.0\"?>";
    os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "<PolyData>\n";
    int non = 0; // number of nodes
    for (const auto& r : roots) {
        non += r->getNumberOfNodes();
    }
    int nol=roots.size(); // number of lines
    os << "<Piece NumberOfLines=\""<< nol << "\" NumberOfPoints=\""<<non<<"\">\n";

    // POINTDATA
    os << "<PointData Scalars=\" PointData\">\n" << "<DataArray type=\"Float32\" Name=\"time\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
    for (const auto& r: times) {
        for (const auto& t : r) {
            os << t << " ";
        }
    }
    os << "\n</DataArray>\n" << "\n</PointData>\n";

    // CELLDATA (live on the polylines)
    os << "<CellData Scalars=\" CellData\">\n";
    const size_t N = 3; // SCALARS
    int types[N] = { st_type, st_order, st_radius };
    std::string scalarTypeNames[N] = {"type", "order", "radius" };
    for (size_t i=0; i<N; i++) {
        os << "<DataArray type=\"Float32\" Name=\"" << scalarTypeNames[i] <<"\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
        auto scalars = getScalar(types[i]);
        for (auto s : scalars) {
            os << s<< " ";
        }
        os << "\n</DataArray>\n";
    }
    os << "\n</CellData>\n";

    // POINTS (=nodes)
    os << "<Points>\n"<<"<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
    for (const auto& r : nodes) {
        for (const auto& n : r) {
            os << n.x << " "<< n.y <<" "<< n.z<< " ";
        }
    }
    os << "\n</DataArray>\n"<< "</Points>\n";

    // LINES (polylines)
    os << "<Lines>\n"<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
    int c=0;
    for (const auto& r : roots) {
        for (size_t i=0; i<r->getNumberOfNodes(); i++) {
            os << c << " ";
            c++;
        }
    }
    os << "\n</DataArray>\n"<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
    c = 0;
    for (const auto& r : roots) {
        c += r->getNumberOfNodes();
        os << c << " ";
    }
    os << "\n</DataArray>\n";
    os << "\n</Lines>\n";

    os << "</Piece>\n";
    os << "</PolyData>\n" << "</VTKFile>\n";
}

/**
 * Writes the current confining geometry (e.g. a plant container) as paraview python script
 * Just adds the initial lines, before calling the method of the sdf.
 *
 * @param os      typically a file out stream
 */
void RootSystem::writeGeometry(std::ostream & os) const
{
    os << "from paraview.simple import *\n";
    os << "paraview.simple._DisableFirstRenderCameraReset()\n";
    os << "renderView1 = GetActiveViewOrCreate('RenderView')\n\n";
    geometry->writePVPScript(os);
}


RootSystemState::RootSystemState(const RootSystem& rs) : simtime(rs.simtime), rid(rs.rid), old_non(rs.old_non), old_nor(rs.old_nor),
    numberOfCrowns(rs.numberOfCrowns), manualSeed(rs.manualSeed), gen(rs.gen), UD(rs.UD), ND(rs.ND)
{
// todo copy from above
//    tf = std::vector<Tropism*>(rs.tf.size()); // deep copy tropisms
//    for (size_t i=0; i<rs.tf.size(); i++) {
//        tf[i]= rs.tf[i]->copy();
//    }
//    gf = std::vector<GrowthFunction*>(rs.gf.size()); // deep copy growth
//    for (size_t i=0; i<rs.gf.size(); i++) {
//        gf[i]= rs.gf[i]->copy();
//    }
    baseRoots = std::vector<RootState>(rs.baseRoots.size()); // store base roots
    for (size_t i=0; i<baseRoots.size(); i++) {
        baseRoots[i] = RootState(*(rs.baseRoots[i]));
    }
}


void RootSystemState::restore(RootSystem& rs)
{
    rs.roots.clear(); // clear buffer
    // rs.rtparam  = rtparam; // TODO parameters are not restored
    rs.simtime = simtime; // copy back everything
    rs.rid = rid;
    rs.nodeId = nid;
    rs.old_non = old_non;
    rs.old_nor = old_nor;
    rs.numberOfCrowns = numberOfCrowns;
    rs.manualSeed = manualSeed;
    rs.gen = gen;
    rs.UD = UD;
    rs.ND = ND;
    for (size_t i=0; i<rs.tf.size(); i++) { // restore tropism functions
        delete rs.tf[i];
    }
    for (size_t i=0; i<rs.tf.size(); i++) {
        rs.tf[i] = tf[i];
    }
    for (size_t i=0; i<rs.gf.size(); i++) { // restore growth functions
        delete rs.gf[i];
    }
    for (size_t i=0; i<rs.gf.size(); i++) {
        rs.gf[i]= gf[i];
    }
    for (size_t i=0; i<baseRoots.size(); i++) { // restore base roots
        baseRoots[i].restore(*(rs.baseRoots[i]));
    }
}

} // end namespace CRootBox
