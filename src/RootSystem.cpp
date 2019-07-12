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
RootSystem::RootSystem(const RootSystem& rs):Organism(rs), rsparam(rs.rsparam), geometry(rs.geometry), soil(rs.soil),
    old_non(rs.old_non), old_nor(rs.old_nor)
{
    std::cout << "Copying root system ("<<rs.baseOrgans.size()<< " base roots) \n";
    roots = std::vector<Root*>(0); // new empty buffer
}

/**
 * Destructor
 */
RootSystem::~RootSystem()
{
    reset();
}


/**
 * Returns the i-th root parameter of sub type @param type.
 */
RootTypeParameter* RootSystem::getRootTypeParameter(int type) const
{
    return (RootTypeParameter*) getOrganTypeParameter(Organism::ot_root, type);
}

/**
 *
 */
std::vector<RootTypeParameter*> RootSystem::getRootTypeParameter() const
{
    std::vector<RootTypeParameter*>  otps = std::vector<RootTypeParameter*>(0);
    for (auto& otp : organParam[Organism::ot_root]) {
        otps.push_back((RootTypeParameter*)otp.second);
    }
    return otps;
}

/**
 * sets the root system parameters
 */
void RootSystem::setRootSystemParameter(const RootSystemParameter& rsp)
{
    rsparam = rsp;
}

/**
 * gets the root system paramete
 */
RootSystemParameter* RootSystem::getRootSystemParameter() {
    return &rsparam;
}

/**
 * Resets the root system: deletes all roots, sets simulation time to 0
 */
void RootSystem::reset()
{
    for(auto b :baseOrgans) {
        delete b;
    }
    baseOrgans.clear();
    simtime=0;
    organId = -1;
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
    baseOrgans.push_back(taproot);

    auto& pmap = organParam[Organism::ot_root];
    // Basal roots
    if (rs.maxB>0) {
        if (pmap.find(basaltype)==pmap.end()) { // if the type is not defined, copy tap root
            std::cout << "Basal root type #" << basaltype << " was not defined, using tap root parameters instead\n" << std::flush;
            RootTypeParameter* brtp = (RootTypeParameter*)getRootTypeParameter(1)->copy(this); // TODO
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
            baseOrgans.push_back(basalroot);
            delay += rs.delayB;
        }
    }
    // Shoot borne roots
    if ((rs.nC>0) && (rs.delaySB<maxT)) { // if the type is not defined, copy basal root
        if (pmap.find(shootbornetype)==pmap.end()) {
            std::cout << "Shootborne root type #" << shootbornetype << " was not defined, using tap root parameters instead\n";
            RootTypeParameter* srtp =  (RootTypeParameter*)getRootTypeParameter(1)->copy(this);
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
            baseOrgans.push_back(shootborne0);
            delay += rs.delaySB;
            for (int j=1; j<rs.nC; j++) {
                Root* shootborne = new Root(this, shootbornetype, iheading ,delay, nullptr, 0, 0);
                // TODO fix the initial radial heading
                shootborne->addNode(shootborne0->getNode(0), shootborne0->getNodeId(0),delay);
                baseOrgans.push_back(shootborne);
                delay += rs.delaySB;
            }
            sbpos.z+=rs.nz;  // move up, for next root crown
            delay = rs.firstSB + i*rs.delayRC; // reset age
        }
    } else {
        numberOfCrowns=0;
    }
    old_non = baseOrgans.size();

    initCallbacks(); //
}

/**
 * Called by RootSystem::initialize, sets up tropism and growth functions call backs
 * using RootSystem::createTropismFunction and RootSystem::createGrowthFunction
 */
void RootSystem::initCallbacks()
{
    // Create tropisms and growth functions per root type
    for (auto& p_otp :organParam[Organism::ot_root]) {
        RootTypeParameter* rtp = (RootTypeParameter*)p_otp.second;
        Tropism* tropism = this->createTropismFunction(rtp->tropismT, rtp->tropismN, rtp->tropismS);
        tropism->setGeometry(geometry);
        delete rtp->f_tf; // delete old tropism
        rtp->f_tf = tropism; // set new one
        GrowthFunction* gf_ = this->createGrowthFunction(rtp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        delete rtp->f_gf;
        rtp->f_gf  = gf_;
    }
}

/**
 * Manually sets a tropism function for a specific or for all root types.
 * Must be called after RootSystem::initialize()
 */
void RootSystem::setTropism(Tropism* tf_, int rt)
{
    if (rt>-1) { // set for a specific root type
        getRootTypeParameter(rt)->f_tf=tf_;
    } else { // set for all root types (default)
        for (auto& p_otp :organParam[Organism::ot_root]) {
            RootTypeParameter* rtp = (RootTypeParameter*)p_otp.second;
            rtp->f_tf = tf_;
        }
    }
}

/**
 * Simulates root system growth for time span dt
 *
 * @param dt    	time step [day]
 * @param silence 	indicates if status is written to the console (cout) (default = false)
 */
void RootSystem::simulate(double dt, bool verbose)
{
    old_non = getNumberOfNodes();
    old_nor = getRoots().size();
    Organism::simulate(dt,verbose);
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

    double ol = getSummed("length");
    int i = 0;

    push();
    se->setScale(1.);
    simulate(dt, silence);
    double l = getSummed("length");
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
            l = getSummed("length");
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
 * Creates a specific tropsim,
 * the function must be extended or overwritten to add more tropisms
 */
Tropism* RootSystem::createTropismFunction(int tt, double N, double sigma) {
    switch (tt) {
    case tt_plagio: return new Plagiotropism(this,N,sigma);
    case tt_gravi: return new Gravitropism(this,N,sigma);
    case tt_exo: return new Exotropism(this,N,sigma);
    case tt_hydro: {
        Tropism* gt =  new Gravitropism(this,N,sigma);
        Tropism* ht= new Hydrotropism(this,N,sigma,soil);
        Tropism* cht = new CombinedTropism(this,N,sigma,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
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
        std::vector<Organ*> organs;
        for (const auto& br : this->baseOrgans) {
            ((Root*)br)->getOrgans(ot_root, organs);
        }
        for (auto& o :organs) {
            roots.push_back((Root*)o);
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
    auto v = Organism::getNodes();
    v.at(0) = Vector3d(0.,0.,0.); // add artificial shoot
    return v;
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
    seg.push_back(Vector2i(0,baseOrgans.at(0)->getNodeId(0)));

    int n1=0, n2=0;
    for (int i=0; i<numberOfCrowns-1; i++) { // connecting root crowns
        int brn = baseOrgans.size()-1;
        n1 = baseOrgans.at(brn-i*rsparam.nC)->getNodeId(0);
        n2 = baseOrgans.at(brn-(i+1)*rsparam.nC)->getNodeId(0);
        seg.push_back(Vector2i(n1,n2));
    }
    if (numberOfCrowns>0) { // connect to basal roots node (seed)
        int ti = baseOrgans.at(0)->getNodeId(0);
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
 * The indices of the nodes that were updated in the last time step
 */
std::vector<int> RootSystem::getUpdatedNodeIndices() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<int> ni = std::vector<int>(0);
    for (const auto& r : roots) {
        if (r->getOldNumberOfNodes()>0){
            ni.push_back(r->getNodeId(r->getOldNumberOfNodes()-1));
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
        if (r->getOldNumberOfNodes()>0){
            nv.push_back(r->getNode(r->getOldNumberOfNodes()-1));
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
        int onon = std::abs(r->getOldNumberOfNodes());
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
    this->getRoots(); // update roots (if necessary)   v.at(0) = Vector3d(0.,0.,0.)
    std::vector<int> nv(this->getNumberOfNewNodes());
    for (const auto& r : roots) {
        int onon = std::abs(r->getOldNumberOfNodes());
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
        int onon = std::abs(r->getOldNumberOfNodes());
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
        int onon = std::abs(r->getOldNumberOfNodes());
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
        int onon = std::abs(r->getOldNumberOfNodes());
        for (size_t i=onon; i<r->getNumberOfNodes(); i++) { // loop over all new nodes
            nodeCTs.at(r->getNodeId(i)-this->old_non) = r->getNodeCT(i); // pray that ids are correct

        }
    }
    return nodeCTs;
}

/**
 * Returns a vector of emergence times of newly created segments since the last call of RootSystem::simulate(dt)
 */
std::vector<double> RootSystem::getNewSegmentsTimes() const
{
    this->getRoots(); // update roots (if necessary)
    std::vector<double> setimes(this->getNumberOfNewNodes());
    int c=0;
    for (const auto& r : roots) {
        int onon = std::abs(r->getOldNumberOfNodes());
        for (size_t i=onon-1; i<r->getNumberOfNodes()-1; i++) {
            setimes.at(c) = r->getNodeCT(i);
            c++;
        }
    }
    return setimes;
}

/**
 *
 */
void RootSystem::push()
{
    stateStack.push(RootSystemState(*this));
}

/**
 *
 */
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
    return Organism::toString();
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
        //writeRSML(fos);
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
 * Writes current simulation results as VTP (VTK polydata file),
 * where each root is represented by a polyline.
 *
 * Use SegmentAnalyser::writeVTP() for a representation based on segments,
 * e.g. for creating a movie (and run the animate.py script), or mapping values to segments
 *
 * todo use tinyxml2, move to Organism
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
    std::string scalarTypeNames[N] = {"type", "order", "radius" };
    for (size_t i=0; i<N; i++) {
        os << "<DataArray type=\"Float32\" Name=\"" << scalarTypeNames[i] <<"\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
        auto scalars = getParameters(scalarTypeNames[i]);
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



/**
 * todo docme
 */
RootSystemState::RootSystemState(const RootSystem& rs) : simtime(rs.simtime), rid(rs.organId), old_non(rs.old_non), old_nor(rs.old_nor),
    numberOfCrowns(rs.numberOfCrowns), gen(rs.gen), UD(rs.UD), ND(rs.ND)
{
    baseRoots = std::vector<RootState>(rs.baseOrgans.size()); // store base roots
    for (size_t i=0; i<baseRoots.size(); i++) {
        baseRoots[i] = RootState(*((Root*)rs.baseOrgans[i]));
    }
}

/**
 * todo docme
 */
void RootSystemState::restore(RootSystem& rs)
{
    rs.roots.clear(); // clear buffer
    // rs.rtparam  = rtparam; // TODO parameters are not restored !!! this might not work right now...
    rs.simtime = simtime; // copy back everything
    rs.organId = rid;
    rs.nodeId = nid;
    rs.old_non = old_non;
    rs.old_nor = old_nor;
    rs.numberOfCrowns = numberOfCrowns;
    rs.gen = gen;
    rs.UD = UD;
    rs.ND = ND;
    for (size_t i=0; i<baseRoots.size(); i++) { // restore base roots
        baseRoots[i].restore(*((Root*)rs.baseOrgans[i]));
    }
}

} // end namespace CRootBox
