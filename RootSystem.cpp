#include "RootSystem.h"

const std::vector<std::string> RootSystem::scalarTypeNames = {"type","radius","order","time","length","surface","1","userdata 1", "userdata 2", "userdata 3", "parent type",
	"basal length", "apical length", "number of branches", "initial growth rate", "insertion angle", "root life time", "mean inter nodal distance", "standard deviation of inter nodal distance", "volume","cumulative uptake",
	 "branch ID"};

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
	for(auto f:gf) {
		delete f;
	}
	for(auto f:tf) {
		delete f;
	}
	baseRoots.clear();
	gf.clear();
	tf.clear();
	simtime=0;
	rid = -1;
	nid = -1;
}

/**
 * Puts default values into the root type parameters vector
 */
void RootSystem::initRTP()
{
	rtparam = std::vector<RootTypeParameter> (maxtypes);
	for (auto& rtp:rtparam) {
		rtp = RootTypeParameter();
	}
}

/**
 * Reads the root parameter from a file. Opens plant parameters with the same filename if available,
 * othterwise assumes a tap root system at position (0,0,-3).
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
	initRTP();
	int c = 0;
	while (cin.good()) {
		RootTypeParameter p;
		p.read(cin);
		setRootTypeParameter(p); // sets the param to the index (p.type-1)
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
	int t = 0;
	for (auto const& rp:rtparam) {
		t++;
		if (rp.type>0) {
			assert(rp.type==t); // check if index is really type-1
			rp.write(os); // only write if defined
		}
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
	//cout << "Root system initialize\n";
	reset(); // just in case

	// Create root system
	const double maxT = 365.; // maximal simulation time
	const double dzB = 0; // distance of basal roots up the mesocotyl [cm] (hardcoded in the orginal version, hardcoded here)
	RootSystemParameter const &rs = rsparam; // rename
	Vector3d iheading(0,0,-1);

	// Taproot
	Root* taproot = new Root(this, 1, iheading ,0, nullptr, 0, 0); // tap root has root type 1
	Vector3d rootCollar(0,0,0);
	taproot->addNode(rootCollar,0);
	taproot->addNode(rs.seedPos,0);
	baseRoots.push_back(taproot);
	//if (rs.maxB==0)
	//{
	//	taproot->addNode(rs.seedPos,0);
	//	baseRoots.push_back(taproot);
	//}

	// Basal roots
	if (rs.maxB>0) {
		if (getRootTypeParameter(basaltype)->type<1) { // if the type is not defined, copy tap root
			std::cout << "Basal root type #" << basaltype << " was not defined, using tap root parameters instead\n";
			RootTypeParameter brtp = RootTypeParameter(*getRootTypeParameter(1));
			brtp.type = basaltype;
			setRootTypeParameter(brtp);
		}
		int maxB = rs.maxB;
		if (rs.delayB>0) {
			maxB = std::min(maxB,int(ceil((maxT-rs.firstB)/rs.delayB))); // maximal for simtime maxT
		}
		double delay = rs.firstB;
		for (int i=0; i<maxB; i++) {
			Root* basalroot = new Root(this, basaltype, iheading ,delay, nullptr, 0, 0);
			Vector3d node = rs.seedPos.plus(Vector3d(0.,0.,dzB));
			basalroot->addNode(node,delay);
			//basalroot->addExistingNode(rs.seedPos,0,0);
			baseRoots.push_back(basalroot);
			//delay += rs.delayB;
			delay = std::log(1-((1-std::exp(rs.growingRateB/rs.maxB*rs.simtime))*(i+1)/rs.maxB))/(rs.growingRateB/rs.maxB); //rs.maxB*(1-std::exp())
		}
	}

	// Shoot borne roots
	if ((rs.nC>0) && (rs.delaySB<maxT)) { // if the type is not defined, copy basal root
		if (getRootTypeParameter(shootbornetype)->type<1) {
			std::cout << "Shootborne root type #" << shootbornetype << " was not defined, using tap root parameters instead\n";
			RootTypeParameter srtp = RootTypeParameter(*getRootTypeParameter(1));
			srtp.type = shootbornetype;
			setRootTypeParameter(srtp);
		}
		Vector3d sbpos = rs.seedPos;
		sbpos.z=sbpos.z/2.; // half way up the mesocotyl
		int maxSB = ceil((maxT-rs.firstSB)/rs.delayRC); // maximal number of root crowns
		double delay = rs.firstSB;
		for (int i=0; i<maxSB; i++) {
			for (int j=0; j<rs.nC; j++) {
				Root* shootborne = new Root(this, shootbornetype, iheading ,delay, nullptr, 0, 0);
				// TODO fix the initial radial heading
				shootborne->addNode(sbpos,delay);
				baseRoots.push_back(shootborne);
				delay += rs.delaySB;
			}
			sbpos.z+=rs.nz;  // move up, for next root crown
			delay = rs.firstSB + i*rs.delayRC; // reset age
		}
	}

	// Create tropisms and growth functions per root type
	for (size_t i=0; i<rtparam.size(); i++) {
		int type = rtparam.at(i).tropismT;
		double N = rtparam.at(i).tropismN;
		double sigma = rtparam.at(i).tropismS;
		TropismFunction* tropism = this->createTropismFunction(type,N,sigma);
		// std::cout << "#" << i << ": type " << type << ", N " << N << ", sigma " << sigma << "\n";
		tf.push_back(new ConfinedTropism(tropism, geometry)); // wrap confinedTropism around baseTropism
		int gft = rtparam.at(i).gf;
		GrowthFunction* gf_ = this->createGrowthFunction(gft);
		gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (ohterwise an exception is thrown)
		gf.push_back(gf_);
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
	for (auto const& r: baseRoots) {
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
 * Sets the seed of the root systems random number generator,
 * and all subclasses using random number generators:
 * @see TropismFunction, @see RootParameter
 *
 * @param seed      random number generator seed
 */
void RootSystem::setSeed(double seed) {
	std::cout << "Setting random seed "<< seed <<"\n";
	gen.seed(1./seed);
	for (auto t : tf) {
		double s  = rand();
		t->setSeed(1./s);
	}
	for (auto rp : rtparam) {
		double s  = rand();
		rp.setSeed(1./s);
	}
}

/**
 * Creates a new lateral root (called by Root:createLateral)
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
	// call Root* lateral = rootsystem->createRoot(lt,  h, delay,  this, length, nodes.size()-1);
	return new Root(this,lt,h,delay,parent,pbl,pni);
}

/**
 * Creates a specific tropsim,
 * the function must be extended or overwritten to add more tropisms
 */
TropismFunction* RootSystem::createTropismFunction(int tt, double N, double sigma) {
	switch (tt) {
	case tt_plagio: return new Plagiotropism(N,sigma);
	case tt_gravi: return new Gravitropism(N,sigma);
	case tt_exo: return new Exotropism(N,sigma);
	case tt_hydro: {
		TropismFunction* gt =  new Gravitropism(N,sigma);
		TropismFunction* ht= new Hydrotropism(N,sigma,soil);
		TropismFunction* cht = new CombinedTropism(N,sigma,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
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
		for (auto const& br:this->baseRoots) {
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
	std::vector<Vector3d> nv = std::vector<Vector3d>(non-rsparam.maxB); // reserve big enough vector
	for (auto const& r: roots) {
		for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
			nv.at(r->getNodeId(i)) = r->getNode(i); // pray that ids are correct
			//std::cout<<i<<" getNodeId "<<r->getNodeId(i)<<" "<<r->getNode(i).x<<"\n";
		}
	}
	//for (auto const& node: nv) {
	//	std::cout<<node.x<<"\n";
	//}
	return nv;
}

/**
 * Copies the nodes of the root systems into a sequential vector,
 * nodes are unique (default). See also RootSystem::getSegments
 */
std::vector<Vector3d> RootSystem::getNodesDGF() const
{
	this->getRoots(); // update roots (if necessary)
	std::cout<<"roots.size() "<<roots.size()<< "\n";
	int non = getNumberOfNodes();
	std::vector<Vector3d> nv = std::vector<Vector3d>(non-rsparam.maxB); // reserve big enough vector
	//std::cout<<"roots.size() "<<roots.size()<< "\n";
	//nv.at(0) = rs.seedPos:
	int No_root = 0;
	for (auto const& r: roots) {
		for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
			//std::cout<<"No_root "<<No_root<< "r->getNodeId(i) "<<r->getNodeId(i)<< "\n";
			if ((i==0) and (r->getNodeId(0)==0))
				nv.at(r->getNodeId(i)) = r->getNode(i); // pray that ids are correct
			else
				if ((i!=0))
					nv.at(r->getNodeId(i)-rsparam.maxB) = r->getNode(i);
		}
		No_root +=1;
	}
	return nv;
}


/**
 * Return the segments of the root system at the current simulation time
 */
std::vector<Vector2i> RootSystem::getSegmentsDGF() const
{
	this->getRoots(); // update roots (if necessary)
	int nos=getNumberOfSegments();
	std::vector<Vector2i> s(nos);
	int c=0;
	for (auto const& r:roots) {
		for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
			Vector2i v;
			if (r->getNodeId(i)<=rsparam.maxB)
			{
				Vector2i v(0,r->getNodeId(i+1)-rsparam.maxB);
				std::cout<<"v.y"<<" "<<v.y<<"\n";
				s.at(c) = v;
			}
			else
			{
				Vector2i v(r->getNodeId(i)-rsparam.maxB,r->getNodeId(i+1)-rsparam.maxB);
				std::cout<<"v.y"<<" "<<v.y<<"\n";
				s.at(c) = v;
			}
			c++;
		}
	}
	return s;
}

/**
 * Return the segments of the root system at the current simulation time
 */
std::vector<double> RootSystem::getLengthFromOrigin() const
{
	this->getRoots(); // update roots (if necessary)
	int nos=getNumberOfSegments();
	auto nodes = RootSystem::getNodes();
	std::vector<double> s(nos);
	int c=0;
	for (auto const& r:roots) {
		double lengthFromOrigin = 0;
		for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
			Vector2i v(r->getNodeId(i),r->getNodeId(i+1));
			Vector3d n1 = nodes.at(v.x);
			Vector3d n2 = nodes.at(v.y);
			double segmentLength = sqrt((n1.x-n2.x)*(n1.x-n2.x)+(n1.y-n2.y)*(n1.y-n2.y)+(n1.z-n2.z)*(n1.z-n2.z));
			lengthFromOrigin += segmentLength;
			//std::cout<<i<<" getNode "<<v.x<<" "<<v.y<<"\n";
			s.at(c) = lengthFromOrigin;
			c++;
		}
	}
	return s;
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
std::vector<Vector2i> RootSystem::getSegments() const
{
	this->getRoots(); // update roots (if necessary)
	int nos=getNumberOfSegments();
	std::vector<Vector2i> s(nos);
	int c=0;
	for (auto const& r:roots) {
		for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
			Vector2i v(r->getNodeId(i),r->getNodeId(i+1));
			//std::cout<<i<<" getNode "<<v.x<<" "<<v.y<<"\n";
			s.at(c) = v;
			c++;
		}
	}
	return s;
}

/**
 * Returns pointers of the roots corresponding to each segment
 */
std::vector<Root*> RootSystem::getSegmentsOrigin() const
{
	this->getRoots(); // update roots (if necessary)
	int nos=getNumberOfSegments();
	std::vector<Root*> s(nos);
	int c=0;
	for (auto const& r:roots) {
		for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
			s.at(c) = r;
			c++;
		}
	}
	return s;
}

/**
 * Copies the node emergence times of the root systems into a sequential vector,
 * see RootSystem::getNodes()
 */
std::vector<double> RootSystem::getNETimes() const
{
	this->getRoots(); // update roots (if necessary)
	int nos=getNumberOfSegments();
	std::vector<double> netv = std::vector<double>(nos); // reserve big enough vector
	int c=0;
	for (auto const& r: roots) {
		for (size_t i=1; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
			netv.at(c) = r->getNodeETime(i); // pray that ids are correct
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
			rt[i] = roots[j]->getNodeETime(i);
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
		double value=0;
		switch(stype) {
		case st_type:  // type
			value = roots[i]->param.type;
			break;
		case st_radius: // root radius
			value = roots[i]->param.a;
			break;
		case st_order: { // root order (calculate)
			value = 0;
			Root* r_ = roots[i];
			while (r_->parent!=nullptr) {
				value++;
				r_=r_->parent;
			}
			break;
		}
		case st_rootID: { // root order (calculate)
			value = roots[i]->id;
			break;
		}
		case st_time: // emergence time of the root
			value = roots[i]->getNodeETime(0);
			break;
		case st_length:
			value = roots[i]->length;
			break;
		case st_surface:
			value =  roots[i]->length*2.*M_PI*roots[i]->param.a;
			break;
		case st_volume:
			value =  roots[i]->length*M_PI*roots[i]->param.a*roots[i]->param.a;
			break;
		case st_one:
			value =  1;
			break;
		case st_parenttype: {
			Root* r_ = roots[i];
			if (r_->parent!=nullptr) {
				value = r_->parent->param.type;
			} else {
				value = 0;
			}
			break;
		}
		case st_lb:
			value = roots[i]->param.lb;
			break;
		case st_la:
			value = roots[i]->param.la;
			break;
		case st_nob:
			value = roots[i]->param.nob;
			break;
		case st_r:
			value = roots[i]->param.r;
			break;
		case st_theta:
			value = roots[i]->param.theta;
			break;
		case st_rlt:
			value = roots[i]->param.rlt;
			break;
		case st_meanln: {
			const std::vector<double>& v_ = roots[i]->param.ln;
			value = std::accumulate(v_.begin(), v_.end(), 0.0) / v_.size();
			break;
		}
		case st_stdln: {
			const std::vector<double>& v_ = roots[i]->param.ln;
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
 * todo test comment
 */
std::vector<int> RootSystem::getNodeUpdateIndices()
{
	this->getRoots(); // update roots (if necessary)
	std::vector<int> ni = std::vector<int> (old_nor);
	int c =0;
	for (auto const& r: roots) {
		if (r->old_non>0){
			ni.at(r->getNodeId(r->old_non));
			c++;
		}
	}
	return ni;
}

/**
 *
 */
std::vector<Vector3d> RootSystem::getUpdatedNodes()
{
	this->getRoots(); // update roots (if necessary)
	int non = getNumberOfNodes();
	std::vector<Vector3d> un(non-old_non); // reserve big enough vector

	return un;
}

/**
 * Returns a vector of newly created nodes since the last call of RootSystem::simulate(dt),
 * to dynamically add to the old node vector, see also RootSystem::getNodes
 */
std::vector<Vector3d> RootSystem::getNewNodes()
{
	this->getRoots(); // update roots (if necessary)
	int non = getNumberOfNodes();
	std::vector<Vector3d> nv(non-old_non); // reserve big enough vector
	for (auto const& r: roots) {
		for (size_t i=r->old_non+1; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
			nv.at(r->getNodeId(i)-old_non) = r->getNode(i); // pray that ids are correct
		}
	}
	return nv;
}

/**
 * Returns a vector of newly created segments since the last call of RootSystem::simulate(dt),
 * to dynamically add to the old segment index vector, see also RootSystem::getSegments
 */
std::vector<Vector2i> RootSystem::getNewSegments()
{
	this->getRoots(); // update roots (if necessary)
	int non=getNumberOfNodes()-1;
	std::vector<Vector2i> si(non);
	int c=0;
	for (auto const& r:roots) {
		for (size_t i=r->old_non; i<r->getNumberOfNodes()-1; i++) {
			Vector2i v(r->getNodeId(i),r->getNodeId(i+1));
			si.at(c) = v;
			c++;
		}
	}
	return si;
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
	} else if (ext.compare("dgf")==0)  {
		std::cout << "writing DGF ... "<< name.c_str() <<"\n";
		writeDGF(fos);
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
	for (auto const& root :baseRoots) {
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
	for (auto const& r : roots) {
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
	const size_t N = 4; // SCALARS
	int types[N] = { st_type, st_order, st_radius, st_rootID };
	for (size_t i=0; i<N; i++) {
		os << "<DataArray type=\"Float32\" Name=\"" << scalarTypeNames[types[i]] <<"\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
		auto scalars = getScalar(types[i]);
		for (auto s : scalars) {
			os << s<< " ";
		}
		os << "\n</DataArray>\n";
	}
	os << "\n</CellData>\n";

	// POINTS (=nodes)
	os << "<Points>\n"<<"<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
	for (auto const& r:nodes) {
		for (auto const& n : r) {
			os << n.x << " "<< n.y <<" "<< n.z<< " ";
		}
	}
	os << "\n</DataArray>\n"<< "</Points>\n";

	// LINES (polylines)
	os << "<Lines>\n"<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	int c=0;
	for (auto const& r:roots) {
		for (size_t i=0; i<r->getNumberOfNodes(); i++) {
			os << c << " ";
			c++;
		}
	}
	os << "\n</DataArray>\n"<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	c = 0;
	for (auto const& r:roots) {
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

void RootSystem::writeDGF(std::ostream & os) const
{
	auto nodes = RootSystem::getNodes();
	auto segments = RootSystem::getSegments();
	auto segO = RootSystem::getSegmentsOrigin();
	auto ctimes = RootSystem::getNETimes();
	os << "DGF \n";
	os << "Vertex \n";
	for (auto& n : nodes) {
		os << n.x/100 << " " << n.y/100 << " " << n.z/100 << " \n";
		//std::cout << n.x/100 << " " << n.y/100 << " " << n.z/100 << " \n";
	}

	os << "# \n";
	os << "SIMPLEX \n";
	os << "parameters 5 \n";
	// node1ID, node2ID, type, branchID, surfaceIdx,  massIdx, creationTimeId
	for (size_t i=0; i<segments.size(); i++) {
		Vector2i s = segments.at(i);
		//std::cout<<"s "<< s.x <<" "<<s.y;
		Vector3d n1 = nodes.at(s.x);
		Vector3d n2 = nodes.at(s.y);
		Root* r = segO.at(i);
		int branchnumber = r->id;
		double radius = r->param.a;
		double length = sqrt((n1.x-n2.x)*(n1.x-n2.x)+(n1.y-n2.y)*(n1.y-n2.y)+(n1.z-n2.z)*(n1.z-n2.z));
		double surface = 2*radius*M_PI*length;
		double time = ctimes.at(i);
		double type = r->param.type;
		os << s.x << " " << s.y << " " << type << " " << branchnumber << " " << surface/10000 <<" 0.00 " << time*3600*24 << "\n";
		//std::cout << s.x << " " << s.y << " " << type << " " << branchnumber << " " << surface/10000 << " " << length/100 <<" " << radius/100 << " " << "0.00" << " " << "0.0001" << " "<< "0.00001" << " " << time*3600*24 << " \n";

	}

	os << "# \n";
	os << "BOUNDARYDOMAIN \n";
	os << "default 1 \n";
	os << "# \n";
}