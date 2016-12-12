#include "RootSystem.h"

const std::vector<std::string> RootSystem::scalarTypeNames = {"type","radius","order","red","green","blue","time","length","surface"};

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
 * a confining geometry, the tropsim functions, and the growth functions.
 *
 * Call this method before simulation and after setting geometry, plant and root parameters
 */
void RootSystem::initialize(int basaltype, int shootbornetype)
{
    //cout << "Root system initialize\n";
    reset(); // just in case

    // Create root system
    const double maxT = 365.; // maximal simulation time
    const double dzB = 0.1; // distance of basal roots up the mesocotyl [cm] (hardcoded in the orginal version, hardcoded here)
    RootSystemParameter const &rs = rsparam; // rename
    Vector3d iheading(0,0,-1);

    // Taproot
    Root* taproot = new Root(this, 1, iheading ,0, nullptr, 0, 0); // tap root has root type 1
    taproot->addNode(rs.seedPos,0);
    baseRoots.push_back(taproot);

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
            Vector3d node = rs.seedPos.minus(Vector3d(0.,0.,dzB));
            basalroot->addNode(node,delay);
            baseRoots.push_back(basalroot);
            delay += rs.delayB;
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
        TropismFunction* tropism = createTropismFunction(type,N,sigma);
        // std::cout << "#" << i << ": type " << type << ", N " << N << ", sigma " << sigma << "\n";
        tf.push_back(new ConfinedTropism(tropism, geometry)); // wrap confinedTropism around baseTropism
        int gft = rtparam.at(i).gf;
        GrowthFunction* gf_ = createGrowthFunction(gft);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (ohterwise an exception is thrown)
        gf.push_back(gf_);
    }

}

/**
 * Simulates root system growth for time span dt
 *
 * @param dt    time step [days]
 */
void RootSystem::simulate(double dt)
{
    std::cout << "RootSystem.simulate(dt) from "<< simtime << " to " << simtime+dt << " days \n";
    simtime+=dt;
    for (auto const& r: baseRoots) {
        r->simulate(dt);
    }
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
    gen.seed(seed);
    for (auto t : tf) {
        double s  = rand();
        t->setSeed(s);
    }
    for (auto rp : rtparam) {
        double s  = rand();
        rp.setSeed(s);
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
 * copies the root only, if it has more than 1 node
 *
 * \return sequential vector of roots with more than 1 node
 */
std::vector<Root*> RootSystem::getRoots() const
{
    std::vector<Root*> v = std::vector<Root*>();
    for (auto const& br:this->baseRoots) {
        br->getRoots(v);
    }
    return v;
}

/**
 * Returns the positions of the root tips
 *
 * @param roots		a vector of roots, if no roots are specified all roots are returned (@see RootSystem::getRoots)
 */
std::vector<Vector3d> RootSystem::getRootTips(std::vector<Root*> roots) const
{
	if (roots.empty()) {
		roots = this->getRoots();
	}
	std::vector<Vector3d> tips;
	for (auto& r : roots) {
		tips.push_back(r->getNode(r->getNumberOfNodes()-1));
	}
	return tips;
}

/**
 * Returns the positions of the root bases
 *
 * @param roots         a vector of roots, if no roots are specified all roots are returned (@see RootSystem::getRoots)
 */
std::vector<Vector3d> RootSystem::getRootBases(std::vector<Root*> roots) const
{
        if (roots.empty()) {
                roots = this->getRoots();
        }
	std::vector<Vector3d> bases;
	for (auto& r : roots) {
	        bases.push_back(r->getNode(0));
	}
	return bases;
}

/**
 * Copies the nodes of the root systems into a sequential vector,
 * there are two different node numberings
 *
 * @param ot        ot_segments: each segment is a line, nodes are unique (default)
 *                  ot_polyline: each root is a line, nodes are not unique
 * @param roots     a vector of roots, if no roots are specified all roots are returned (@see RootSystem::getRoots)
 */
std::vector<Vector3d> RootSystem::getNodes(int ot, std::vector<Root*> roots) const
{
	if (roots.empty()) {
		roots = this->getRoots();
	}
    switch (ot) {
    case ot_segments: {
        int non = getNumberOfNodes();
        std::vector<Vector3d> nv = std::vector<Vector3d>(non); // reserve big enough vector
        for (auto const& r: roots) {
            for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
                nv.at(r->getNodeId(i)) = r->getNode(i); // pray that ids are correct
            }
        }
        return nv;
    }
    case ot_polylines: {
        int non = 0;
        for (auto const& r : roots) {
            non += r->getNumberOfNodes();
        }
        std::vector<Vector3d> nv = std::vector<Vector3d>(non); // reserve big enough vector
        int c=0;
        for (auto const& r: roots) {
            for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
                nv.at(c) = r->getNode(i);
                c++;
            }
        }
        return nv;
    }
    default: throw std::invalid_argument( "RootSystem::getNodes() output type not implemented" );
    }
}

/**
 * Return the segments of the root system at the current simulation time
 *
 * @param ot        ot_segments: each segment connects two unique nodes
 *                  ot_polyline: not implemented
 * @param roots     a vector of roots, if no roots are specified all roots are returned (@see RootSystem::getRoots)
 */
std::vector<Vector2i> RootSystem::getSegments(int ot, std::vector<Root*> roots) const
{
	if (roots.empty()) {
		roots = this->getRoots();
	}
    switch (ot) {
    case ot_segments: {
        int non=0;
        for (auto const& r:roots) {
            non += r->getNumberOfNodes()-1;
        }
        std::vector<Vector2i> s(non);
        int c=0;
        for (auto const& r:roots) {
            for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
                Vector2i v(r->getNodeId(i),r->getNodeId(i+1));
                s.at(c) = v;
                c++;
            }
        }
        return s;
    }
    default: throw std::invalid_argument( "RootSystem::getSegments() output type not implemented" );
    }
}

/**
 * Returns pointers of the roots corresponding to each segment
 * @param ot        ot_segments: each segment connects two unique nodes
 *                  ot_polyline: not implemented
 * @param roots     a vector of roots, if no roots are specified all roots are returned (@see RootSystem::getRoots)
 */
std::vector<Root*> RootSystem::getSegmentsOrigin(int ot, std::vector<Root*> roots) const
{
	if (roots.empty()) {
		roots = this->getRoots();
	}
    switch (ot) {
    case ot_segments: {
        int non=0;
        for (auto const& r:roots) {
            non += r->getNumberOfNodes()-1;
        }
        std::vector<Root*> s(non);
        int c=0;
        for (auto const& r:roots) {
            for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
                s.at(c) = r;
                c++;
            }
        }
        return s;
    }
    default: throw std::invalid_argument( "RootSystem::getSegmentsOrigin() output type not implemented" );
    }
}

/**
 * Copies the node emergence times of the root systems into a sequential vector,
 * there are two different node numberings
 *
 * Gives either cell data (in case of segments) or point data (in case of polyline)
 *
 * @param ot        out_segments: each segment is a line, nodes are unique
 *                  out_polyline: each root is a line (and cell), nodes are not unique
 * @param roots     a vector of roots, if no roots are specified all roots are returned (@see RootSystem::getRoots)
 */
std::vector<double> RootSystem::getNETimes(int ot, std::vector<Root*> roots) const
{
	if (roots.empty()) {
		roots = this->getRoots();
	}
    switch (ot) {
    case ot_segments: {
        int non=0;
        for (auto const& r:roots) {
            non += r->getNumberOfNodes()-1;
        }
        std::vector<double> netv = std::vector<double>(non); // reserve big enough vector
        int c=0;
        for (auto const& r: roots) {
            for (size_t i=1; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
                netv.at(c) = r->getNodeETime(i); // pray that ids are correct
                c++;
            }
        }
        return netv;
    }
    case ot_polylines: {
        int non = 0;
        for (auto const& r : roots) {
            non += r->getNumberOfNodes();
        }
        std::vector<double> netv = std::vector<double>(non); // reserve big enough vector
        int c=0;
        for (auto const& r: roots) {
            for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
                netv.at(c) = r->getNodeETime(i);
                c++;
            }
        }
        return netv;
    }
    default: throw std::invalid_argument( "RootSystem::getNETimes() output type not implemented" );
    }
}

/**
 * Copies a scalar that is constant per root to a sequential vector (data per cell).
 * st_time is the creation time of the root
 *
 * @param ot        out_segments: scalar for each segment
 *                  out_polyline: scalar for each root
 * @param stype     a scalar type: type, radius, order, red, green, blue,... (@see RootSystem::ScalarTypes)
 * @param roots     a vector of roots, if no roots are specified all roots are returned (@see RootSystem::getRoots)
 */
std::vector<double> RootSystem::getScalar(int ot, int stype, std::vector<Root*> roots) const
{
    if (roots.empty()) {
        roots = this->getRoots();
    }
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
        case st_red:
            value = rtparam.at(roots[i]->param.type-1).colorR;
            break;
        case st_green:
            value = rtparam.at(roots[i]->param.type-1).colorG;
            break;
        case st_blue:
            value = rtparam.at(roots[i]->param.type-1).colorB;
            break;
        case RootSystem::st_time:
            value = roots[i]->getNodeETime(0);
            break;
        case RootSystem::st_length:
            value = roots[i]->length;
            break;
        case RootSystem::st_surface:
            value =  roots[i]->length*2.*M_PI*roots[i]->param.a;
            break;
        default:
            throw std::invalid_argument( "RootSystem::copyRootParam() type not implemented" );
        }
        scalars[i]=value;
    }
    if (ot==ot_polylines) {
        return scalars;
    } else if (ot==ot_segments) {
        std::vector<double> values(getNumberOfNodes()-baseRoots.size());
        int c = 0;
        for (size_t i=0; i<roots.size(); i++) {
            for (size_t j=1; j<roots[i]->getNumberOfNodes(); j++) { // loop over all nodes of the root
                values[c] = scalars[i];
                c++;
            }
        }
        return values;
    } else {
        throw std::invalid_argument( "RootSystem::getScalar() type not implemented" );
    }
}

/**
 * Exports the simulation results with the type from the extension in name
 * (that must be lower case)
 *
 * @param name      file name e.g. output.vtp
 */
void RootSystem::write(std::string name, int type) const
{
    std::ofstream fos;
    fos.open(name.c_str());
    std::string ext = name.substr(name.size()-3,name.size()); // pick the right writer
    if (ext.compare("sml")==0) {
        std::cout << "writing RSML... "<< name.c_str() <<"\n";
        writeRSML(fos);
    } else if (ext.compare("vtp")==0) {
        std::cout << "writing VTP... "<< name.c_str() <<"\n";
        writeVTP(fos,type);
    } else if (ext.compare(".py")==0)  {
        std::cout << "writing Geometry ... "<< name.c_str() <<"\n";
        writeGeometry(fos);
    } else if (ext.compare("dgf")==0) {
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
    std::cout << "writing metadata \n";
    writeRSMLMeta(os);
    std::cout << "writing plant \n";
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
    os << "<version>" << 1 << "</version>\n";
    os << "<unit>" << "cm" << "</unit>\n";
    os << "<resolution>" << 1 << "</resolution>\n";
    // fetch time
    //    os << "<last-modified>";
    //    auto t = std::time(nullptr);
    //    auto tm = *std::localtime(&t);
    //    os << std::put_time(&tm, "%d-%m-%Y"); // %H-%M-%S" would do the job for gcc 5.0
    //    os << "</last-modified>\n";
    os << "<software>CRootBox</software>\n";
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
 * use vtp_segments for creating a movie (and run the animate.py script)
 *
 * @param os      typically a file out stream
 * @param type    vtp_segments: each segment is a line, vtp_polyline: each root is a line
 */
void RootSystem::writeVTP(std::ostream & os, int type) const
{
    auto roots = getRoots();
    auto nodes = getNodes(type,roots);
    auto times = getNETimes(type,roots);

    os << "<?xml version=\"1.0\"?>";
    os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    os << "<PolyData>\n";
    int non=nodes.size(); // number of nodes
    int nol=0; // number of lines
    switch (type) {
    case ot_segments:
        nol = non-baseRoots.size(); break; // each tree has number of node -1 segments (the base roots are not connected)
    case ot_polylines:
        nol = roots.size();  break; // number of lines = number of roots
    }
    os << "<Piece NumberOfLines=\""<< nol << "\" NumberOfPoints=\""<<non<<"\">\n";

    // POINTDATA
    if (type==ot_polylines) {
        os << "<PointData Scalars=\" PointData\">\n" << "<DataArray type=\"Float32\" Name=\"time\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
        for (auto t : times) {
            os << t << " ";
        }
        os << "\n</DataArray>\n" << "\n</PointData>\n";
    }
    // CELLDATA (live on the segments or polylines)
    os << "<CellData Scalars=\" CellData\">\n";
    if (type==ot_segments) { // TIME
        os << "<DataArray type=\"Float32\" Name=\"time\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
        for (auto t : times) {
            os << t << " ";
        }
        os << "\n</DataArray>\n";
    }
    const size_t N = 3; // SCALARS
    int types[N] = { st_type, st_order, st_radius };
    for (size_t i=0; i<N; i++) {
        os << "<DataArray type=\"Float32\" Name=\"" << scalarTypeNames[types[i]] <<"\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
        auto scalars = getScalar(type,types[i],roots);
        for (auto s : scalars) {
            os << s<< " ";
        }
        os << "\n</DataArray>\n";
    }
    os << "\n</CellData>\n";

    // POINTS (=nodes)
    os << "<Points>\n"<<"<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
    for (auto const& n:nodes) {
        os << n.x << " "<< n.y <<" "<< n.z<< " ";
    }
    os << "\n</DataArray>\n"<< "</Points>\n";

    // LINES (=segments or polylines)
    os << "<Lines>\n"<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
    switch (type) {
    case ot_segments:
        for (auto const& r:roots) {
            for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
                os << r->getNodeId(i) << " " << r->getNodeId(i+1) << " ";
            }
        }
        os << "\n</DataArray>\n"<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
        for (int i=0; i<(nol); i++) {
            os << 2*i+2 << " ";
        }
        break;
    case ot_polylines:
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
        break;
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
 * TODO finish and test
 */
void RootSystem::writeDGF(std::ostream & os) const
{
    os << "DGF \n";
    os << "Vertex \n";
    //os << "x1 y1 z1 x2 y2 z2 radius R G B time type \n";

    auto roots = getRoots();
    std::vector<int> x1;
    std::vector<int> x2;
    std::vector<int> branchnumber;
    int k=0;
    for (auto const& r:roots) {
        k=k+1;
        for (size_t i=0; i<r->getNumberOfNodes()-1; i++) {
            x1.push_back(r->getNodeId(i)); // from x1
            x2.push_back(r->getNodeId(i+1)); // to x2
            branchnumber.push_back(k);
        }
    }
    auto radius_ = getScalar(ot_segments,st_radius,roots);
    auto time_ = getNETimes(ot_segments,roots);
    auto type_ = getScalar(ot_segments,st_type,roots);
    auto nodes = getNodes(ot_segments,roots);

    for (size_t i =0; i<x1.size(); i++) {
        Vector3d n1 = nodes.at(x1[i]);
        os << n1.x << " " << n1.y << " " << n1.z << " \n";
    }
    Vector3d n2 = nodes.at(x2[x2.size()-1]);
    os << n2.x << " " << n2.y << " " << n2.z << " \n";

    os << "# \n";
    os << "SIMPLEX \n";
    os << "parameters 4 \n";

    //    for (size_t i =0; i<x1.size(); i++) {
    //        Vector3d n1 = nodes.at(x1[i]);
    //        Vector3d n2 = nodes.at(x2[i]);
    //        os << x1[i] << " " << x2[i] << " " << type_.at(i) <<              " "<< branchnumber[i] <<" " << radius_.at(i) << " " << "0.00" << " " << "0.00"<< " \n";
    //    }

    os << "# \n";
    os << "BOUNDARYDOMAIN \n";
    os << "default 1 \n";
    os << "# \n";

}

