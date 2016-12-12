#include "analysis.h"
#include <iostream>
#include <iomanip>

/**
 * Copies the segments of the roots system into the analysis class
 *
 * @param rs      the root system that is analysed
 */
AnalysisSDF::AnalysisSDF(const RootSystem& rs)
{
	auto roots = rs.getRoots();
	nodes = rs.getNodes(RootSystem::ot_segments,roots);
	segments = rs.getSegments(RootSystem::ot_segments,roots);
	ctimes = rs.getNETimes(RootSystem::ot_segments,roots);
	segO = rs.getSegmentsOrigin(RootSystem::ot_segments,roots);
	assert(segments.size()==ctimes.size());
	assert(segments.size()==segO.size());
}

/**
 * Adds all segmentes from root system @param rs to the analysis.
 */
void AnalysisSDF::addSegments(const RootSystem& rs)
{
	addSegments(AnalysisSDF(rs));
}

/**
 * Adds all segmentes from the analyser @param a to this analysis.
 */
void AnalysisSDF::addSegments(const AnalysisSDF& a)
{
	int offset = nodes.size();
	nodes.insert(nodes.end(),a.nodes.begin(),a.nodes.end()); // copy nodes
	auto ns = a.segments;
	for (auto& s : ns) { // shift indices
		s.x += offset;
		s.y += offset;
	}
	segments.insert(segments.end(),ns.begin(),ns.end()); // copy segments
	ctimes.insert(ctimes.end(),a.ctimes.begin(),a.ctimes.end()); // copy times
	segO.insert(segO.end(),a.segO.begin(),a.segO.end());// copy origins
	assert(segments.size()==ctimes.size());
	assert(segments.size()==segO.size());
}

/**
 * Returns a specific parameter per root segment
 *
 * @param st    parameter type @see RootSystem::ScalarType
 * \return      vector containing parameter value per segment
 */
std::vector<double> AnalysisSDF::getScalar(int st) const
{
	// TODO since segments have know their origin root: segO, pure root segment mapping could be done by using RootSystem::getScalar
	// to avoid redundant code,
	// scalars that are segment wise must be treated extra
	std::vector<double> data(segO.size());
	double v = 0; // value
	for (size_t i=0; i<segO.size(); i++) {
		const auto& r = segO.at(i);
		switch (st) {
		case RootSystem::st_type:
			v=r->param.type;
			break;
		case RootSystem::st_radius:
			v=r->param.a;
			break;
		case RootSystem::st_order: {
			Root* r_ = r;
			while (r_->parent!=nullptr) {
				v++;
				r_=r_->parent;
			}
		}
		break;
		case RootSystem::st_time: // at the end of the method
			break;
		case RootSystem::st_length: {
			Vector2i s = segments.at(i);
			Vector3d x = nodes.at(s.x);
			Vector3d y = nodes.at(s.y);
			v = x.minus(y).length();
		}
		break;
		case RootSystem::st_surface: {
			Vector2i s = segments.at(i);
			Vector3d x = nodes.at(s.x);
			Vector3d y = nodes.at(s.y);
			v = x.minus(y).length()*2*M_PI*r->param.a;
		}
		break;
		default:
			throw std::invalid_argument( "AnalysisSDF::getData() type not implemented" );
		}
		data.at(i) = v;
	}
	if (st==RootSystem::st_time) {
		data = ctimes;
	}
	return data;
}

/**
 * Crops the segments with some geometry
 *
 * @param geometry      signed distance function of the geometry
 */
void AnalysisSDF::crop(SignedDistanceFunction* geometry)
{
	//std::cout << "cropping " << segments.size() << " segments...";
	std::vector<Vector2i> seg;
	std::vector<Root*> sO;
	std::vector<double> ntimes;
	for (size_t i=0; i<segments.size(); i++) {
		auto s = segments.at(i);
		Vector3d x = nodes.at(s.x);
		Vector3d y = nodes.at(s.y);
		bool x_ = geometry->getDist(x)<=0; // in?
		bool y_ = geometry->getDist(y)<=0; // in?
		if ((x_==true) && (y_==true)) { //segment is inside
			seg.push_back(s);
			sO.push_back(segO.at(i));
			ntimes.push_back(ctimes.at(i));
		} else if ((x_==false) && (y_==false)) { // segment is outside

		} else { // one node is inside, one outside
			// sort
			Vector3d in;
			Vector3d out;
			int ini;
			if (x_==true) {
				in = x;
				ini = s.x;
				out = y;
			} else {
				in = y;
				ini = s.y;
				out = x;
			}
			// cut
			Vector3d newnode = cut(in, out, geometry);
			// add new segment
			nodes.push_back(newnode);
			Vector2i newseg(ini,nodes.size()-1);
			seg.push_back(newseg);
			sO.push_back(segO.at(i));
			ntimes.push_back(ctimes.at(i));
		}

	}
	segments = seg;
	segO  = sO;
	ctimes = ntimes;
	std::cout << " cropped to " << segments.size() << " segments " << "\n";

}

/**
 * Filters the segments to the ones, where data is within [min,max], @see AnalysisSDF::getData,
 * i.e. all other segments are deleted.
 *
 * @param st    parameter type @see RootSystem::ScalarType
 * @param min   minimal value
 * @param max   maximal value
 */
void AnalysisSDF::filter(int st, double min, double max)
{
  std::vector<double> data = getScalar(st);
  std::vector<Vector2i> seg;
  std::vector<Root*> sO;
  std::vector<double> ntimes;
  for (size_t i=0; i<segments.size(); i++) {
      if ((data.at(i)>=min) && (data.at(i)<=max)) {
          seg.push_back(segments.at(i));
          sO.push_back(segO.at(i));
          ntimes.push_back(ctimes.at(i));
      }
  }
  segments = seg;
  segO  = sO;
  ctimes = ntimes;
}

/**
 * Filters the segments to the ones, where data eqals value, @see AnalysisSDF::getData,
 * i.e. all other segments are deleted.
 *
 * @param st        parameter type @see RootSystem::ScalarType
 * @param value     parameter value of the segments that are kept
 */
void AnalysisSDF::filter(int st, double value)
{
	std::vector<double> data = getScalar(st);
	std::vector<Vector2i> seg;
	std::vector<Root*> sO;
	std::vector<double> ntimes;
	for (size_t i=0; i<segments.size(); i++) {
		if (data.at(i)==value) {
			seg.push_back(segments.at(i));
			sO.push_back(segO.at(i));
			ntimes.push_back(ctimes.at(i));
		}
	}
	segments = seg;
	segO  = sO;
	ctimes = ntimes;
}

/**
 * Sorts nodes and deletes unused nodes.
 * This can save a lot of memory, since AnalysisSDF::crop and AnalysisSDF::filter only delete segments, not unused nodes
 */
void AnalysisSDF::pack() {
	std::vector<double> ni(nodes.size());
	std::fill(ni.begin(),ni.end(), 0.);
	std::vector<Vector3d> newnodes;
	for (auto& s:segments) {
		if (ni.at(s.x)==0) { // the node is new
			newnodes.push_back(nodes.at(s.x));
			ni.at(s.x) = newnodes.size()-1; // set index of the new node
		}
		s.x = ni.at(s.x);
		if (ni.at(s.y)==0) { // the node is new
			newnodes.push_back(nodes.at(s.y));
			ni.at(s.y) = newnodes.size()-1; // set index of the new node
		}
		s.y = ni.at(s.y);
	}
	// std::cout << "pack(): nodes: " << nodes.size() << " -> " << newnodes.size() << ", " << double(newnodes.size())/double(nodes.size()) << " \n";
	nodes = newnodes; // kabum!
}

/**
 *  Numerically computes the intersetion point
 *
 * @param in       the node within the domain
 * @param out      the node outside of the domain
 * @param geometry signed distance function of the geometry
 * \return         the intersection point
 */
Vector3d AnalysisSDF::cut(Vector3d in, Vector3d out, SignedDistanceFunction* geometry)
{
	assert(geometry->getDist(in)<=0);
	assert(geometry->getDist(out)>=0);
	if (std::abs(geometry->getDist(out))>1e-6) {
		Vector3d c =  in.plus(out).times(0.5); // mid
		if (geometry->getDist(c)<0) { // in
			return cut(c,out,geometry);
		} else { // out
			return cut(in,c,geometry);
		}
	} else {
		return out;
	}
}

/**
 * \return The summed parameter of type @param st (@see RootSystem::ScalarType)
 */
double AnalysisSDF::getSummed(int st) const {
	std::vector<double> data = getScalar(st);
	double v = 0;
	for (auto const& d : data) {
		v += d;
	}
	return v;
}

/**
 * \return The summed parameter of type @param st (@see RootSystem::ScalarType), that is within geometry @param g ,
 * based on the segment mid point (i.e. not exact)
 */
double AnalysisSDF::getSummed(int st, SignedDistanceFunction* g) const {
	std::vector<double> data = getScalar(st);
	double v = 0;
	for (size_t i=0; i<segments.size(); i++) {
		double d = data.at(i);
		Vector2i s = segments.at(i);
		Vector3d n1 = nodes.at(s.x);
		Vector3d n2 = nodes.at(s.y);
		Vector3d mid = n1.plus(n2).times(0.5);
		if (g->getDist(mid)<0) {
			v += d;
		}
	}
	return v;
}

/**
 * \return The number of roots
 */
int AnalysisSDF::getNumberOfRoots() const
{
	std::set<Root*> rootset;  // praise the stl
	for (Root* r : segO) {
		rootset.insert(r);
	}
	return rootset.size();
}

/**
 * Projects the segments to an image plane
 *
 * @param pos       position of camera
 * @parma ons       orthonormal system, row 1 is orthogonal to the image plane given by [row 2,row 3]
 * @param fl        focal length, alpha = 2*arctan(d/(2*fl)), were alpha is the angle of field, and d the image diagonal
 *
 * \return The image segments in the x-y plane (z=0)
 */
AnalysisSDF AnalysisSDF::foto(const Vector3d& pos, const Matrix3d& ons, double fl) const
{
	AnalysisSDF f(*this); // copy
	for (auto& n : f.nodes) { // translate
		n = n.minus(pos);
	}
	Matrix3d m = ons.inverse(); // rotate
	for (auto& n : f.nodes) {
		n = m.times(n);
	}
//	// crop to objects in front of the camera
	Vector3d o(0.,0.,0.);
	Vector3d plane(0.,0., -1);
	SDF_HalfPlane sdf = SDF_HalfPlane(o,plane);
	f.crop(&sdf);
	f.pack();
	// project
	for (auto& a : f.nodes) {
		a = a.times(fl/(-plane.times(a)));
		a.z = 0;
	}
	// final image crop
	SDF_PlantBox box(1.,1.,2.); // TODO --> d = sqrt(2) cm
	f.crop(&box);
	f.pack();
	return f;
}

/**
 * Cuts the segments with a plane
 */
AnalysisSDF AnalysisSDF::cut(const SDF_HalfPlane& plane) const
{
	AnalysisSDF f;
	f.nodes = nodes; // copy all nodes
	for (size_t i=0; i<segments.size(); i++) {
		Vector2i s = segments.at(i);
		Vector3d n1 = nodes.at(s.x);
		Vector3d n2 = nodes.at(s.y);
		double d = plane.getDist(n1)*plane.getDist(n2);
		if (d<=0) { // one is inside, one is outside
			f.segments.push_back(s);
			f.ctimes.push_back(ctimes.at(i));
			f.segO.push_back(segO.at(i));
		}
	}
	f.pack(); // delete unused nodes
	return f;
}

/**
 *  Creates a vertical distribution of the parameter of type @param st (@see RootSystem::ScalarType)
 *
 * @param st        parameter type @see RootSystem::ScalarType
 * @param top       vertical top position (cm)
 * @param bot       vertical bot position (cm)
 * @param n         number of layers (each with a height of (bot-top)/n )
 * @param exact     calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)
 * \return Vector of size @param n containing the summed parameter in this layer
 */
std::vector<double> AnalysisSDF::distribution(int st, double top, double bot, int n, bool exact) const
{
	std::vector<double> d(n);
	double dz = (bot-top)/double(n);
	SDF_PlantBox* layer = new SDF_PlantBox(1e100,1e100,dz);
	for (int i=0; i<n; i++) {
		Vector3d t(0,0,-dz/2.-i*dz);
		SDF_RotateTranslate g(layer,t);
		if (exact) {
			AnalysisSDF a(*this); // copy everything
			a.crop(&g); // crop exactly
			d.at(i) = a.getSummed(st);
		} else {
			d.at(i) = this->getSummed(st, &g);
		}
	}
	delete layer;
	return d;
}

/**
 *  Creates a vertical distribution
 *
 * @param st        parameter type @see RootSystem::ScalarType
 * @param top       vertical top position (cm)
 * @param bot       vertical bot position (cm)
 * @param n         number of layers (each with a height of (bot-top)/n )
 * @param exact     calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)
 * \return Vector of size @param n containing an Analysis object of the layer (cropped exactly)
 */
std::vector<AnalysisSDF> AnalysisSDF::distribution(double top, double bot, int n) const
{
	std::vector<AnalysisSDF> d(n);
	double dz = (bot-top)/double(n);
	SDF_PlantBox* layer = new SDF_PlantBox(1e100,1e100,dz);
	for (int i=0; i<n; i++) {
		Vector3d t(0,0,-dz/2.-i*dz);
		SDF_RotateTranslate g(layer,t);
		AnalysisSDF a = AnalysisSDF(*this); // copy everything
		a.crop(&g); // crop exactly
		d.at(i) = a;
	}
	delete layer;
	return d;
}

/**
 *  Creates a two-dimensional distribution of the parameter of type @param st (@see RootSystem::ScalarType)
 *
 * @param st        parameter type @see RootSystem::ScalarType
 * @param top       vertical top position (cm)
 * @param bot       vertical bot position (cm)
 * @param left      left along x-axis (cm)
 * @param rigth     right along x-axis (cm)
 * @param n         number of vertical grid elements (each with height of (bot-top)/n )
 * @param m 		number of horizontal grid elements (each with length of (right-left)/m)
 * @param exact     calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)
 * \return Vector of size @param n containing the summed parameter in this layer
 */
std::vector<std::vector<double>> AnalysisSDF::distribution2(int st, double top, double bot, double left, double right, int n, int m, bool exact) const
{
	std::vector<std::vector<double>> d(n);
	double dz = (bot-top)/double(n);
	double dx = (right-left)/double(m);
	SDF_PlantBox* layer = new SDF_PlantBox(dx,1e100,dz);

	for (int i=0; i<n; i++) {

		std::vector<double> row(m); // m columns
		for (int j=0; j<m; j++) {

			Vector3d t(dx/2.+j*dx,0,-dz/2.-i*dz);
			SDF_RotateTranslate g(layer,t);

			if (exact) {
				AnalysisSDF a(*this); // copy everything
				a.crop(&g); // crop exactly
				row.at(j) = a.getSummed(st);
			} else {
				row.at(j) = this->getSummed(st, &g);
			}

		}

		d.at(i)=row; // store the row (n rows)
	}
	delete layer;
	return d;
}

/**
 *  Creates a vertical distribution of the parameter of type @param st (@see RootSystem::ScalarType)
 *
 * @param st        parameter type @see RootSystem::ScalarType
 * @param top       vertical top position (cm)
 * @param bot       vertical bot position (cm)
 * @param n         number of layers (each with a height of (bot-top)/n )
 * \return Vector of size @param n containing the summed parameter in this layser
 */
std::vector<std::vector<AnalysisSDF>> AnalysisSDF::distribution2(double top, double bot, double left, double right, int n, int m) const
{
	std::vector<std::vector<AnalysisSDF>> d(n);
	double dz = (top-bot)/double(n);
	double dx = (right-left)/double(m);
	SDF_PlantBox* layer = new SDF_PlantBox(dx,1.e9,dz);
	// std::cout << "dx " << dx  <<", dz "<< dz << "\n";
	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			Vector3d t(left+dx/2.+j*dx,0.,top-i*dz);
			//std::cout << i<< "," << j << ", " <<t.getString() << "; ";
			SDF_RotateTranslate g(layer,t);
			AnalysisSDF a(*this); // copy everything
			a.crop(&g); // crop exactly
			d.at(i).push_back(a);

		}
	}
	delete layer;
	return d;
}

/**
 * Exports the simulation results with the type from the extension in name
 * (that must be lower case)
 *
 * @param name      file name e.g. output.vtp
 */
void AnalysisSDF::write(std::string name) const
{
	std::ofstream fos;
	fos.open(name.c_str());
	std::string ext = name.substr(name.size()-3,name.size()); // pick the right writer
	if (ext.compare("vtp")==0) {
		std::cout << "writing VTP: " << name << "\n";
		this->writeVTP(fos,{ RootSystem::st_radius, RootSystem::st_type, RootSystem::st_time });
	} else if (ext.compare("txt")==0)  {
		std::cout << "writing text file for Matlab import: "<< name << "\n";
		writeRBSegments(fos);
	} else {
		throw std::invalid_argument("RootSystem::write(): Unkwown file type");
	}
	fos.close();
}

/**
 * Writes a VTP file with @param types data per segment.
 *
 * @param os        typically a file out stream
 * @param types     multiple parameter types (@see RootSystem::ScalarType) that are saved in the VTP file
 */
void AnalysisSDF::writeVTP(std::ostream & os, std::vector<int> types) const
{
	assert(segments.size() == segO.size());
	assert(segments.size() == ctimes.size());
	os << "<?xml version=\"1.0\"?>";
	os << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	os << "<PolyData>\n";
	os << "<Piece NumberOfLines=\""<< segments.size() << "\" NumberOfPoints=\""<< nodes.size()<< "\">\n";
	// data (CellData)
	os << "<CellData Scalars=\" CellData\">\n";
	for (auto i : types) {
		std::vector<double> data = getScalar(i);
		os << "<DataArray type=\"Float32\" Name=\"" << RootSystem::scalarTypeNames.at(i) << "\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
		for (auto const& t : data) {
			os << t << " ";
		}
		os << "\n</DataArray>\n";
	}
	os << "\n</CellData>\n";
	// nodes (Points)
	os << "<Points>\n"<<"<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\" >\n";
	for (auto const& n:nodes) {
		os << n.x << " "<< n.y <<" "<< n.z<< " ";
	}
	os << "\n</DataArray>\n"<< "</Points>\n";
	// segments (Lines)
	os << "<Lines>\n"<<"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (auto const& s:segments) {
		os << s.x << " " << s.y << " ";
	}
	os << "\n</DataArray>\n"<<"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\" >\n";
	for (size_t i=0; i<segments.size(); i++) {
		os << 2*i+2 << " ";
	}
	os << "\n</DataArray>\n";
	os << "\n</Lines>\n";
	//
	os << "</Piece>\n";
	os << "</PolyData>\n" << "</VTKFile>\n";
}

/**
 * Writes the (line)segments of the root system, and
 * mimics the Matlab script getSegments() of RootBox
 *
 * @param os      typically a file out stream
 */
void AnalysisSDF::writeRBSegments(std::ostream & os) const
{
	os << "x1 y1 z1 x2 y2 z2 radius R G B time type \n";
	for (size_t i=0; i<segments.size(); i++) {
		Vector2i s = segments.at(i);
		Vector3d n1 = nodes.at(s.x);
		Vector3d n2 = nodes.at(s.y);
		Root* r = segO.at(i);
		double radius = r->param.a;
		double red = r->getRootTypeParameter()->colorR;
		double green = r->getRootTypeParameter()->colorG;
		double blue = r->getRootTypeParameter()->colorB;
		double time = ctimes.at(i);
		double type = r->param.type;
		os << std::fixed << std::setprecision(4)<< n1.x << " " << n1.y << " " << n1.z << " " << n2.x << " " << n2.y << " " << n2.z << " " <<
				radius << " " << red << " " << green << " " << blue << " " << time<< " " << type << " \n";
	}
}
