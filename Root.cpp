#include "Root.h"

/**
 * Constructor
 *
 * Typically called by the RootSystem::RootSystem(), or Root::createNewRoot().
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
Root::Root(RootSystem* rs, int type, Vector3d pheading, double delay,  Root* parent, double pbl, int pni)
{
  //cout << "Root: constructor \n";
  rootsystem=rs; // remember
  param = rs->getRootTypeParameter(type)->realize(); // throw the dice
  // initial rotation
  double beta = 2*M_PI*rs->rand();
  Matrix3d ons = Matrix3d::ons(pheading);
  ons.times(Matrix3d::rotX(beta));
  double theta = param.theta;
//  if (parent!=nullptr) { // scale if not a baseRoot
//      double scale = rs->getRootParameter(type)->saf->getRelativeValue(parent->getNode(pni),this);
//      theta*=scale;
//  }
  ons.times(Matrix3d::rotZ(theta));
  this->iheading = ons.column(0);  // new initial heading
  //
  age = -delay; // the root starts growing when age>0
  alive = 1; // alive per default
  id = rs->getRootIndex(); // root id
  this->parent = parent;
  parent_base_length=pbl;
  parent_ni=pni;
  length = 0;
  dx = getRootTypeParameter()->dx;
  // initial node
  if (parent!=nullptr) { // the first node of the base roots must be created in initialize
      // dont use addNode for the first node of the root,
      // since this node exists already and does not need a new identifier
      nodes.push_back(parent->getNode(pni));
      nodeIds.push_back(parent->getNodeId(pni));
      netimes.push_back(parent->getNodeETime(pni)+delay);
  }
}

/**
 * Destructor, spread the word
 */
Root::~Root()
{
  for(auto l : laterals) {
      delete l;
  }
}

/**
 * Exact length of the root at a certain age
 *
 * @param age          age of the root [day]
 */
double Root::getLength(double age)
{
  assert(age>=0);
  return rootsystem->gf.at(param.type-1)->getLength(age,param.r,param.getK(),this);
}

/**
 * Creation time of a node at a certain length
 *
 * @param length   length of the root [cm]
 */
double Root::getCreationTime(double length)
{
  assert(length>=0);
  double rootage = getAge(length);
  if (rootage<0) {
      std::cout << "negative root age "<<rootage<<" at length "<< length;
      std::cout.flush();
      throw std::invalid_argument( "bugbugbug" );
  }
  if (parent!=nullptr) {
      double pl = parent_base_length+parent->param.la; // parent length, when this root was created
      double page=parent->getCreationTime(pl);
      assert(page>=0);
      return rootage+page;
  } else {
      return rootage;
  }
}

/**
 * Exact age of the root at a ceratain length
 *
 * @param length   length of the root [cm]
 */
double Root::getAge(double length)
{
  assert(length>=0);
  return rootsystem->gf.at(param.type-1)->getAge(length,param.r,param.getK(),this);
}

RootTypeParameter* Root::getRootTypeParameter()
{
  return rootsystem->getRootTypeParameter(param.type);
}

/**
 * Simulates growth of this root for a time span dt
 *
 * @param dt       time step [day]
 */
void Root::simulate(double dt)
{
  // cout << "Root.simulate(dt) " <<rp.name<< "\n";

  const RootParameter &p = param; // rename

  // increase age
  if (age+dt>p.rlt) { // root life time
      dt=p.rlt-age; // dt might be 0
      alive = false; // this root is dead TODO
  }
  age+=dt;

  if (alive) { // dead roots wont grow

      if (age>0) {

          // children first (always care for the children, even if inactive)
          for (auto l:laterals) {
              l->simulate(dt);
          }

          if (active) {
              // length increment
              double targetlength = getLength(age);
              //cout << "targetlength: "<< targetlength << " at age " << age <<", r="<<p.r<< ", k="<<p.getK()<<"\n";
              double dl = std::max(targetlength-length, double(0)); // length increment

              // create geometry
              if (p.nob>0) { // root has laterals
                  // basal zone
                  if ((dl>0)&&(length<p.lb)) { // length is the current length of the root
                      if (length+dl<=p.lb) {
                          createSegments(dl);
                          length+=dl;
                          dl=0;
                      } else {
                          double ddx = p.lb-length;
                          createSegments(ddx);
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
                              if (i==laterals.size()) { // new lateral
                                  createLateral();
                              }
                              if (length+dl<=s) { // finish within inter-lateral distance i
                                  createSegments(dl);
                                  length+=dl;
                                  dl=0;
                              } else { // grow over inter-lateral distance i
                                  double ddx = s-length;
                                  createSegments(ddx);
                                  dl-=ddx;
                                  length=s;
                              }
                          }
                      }
                      if (dl>0) {
                          if (p.ln.size()==laterals.size()) { // new lateral (the last one)
                              createLateral();
                          }
                      }
                  }
                  // apical zone
                  if (dl>0) {
                      createSegments(dl);
                      length+=dl;
                  }
              } else { // no laterals
                  if (dl>0) {
                      createSegments(dl);
                      length+=dl;
                  }
              } // if laterals
          } // if active
          active = getLength(std::max(age,0.))<(p.getK()-dx/10); // become inactive, if final length is nearly reached
      }
  } // if alive
}

/**
 * Creates a new lateral by calling RootSystem::createNewRoot().
 *
 * Overwrite this method to implement more sezialized root classes.
 */
void Root::createLateral()
{
  const RootParameter &p = param; // rename

  int lt = rootsystem->getRootTypeParameter(p.type)->getLateralType(nodes.back());
  if (lt>0) {

      Vector3d h; // old heading
      if (nodes.size()>1) {
          h = nodes.back().minus(nodes.at(nodes.size()-2)); // getHeading(b-a)
      } else {
          h= iheading;
      }

      double ageLN = this->getAge(length); // age of root when lateral node is created
      double ageLG = this->getAge(length+p.la); // age of the root, when the lateral starts growing (i.e when the apical zone is developed)
      double delay = ageLG-ageLN; // time the lateral has to wait

      Root* lateral = rootsystem->createRoot(lt,  h, delay,  this, length, nodes.size()-1);
      //cout << "Root.createLateral() type "<< lt << "\n";

      laterals.push_back(lateral);
      lateral->simulate(age-ageLN); // pass time overhead (age we want to achieve minus current age)
      //cout << "time overhead " << age-ageLN << "\n";
  }
}

/**
 *  Creates nodes and node emergence times for length l,
 *  and updates the root heading
 *
 *  Cecks that each new segments length is <= dx but >= ddx
 *
 *  @param l       length the root growth [cm]
 */
void Root::createSegments(double l)
{
  assert(l>0);
  double scale = rootsystem->getRootTypeParameter(param.type)->sef->getRelativeValue(nodes.back(),this); // hope this optimized out if not set
  l = l*scale;

  double sl=0; // summed length of created segment
  int n = floor(l/dx);

  for (int i=0; i<n+1; i++) {

      Vector3d h; // current heading
      if (nodes.size()>1) {
          h = nodes.back().minus(nodes.at(nodes.size()-2));
          h.normalize();
      } else {
          h = iheading;
      }

      double sdx; // segment length (<=dx)
      if (i<n) {  // normal case
          sdx = dx;
      } else { // last segment
          sdx = l-n*dx;
          if (sdx<ddx) { // ensure minimal length
              // std::cout << "Root::createSegments() dropped very small segment\n";
              sdx=0;
              break;
          }
      }
      sl+=sdx;

      //std::cout << "heading: " << h.toString() << "  ";
      Matrix3d ons = Matrix3d::ons(h);
      //cout << "\n" << ons.getString() << "\n";
      Vector2d ab = rootsystem->tf.at(param.type-1)->getHeading(nodes.back(),ons,sdx,this);
      ons.times(Matrix3d::rotX(ab.y));
      ons.times(Matrix3d::rotZ(ab.x));
      //std::cout << "norm heading " << ons.column(0).toString() << "\n"; // == 1
      Vector3d newdx = Vector3d(ons.column(0).times(sdx));
      //std::cout << "   node " << newdx.toString() <<"\n";

      Vector3d newnode = Vector3d(nodes.back().plus(newdx));
      double et = this->getCreationTime(length+sl);
      if (std::isnan(et)) {
          std::cout << "Creation time is nan \n";
          std::cout.flush();
      }
      addNode(newnode,et);
  } // for
}

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
  for (auto const& l:this->laterals) {
      l->getRoots(v);
  }
}

/**
 * Adds the next node to the root.
 *
 * Add nodes only with this function! For simplicity nodes can not be deleted, and roots can only become deactivated by dying
 *
 * @param n        the new node
 * @param t        exact creation time of the node
 */
void Root::addNode(Vector3d n, double t)
{
  assert(t>=0.);
  nodes.push_back(n); // node
  nodeIds.push_back(rootsystem->getNodeIndex()); // new unique id
  netimes.push_back(t); // exact creation time
}

/**
 * writes RSML root tag
 *
 * @param cout      typically a file out stream
 * @param indent    we care for looks
 */
void Root::writeRSML(std::ostream & cout, std::string indent) const
{
  cout << indent << "<root id=\"" <<  id << "\">\n";  // open root
  /* geometry tag */
  cout << indent << "\t<geometry>\n"; // open geometry
  cout << indent << "\t\t<polyline>\n"; // open polyline
  for (size_t i = 0; i<nodes.size(); i++) {
      cout << indent << "\t\t\t" << "<point ";
      Vector3d v = nodes.at(i);
      cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
  }
  cout << indent << "\t\t</polyline>\n"; // close polyline
  cout << indent << "\t</geometry>\n"; // close geometry
  /* properties */
  cout << indent <<"\t<properties>\n"; // open propertie
  // TODO
  cout << indent << "\t</properties>\n"; // close properties
  /* laterals roots */
  for (size_t i = 0; i<laterals.size(); i++) {
      laterals[i]->writeRSML(cout,indent.append("\t"));
  }
  cout << indent << "\t<functions>\n"; // open functions
  // TODO
  cout << indent << "\t</functions>\n"; // close functions
  cout << indent << "</root>\n"; // close root
}
