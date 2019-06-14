// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Organ.h"

#include "Organism.h"
#include "OrganParameter.h"

namespace CRootBox {

/**
 * Constructor
 */
Organ::Organ(Organism* plant, Organ* parent, int subtype, double delay) : plant(plant), parent(parent), age(-delay)
{
    id = plant->getOrganIndex(); // unique id from the plant
    param_ = plant->getOrganTypeParameter(organType(), subtype)->realize();
}

/**
 * Destructor, tell the kids (bad news)
 */
Organ::~Organ()
{
    for(auto c : children) {
        delete c;
    }
    delete param_; // organ parameters
}

int Organ::organType() const
{
    return Organism::ot_organ;
}

/**
 * Asks the plant for the organ's sub type parameter
 */
OrganTypeParameter* Organ::getOrganTypeParameter() const
{
    return plant->getOrganTypeParameter(this->organType(), param_->subType);
}

/**
 * Calls sub organs (children)
 */
void Organ::simulate(double dt, bool silence)
{
//    if (alive) {
//        age += dt;
//        if (active && (age>0)) {
//            for (auto& c : children)  {
//                c->simulate(dt);
//            }
//        }
//    }
}

/**
 * Adds a node to the organ.
 *
 * For simplicity nodes can not be deleted, organs can only become deactivated or die
 *
 * @param n        new node
 * @param id       global node index
 * @param t        exact creation time of the node
 */
void Organ::addNode(Vector3d n, int id, double t)
{
    nodes.push_back(n); // node
    nodeIds.push_back(id); // new unique id
    nodeCTs.push_back(t); // exact creation time
}

/**
 * Adds the node with the next global index to the root.
 *
 * For simplicity nodes can not be deleted, organs can only become deactivated or die
 *
 * @param n        the new node
 * @param t        exact creation time of the node
 */
void Organ::addNode(Vector3d n, double t)
{
    addNode(n,plant->getNodeIndex(),t);
}

/**
 *
 */
std::vector<Vector2i> Organ::getSegments(int otype) const
{
    if (this->nodes.size()>1) {
        int ot = this->organType();
        if ((otype<0) || (otype==ot)) {
            std::vector<Vector2i> segs = std::vector<Vector2i>(nodes.size()-1);
            for (size_t i=0; i<nodes.size()-1; i++) {
                Vector2i s(getNodeId(i),getNodeId(i+1));
                segs[i] = s;
            }
            return segs;
        }
    }
    return std::vector<Vector2i>(0);
}

/**
 * Returns the organs as sequential list,
 * copies only organs with more than 1 node.
 *
 * \return sequential list of organs
 */
std::vector<Organ*> Organ::getOrgans(int otype)
{
    std::vector<Organ*> v = std::vector<Organ*>();
    this->getOrgans(otype, v);
    return v;
}

/**
 * Returns the organs as sequential list,
 * copies only organs with more than 1 node.
 *
 * @param v     adds the organ sub tree to this vector
 */
void Organ::getOrgans(int otype, std::vector<Organ*>& v)
{
    if (this->nodes.size()>1) {
        int ot = this->organType();
        if ((otype<0) || (otype==ot)) {
            v.push_back(this);
        }
    }
    for (const auto& c : this->children) {
        c->getOrgans(otype,v);
    }
}

/**
 * Returns the parameter called @param name
 */
double Organ::getParameter(std::string name) const {
    double r = std::numeric_limits<double>::quiet_NaN(); // default if name is unknown
    if (name=="one") { r = 1; } // e.g. for counting the organs
    if (name=="id") { r = id; }
    if (name=="organType") { r = this->organType(); }
    if (name=="subType") { r = this->param_->subType; }
    if (name=="alive") { r = alive; }
    if (name=="active") { r = active; }
    if (name=="age") { r = age; }
    if (name=="length") { r = length; }
    if (name=="order") {
        r = 0;
        const Organ* p = this;
        while (p->parent != nullptr) {
            r++;
            p = p->parent; // up the organ tree
        }
    }
    if (name=="nubmer_of_children") { r = children.size(); }
    //    if ((name=="creation_time") && (nodeCTs.size()>0)) { return nodeCTs.at(0); } // check if they make sense
    //    if ((name=="emergence_time") && (nodeCTs.size()>1)) { return nodeCTs.at(1); }
    return r;
}

/**
 * Writes RSML root tag
 *
 * @param cout      typically a file out stream
 * @param indent    we care for looks
 */
void Organ::writeRSML(std::ostream & cout, std::string indent) const
{
    if (this->nodes.size()>1) {
        cout << indent << "<root id=\"" <<  id << "\">\n";  // open root
        /* geometry tag */
        cout << indent << "\t<geometry>\n"; // open geometry
        cout << indent << "\t\t<polyline>\n"; // open polyline
        cout << indent << "\t\t\t" << "<point ";
        Vector3d v = this->nodes.at(0);
        cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
        int n = 1; //this->plant->rsmlReduction;
        for (size_t i = 1; i<nodes.size()-1; i+=n) {
            cout << indent << "\t\t\t" << "<point ";
            Vector3d v = nodes.at(i);
            cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
        }
        cout << indent << "\t\t\t" << "<point ";
        v = nodes.back();
        cout << "x=\"" << v.x << "\" y=\"" << v.y << "\" z=\"" << v.z << "\"/>\n";
        cout << indent << "\t\t</polyline>\n"; // close polyline
        cout << indent << "\t</geometry>\n"; // close geometry
        /* properties */
        cout << indent <<"\t<properties>\n"; // open properties
        // TODO
        cout << indent << "\t</properties>\n"; // close properties
        /* functions */
        cout << indent << "\t<functions>\n"; // open functions
        cout << indent << "\t\t<function name='emergence_time' domain='polyline'>\n"; // open function
        cout << indent << "\t\t\t" << "<sample>" << nodeCTs.at(0) << "</sample>\n";
        for (size_t i = 1; i<nodeCTs.size()-1; i+=n) {
            cout << indent << "\t\t\t" << "<sample>" << nodeCTs.at(i) << "</sample>\n";

        }
        cout << indent << "\t\t\t" << "<sample>" << nodeCTs.back() << "</sample>\n";
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
std::string Organ::toString() const
{
    std::stringstream str;
    str << "Organ #"<< id <<": sub type "<< param_->subType << ", length "<< length << " cm, age " << age
        << " days, alive " << alive << ", active " << active << ", number of nodes" << this->getNumberOfNodes()
        << ", with "<< children.size() << " successors\n";
    return str.str();
}

}
