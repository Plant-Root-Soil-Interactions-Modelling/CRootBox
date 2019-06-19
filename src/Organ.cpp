// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Organ.h"

#include "Organism.h"
#include "OrganParameter.h"

#include <iostream>

namespace CRootBox {

/*
 * Deep copies the organ @param o into a new plant @param plant.
 * All children are deep copied, plant and parent pointers are updated.
 *
 * @param o         the organ to be copied
 * @param plant     the plant the copied organ will be part of
 */
Organ::Organ(const Organ& o, Organism* plant):
        plant(plant), // different plant then o.plant
        parent(nullptr),
        id(o.id),
        param_(o.param_),
        alive(o.alive),
        active(o.active),
        age(o.age),
        length(o.length),
        nodes(o.nodes),
        nodeIds(o.nodeIds),
        nodeCTs(o.nodeCTs)
{
    children = std::vector<Organ*>(o.children.size());
    for (size_t i=0; i< o.children.size(); i++) {
        children[i] = new Organ(*o.children[i], plant); // copy lateral
        children[i]->setParent(this);
    }
}

/**
 * The constructor is used for simulation.
 * The organ parameters are chosen from random distributions within the the OrganTypeParameter class.
 * The next organ id is retrieved from the plant,
 * and the organ starts growing after a delay (starts with age = -delay).
 *
 * @param plant     the plant the new organ will be part of
 * @param parent    the parent organ, equals nullptr if there is no parent
 * @param ot        organ type
 * @param st        sub type of the organ type, e.g. different root types
 * @param delay     time delay in days when the organ will start to grow
 */
Organ::Organ(Organism* plant, Organ* parent, int ot, int st, double delay): plant(plant), parent(parent),
    id(plant->getOrganIndex()),  // unique id from the plant
    param_(plant->getOrganTypeParameter(ot, st)->realize()), // draw specific parameters from random distributions
    age(-delay)
{ }

/**
 * Destructor deletes all children, and its parameter class
 */
Organ::~Organ()
{
    for(auto c : children) {
        delete c;
    }
    delete param_; // organ parameters
}

/**
 * \return The organ type, which is a coarse classification of the organs.
 * Currently there are: ot_organ (for unspecified organs), ot_seed, ot_root, ot_stem, and ot_leaf.
 * There can be different classes with the same organ type.
 */
int Organ::organType() const
{
    return Organism::ot_organ;
}

/**
 * \return The organ type parameter is retrieved from the plant
 */
OrganTypeParameter* Organ::getOrganTypeParameter() const
{
    return plant->getOrganTypeParameter(this->organType(), param_->subType);
}

/**
 * The member function should be overwritten to describe the organ development in a duration of @param dt days.
 *
 * desciption what has to be implemented TODO
 *
 */
void Organ::simulate(double dt, bool verbose)
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
 * By default the organ is represented by a polyline,
 * i.e. the segments of the nodes {n1, n2, n3, n4}, are the indices { [i1,i2], [i2,i3], [i3,i4] }.
 * For other geometries this member function must be overwritten.
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 *
 * \return A vector of line segments, where each line segment is described as two global node indices.
 * If there are less than two nodes, or another organ type is expected, an empty vector is returned.
 */
std::vector<Vector2i> Organ::getSegments(int ot) const
{
    if (this->nodes.size()>1) {
        if ((ot<0) || (ot==this->organType())) {
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
 * Returns the organs as sequential list, copies only organs with more than one node.
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 *
 * \return Sequential list of organs. If there is less than one node,
 * or another organ type is expected, an empty vector is returned.
 */
std::vector<Organ*> Organ::getOrgans(int ot)
{
    std::vector<Organ*> v = std::vector<Organ*>();
    this->getOrgans(ot, v);
    return v;
}

/**
 * Returns the organs as sequential list, copies only organs with more than one node.
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 * @param v         vector of organs where the subtree is added,
 *                  only expected organ types with more than one nodes are added.
 */
void Organ::getOrgans(int ot, std::vector<Organ*>& v)
{
    if (this->nodes.size()>1) {
        if ((ot<0) || (ot==this->organType())) {
            v.push_back(this);
        }
    }
    for (const auto& c : this->children) {
        c->getOrgans(ot,v);
    }
}

/**
 * Returns a single scalar parameter called @param name of the organ.
 * This method is for post processing, since it is flexible but slow.
 * Overwrite to add more parameters for specific organs.
 *
 * \return The parameter value, if unknown NaN
 */
double Organ::getParameter(std::string name) const {
    if (name=="length") { return getLength(); }
    if (name=="age") { return getAge(); }
    if (name=="order") { // count how often it is possible to move up
        int r = 0;
        const Organ* p = this;
        while (p->parent != nullptr) {
            r++;
            p = p->parent; // up the organ tree
        }
        return r;
    }
    if (name=="one") { return 1; } // e.g. for counting the organs
    if (name=="id") { return getId(); }
    if (name=="organType") { return this->organType(); }
    if (name=="subType") { return this->param_->subType; }
    if (name=="alive") { return isAlive(); }
    if (name=="active") { return isActive(); }
    if (name=="nubmerOfChildren") { return children.size(); }

    // if ((name=="creationTime") && (nodeCTs.size()>0)) { return nodeCTs.at(0); }
    // if ((name=="emergence_time") && (nodeCTs.size()>1)) { return nodeCTs.at(1); }

    double r = std::numeric_limits<double>::quiet_NaN(); // default if name is unknown

    // TODO pass to organ type parameters

    return r;
}

/**
 * Writes RSML root tag TODO switch to tiny xml
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
 * \return Quick info about the object for debugging
 */
std::string Organ::toString() const
{
    std::stringstream str;
    str << "Organ #"<< getId() << ": organ type "<< organType() << " sub type "<< param_->subType << ", length " <<
        getLength() << " cm, age " << getAge() << " days, alive " << isAlive() << ", active " << isActive()
        << ", number of nodes" << this->getNumberOfNodes() << ", with "<< children.size() << " children\n";
    return str.str();
}

}
