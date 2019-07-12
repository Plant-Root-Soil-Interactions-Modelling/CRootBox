// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "Organism.h"

#include "OrganParameter.h"
#include "Organ.h"

#include <stdexcept>
#include <iostream>
#include <ctime>
#include <numeric>

namespace CRootBox {

/**
 * Copy constructor
 */
Organism::Organism(const Organism& o): organParam(o.organParam), simtime(o.simtime),
    organId(o.organId), nodeId(o.nodeId), gen(o.gen), UD(o.UD), ND(o.ND)
{
    std::cout << "copy organism \n" << std::flush;
    baseOrgans.resize(o.baseOrgans.size());  // copy base organs
    for (int i=0; i<baseOrgans.size(); i++) {
        baseOrgans[i] = o.baseOrgans[i]->copy(this);
    }
    for (int ot = 0; ot < numberOfOrganTypes; ot++) { // copy organ type parameters
        for (auto& otp : organParam[ot]) {
            otp.second = otp.second->copy(this);
        }
    }
}

/*
 * Destructor: deletes all base organs, and organ type parameters
 */
Organism::~Organism()
{
    for(auto o :baseOrgans) { // delete base organs
        delete o;
    }
    for (int ot = 0; ot < numberOfOrganTypes; ot++) {  // delete organ type parameters
        for (auto& otp : organParam[ot]) {
            delete otp.second;
        }
    }
}

/**
 * Copies the organ type parameters of a specific organ type into a vector
 *
 * @param ot    the organ type
 */
std::vector<OrganTypeParameter*> Organism::getOrganTypeParameter(int ot) const
{
    std::vector<OrganTypeParameter*>  otps = std::vector<OrganTypeParameter*>(0);
    for (auto& otp : organParam[ot]) {
        otps.push_back(otp.second);
    }
    return otps;
}

/**
 * Returns an organ type parameter of a specific organ type and sub type
 *
  * @param ot       the organ type (e.g. ot_root)
  * @param subType  the sub type (e.g. root type)
  * @return         the respective type parameter
 */
OrganTypeParameter* Organism::getOrganTypeParameter(int ot, int subtype) const
{
    try {
        //        std::cout << "reading organ type " << ot << " sub type " << subtype <<": ";
        //        for (auto& p : organParam[ot]) {
        //            std::cout << p.first;
        //        }
        //        std::cout << "\n" << std::flush;
        return organParam[ot].at(subtype);
    } catch(const std::out_of_range& oor) {
        std::cout << "Organism::getOrganTypeParameter: Organ type parameter of sub type " << subtype << " was not set \n" << std::flush;
        throw;
    }
}

/**
 *  Sets the  type parameter, subType and organType defined within p
 *  Deletes the old parameter if there is one, takes ownership of the new one
 *
 *  @param p    the organ type parameter
 */
void Organism::setOrganTypeParameter(OrganTypeParameter* p)
{
    int otype = p->organType;
    int subtype = p->subType;
    try {
        delete organParam[otype][subtype];
    } catch (std::exception& e) {
        // did not exist, nothing to delete
    }
    organParam[otype][subtype] = p;
    // std::cout << "setting organ type " << otype << " sub type " << subtype << "\n";
}

/**
 * Overwrite if there is the need for additional initializations,
 * before simulation starts.
 *
 * e.g. initialization of GrowthFunctions, TropismFunctions, set up base Organs
 */
void Organism::initialize()
{ }

/**
 * Simulates the development of the organism in a time span of @param dt days.
 *
 * @param dt        time step [day]
 * @param verbose   turns console output on or off
 */
void Organism::simulate(double dt, bool verbose)
{
    if (verbose) {
        std::cout << "Organism::simulate: from "<< simtime << " to " << simtime+dt << " days" << std::endl;
    }
    for (const auto& r : baseOrgans) {
        r->simulate(dt, verbose);
    }
    simtime+=dt;
}

/**
 * Creates a sequential list of organs. Considers only organs with more than 1 node.
 *
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 * @return Sequential list of organs. If there is less than one node,
 * or another organ type is expected, an empty vector is returned.
 */
std::vector<Organ*> Organism::getOrgans(int ot) const
{
    std::vector<Organ*> organs = std::vector<Organ*>(0);
    for (const auto& o : this->baseOrgans) {
        o->getOrgans(ot, organs);
    }
    return organs;
}

/**
 * Returns a single scalar parameter for each organ as sequential list,
 * corresponding to the sequential organ list, see Organism::getOrgans.
 *
 * This method is mostly for post processing, since it is flexible but slow.
 *
 * @param name      name of the parameter
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 * @param organs    optional a predefined sequential organ list can be used (@param ot is ignored in this case)
 * @return A vector of one parameter values per each organ, if unknown NaN
 */
std::vector<double> Organism::getParameters(std::string name, int ot, std::vector<Organ*> organs) const
{
    if (organs.empty()) {
        organs = getOrgans(ot);
    }
    std::vector<double> p = std::vector<double>(organs.size());
    for (int i=0; i<organs.size(); i++) {
        p[i] = organs[i]->getParameter(name);
    }
    return p;
}

/**
 * Returns the summed parameter, obtained by Organism::getParameters (e.g. getSummed("length"))
 *
 * @param name      name of the parameter
 * @param ot        the expected organ type, where -1 denotes all organ types (default).
 * @return          the summed up value
 */
double Organism::getSummed(std::string name, int ot) const
{
    auto v = getParameters(name, ot);
    return std::accumulate(v.begin(), v.end(), 0.0);
}

/**
 *
 */
std::vector<std::vector<Vector3d>> Organism::getPolylines(int ot) const
{
    auto organs = getOrgans(ot);
    std::vector<std::vector<Vector3d>> nodes = std::vector<std::vector<Vector3d>>(organs.size()); // reserve big enough vector
    for (size_t j=0; j<organs.size(); j++) {
        std::vector<Vector3d>  rn = std::vector<Vector3d>(organs[j]->getNumberOfNodes());
        for (size_t i=0; i<organs[j]->getNumberOfNodes(); i++) { // loop over all nodes of all organs
            rn.at(i) = organs[j]->getNode(i);
        }
        nodes[j] = rn;
    }
    return nodes;
}

/**
 *
 */
std::vector<std::vector<double>> Organism::getPolylinesIds(int ot) const
{
}

/**
 *
 */
std::vector<std::vector<double>> Organism::getPolylinesCTs(int ot) const
{
}

/**
 *
 */
std::vector<Vector3d> Organism::getNodes() const
{
    auto organs = getOrgans();
    std::vector<Vector3d> nv = std::vector<Vector3d>(getNumberOfNodes()); // reserve big enough vector
    // copy initial nodes (organs might not have developed)
    for (const auto& r : baseOrgans) {
        nv.at(r->getNodeId(0)) = r->getNode(0);
    }
    // copy organ nodes
    for (const auto& r : organs) {
        for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
            nv.at(r->getNodeId(i)) = r->getNode(i); // pray that ids are correct
        }
    }
    return nv;
}

/**
 *
 */
std::vector<double> Organism::getNodeCTs() const
{
    auto organs = getOrgans();
    std::vector<double> nv = std::vector<double>(getNumberOfNodes()); // reserve big enough vector
    // copy initial nodes (organs might not have developed)
    for (const auto& r : baseOrgans) {
        nv.at(r->getNodeId(0)) = r->getNodeCT(0);
    }
    // copy organ nodes
    for (const auto& r : organs) {
        for (size_t i=0; i<r->getNumberOfNodes(); i++) { // loop over all nodes of all roots
            nv.at(r->getNodeId(i)) = r->getNodeCT(i); // pray that ids are correct
        }
    }
    return nv;
}

/**
 *
 */
std::vector<Vector2i> Organism::getSegments(int ot) const
{
    auto organs = getOrgans(ot);
    std::vector<Vector2i> segs = std::vector<Vector2i>(0);
    for (const auto& o : organs) {
        auto s = o->getSegments();
        segs.insert(segs.end(), s.begin(), s.end()); // append s;
    }
    return segs;
}

/**
 *
 */
std::vector<double> Organism::getSegmentCTs(int ot) const
{
    auto nodeCT = getNodeCTs();
    auto segs = getSegments(ot);
    std::vector<double> cts = std::vector<double>(segs.size());
    for (int i=0; i<cts.size(); i++) {
        cts[i] = nodeCT[segs[i].y]; // segment creation time is the node creation time of the second node
    }
    return cts;
}

/**
 *
 */
std::vector<Organ*> Organism::getSegmentOrigins(int ot) const
{
    auto organs = getOrgans(ot);
    std::vector<Organ*> segs = std::vector<Organ*>(0);
    for (const auto& o : organs) {
        auto s = o->getSegments();
        for (int i=0; i<s.size(); i++) {
            segs.push_back(o);
        }
    }
    return segs;
}

/**
 * @return Quick info about the object for debugging
 */
std::string Organism::toString() const
{
    std::stringstream str;
    str << "Organism with "<< baseOrgans.size() <<" base organs, " << getNumberOfNodes()
                << " nodes, and a total of " << getNumberOfOrgans() << " organs, after " << getSimTime() << " days";
    return str.str();
}


/**
 * Creates a rsml file with filename @param name.
 *
 * @param name      name of the rsml file
 */
void Organism::writeRSML(std::string name) const
{
    tinyxml2::XMLDocument xmlDoc;
    tinyxml2:: XMLElement* rsml = xmlDoc.NewElement("rsml"); // RSML
    tinyxml2:: XMLElement* meta = getRSMLMetadata(xmlDoc);
    tinyxml2:: XMLElement* scene = getRSMLScene(xmlDoc);
    rsml->InsertEndChild(meta);
    rsml->InsertEndChild(scene);
    xmlDoc.InsertEndChild(rsml);
    xmlDoc.SaveFile(name.c_str());
}

/**
 *
 */
tinyxml2:: XMLElement* Organism::getRSMLMetadata(tinyxml2::XMLDocument& xmlDoc) const
{
    tinyxml2:: XMLElement* metadata = xmlDoc.NewElement("metadata"); // META
    tinyxml2:: XMLElement* version = xmlDoc.NewElement("version");
    version->SetText(1);
    tinyxml2:: XMLElement* unit = xmlDoc.NewElement("unit");
    unit->SetText("cm");
    tinyxml2:: XMLElement* resolution = xmlDoc.NewElement("resolution");
    resolution->SetText(1);
    tinyxml2:: XMLElement* last = xmlDoc.NewElement("last-modified");
    std::time_t t = std::time(0);
    std::tm* now = std::localtime(&t);
    std::string s = std::to_string(now->tm_mday)+"-"+std::to_string(now->tm_mon+1)+"-"+std::to_string(now->tm_year + 1900);
    last->SetText(s.c_str());
    tinyxml2:: XMLElement* software = xmlDoc.NewElement("software");
    software->SetText("OrganicBox");
    // todo no image tag (?)
    // todo property-definitions
    // todo time sequence (?)
    metadata->InsertEndChild(version);
    metadata->InsertEndChild(unit);
    metadata->InsertEndChild(resolution);
    metadata->InsertEndChild(last);
    metadata->InsertEndChild(software);
    // todo insert remaining tags
    return metadata;
}

/**
 *
 */
tinyxml2:: XMLElement* Organism::getRSMLScene(tinyxml2::XMLDocument& xmlDoc) const
{
    tinyxml2:: XMLElement* scene = xmlDoc.NewElement("scene");
    tinyxml2:: XMLElement* plant = xmlDoc.NewElement("plant");
    scene->InsertEndChild(plant);
    for (auto& o: baseOrgans) {
        o->writeRSML(xmlDoc, plant);
    }
    return scene;
}

/**
 * Sets the seed of the organisms random number generator.
 * In order to obtain two exact same organisms call before Organism::initialize().
 *
 * @param seed      the random number generator seed
 */
void Organism::setSeed(unsigned int seed)
{
    this->gen = std::mt19937(seed);
}


} // namespace
