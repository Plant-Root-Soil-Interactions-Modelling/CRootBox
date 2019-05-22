// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PY_ROOTBOX_H_
#define PY_ROOTBOX_H_

/**
 *  A Python module for CRootbox based on boost.python
 *
 *  build a shared library from this file
 *  put comment to line 16 to ignore this file
 */

#define PYTHON_WRAPPER // UNCOMMENT TO BUILD SHARED LIBRARY

#ifdef PYTHON_WRAPPER

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/call.hpp>

#include "mymath.h"
#include "sdf.h"
#include "RootSystem.h"
#include "sdf_rs.h"
#include "analysis.h"
#include "../examples/Exudation/example_exudation.h"

namespace CRootBox {

using namespace boost::python;

/*
 * Functions overloading (by hand, there are also macros available)
 *
 * Tutorial example:
 * bool    (X::*fx1)(int)              = &X::f;
 * bool    (X::*fx2)(int, double)      = &X::f;
 * bool    (X::*fx3)(int, double, char)= &X::f;
 * int     (X::*fx4)(int, int, int)    = &X::f;
 *
 */
Vector3d (Vector3d::*times1)(const double) const = &Vector3d::times;
double (Vector3d::*times2)(const Vector3d&) const = &Vector3d::times;
void (Matrix3d::*times3)(const Matrix3d&) = &Matrix3d::times;
Vector3d (Matrix3d::*times4)(const Vector3d&) const = &Matrix3d::times;

std::string (SignedDistanceFunction::*writePVPScript)() const = &SignedDistanceFunction::writePVPScript; // because of default value

//void (RootSystem::*simulate1)(double dt, bool silence) = &RootSystem::simulate;
//void (RootSystem::*simulate2)() = &RootSystem::simulate;
//void (RootSystem::*simulate3)(double dt, double maxinc, ProportionalElongation* se, bool silence) = &RootSystem::simulate;
//void (RootSystem::*initialize1)() = &RootSystem::initialize;
//void (RootSystem::*initialize2)(int basal, int shootborne) = &RootSystem::initialize;

void (SegmentAnalyser::*addSegments1)(const PlantBase& plant) = &SegmentAnalyser::addSegments;
void (SegmentAnalyser::*addSegments2)(const SegmentAnalyser& a) = &SegmentAnalyser::addSegments;
void (SegmentAnalyser::*filter1)(std::string name, double min, double max) = &SegmentAnalyser::filter;
void (SegmentAnalyser::*filter2)(std::string name, double value) = &SegmentAnalyser::filter;
double (SegmentAnalyser::*getSummed1)(std::string name) const = &SegmentAnalyser::getSummed;
double (SegmentAnalyser::*getSummed2)(std::string name, SignedDistanceFunction* geometry) const = &SegmentAnalyser::getSummed;
std::vector<double> (SegmentAnalyser::*distribution_1)(std::string name, double top, double bot, int n, bool exact) const = &SegmentAnalyser::distribution;
std::vector<SegmentAnalyser> (SegmentAnalyser::*distribution_2)(double top, double bot, int n) const = &SegmentAnalyser::distribution;
std::vector<std::vector<double>> (SegmentAnalyser::*distribution2_1)(std::string name, double top, double bot, double left, double right, int n, int m, bool exact) const = &SegmentAnalyser::distribution2;
std::vector<std::vector<SegmentAnalyser>> (SegmentAnalyser::*distribution2_2)(double top, double bot, double left, double right, int n, int m) const = &SegmentAnalyser::distribution2;
SegmentAnalyser (SegmentAnalyser::*cut1)(const SDF_HalfPlane& plane) const = &SegmentAnalyser::cut;



/**
 * Default arguments: no idea how to do it by hand, magic everywhere...
 */
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(openFile_overloads,openFile,1,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(simulate1_overloads,simulate,1,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(simulate3_overloads,simulate,3,4);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getValue_overloads,getValue,1,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(tropismObjective_overloads,tropismObjective,5,6);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getNumberOfRoots_overloads,getNumberOfRoots,0,1);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(getNETimes_overloads, getNETimes, 0, 1);


/**
 * Virtual functions
 */
class SoilLookUp_Wrap : public SoilLookUp, public wrapper<SoilLookUp> {
public:

    virtual double getValue(const Vector3d& pos, const Root* root = nullptr) const override {
        return this->get_override("getValue")(pos, root);
    }

    virtual std::string toString() const override {
        return this->get_override("toString")();
    }

};

class Tropism_Wrap : public Tropism, public wrapper<Tropism> {
public:

    //	Tropism_Wrap(): Tropism() { }
    //	Tropism_Wrap(double n,double sigma): Tropism(n,sigma) { } // todo cant get it working with constructors other than ()

    virtual double tropismObjective(const Vector3d& pos, Matrix3d old, double a, double b, double dx, const Root* root = nullptr) override {
        return this->get_override("tropismObjective")(pos, old, a, b, dx, root);
    }

    virtual Tropism* copy() override {
        return this->get_override("copy")();
    }

};



/**
 * Expose classes to Python module
 */
BOOST_PYTHON_MODULE(py_rootbox)
{
    /*
     * general
     */
    class_<std::vector<double>>("std_vector_double_")
        .def(vector_indexing_suite<std::vector<double>>() )
        ;

    class_<std::vector<int>>("std_vector_int_")
        .def(vector_indexing_suite<std::vector<int>>() )
        ;
    /*
     * mymath.h
     */
    class_<Vector2i>("Vector2i", init<>())
        .def(init<int,int>())
        .def(init<Vector2i&>())
        .def_readwrite("x",&Vector2i::x)
        .def_readwrite("y",&Vector2i::y)
        .def("__str__",&Vector2i::toString)
    ;
    class_<std::vector<Vector2i>>("std_vector_Vector2i_")
        .def(vector_indexing_suite<std::vector<Vector2i>>() )
    ;
    class_<Vector3d>("Vector3d", init<>())
        .def(init<double,double,double>())
        .def(init<Vector3d&>())
        .def_readwrite("x",&Vector3d::x)
        .def_readwrite("y",&Vector3d::y)
        .def_readwrite("z",&Vector3d::z)
        .def("normalize",&Vector3d::normalize)
        .def("times",times1)
        .def("times",times2)
        .def("length",&Vector3d::length)
        .def("plus",&Vector3d::plus)
        .def("minus",&Vector3d::minus)
        .def("cross",&Vector3d::cross)
        .def("__str__",&Vector3d::toString)
        .def("__rep__",&Vector3d::toString)
    ;
    class_<std::vector<Vector3d>>("std_vector_Vector3d_")
        .def(vector_indexing_suite<std::vector<Vector3d>>() )
    ;
    class_<Matrix3d>("Matrix3d", init<>())
        .def(init<double,double,double,double,double,double,double,double,double>())
        .def(init<Vector3d&, Vector3d&, Vector3d&>())
        .def(init<Matrix3d&>())
        .def_readwrite("r0",&Matrix3d::r0)
        .def_readwrite("r1",&Matrix3d::r1)
        .def_readwrite("r2",&Matrix3d::r2)
        .def("rotX",&Matrix3d::rotX)
        .def("rotY",&Matrix3d::rotY)
        .def("rotZ",&Matrix3d::rotZ)
        .def("ons",&Matrix3d::ons)
        .def("det",&Matrix3d::det)
        .def("inverse",&Matrix3d::inverse)
        .def("column",&Matrix3d::column)
        .def("row",&Matrix3d::row)
        .def("times",times3)
        .def("times",times4)
        .def("__str__",&Matrix3d::toString)
        .def("__rep__",&Matrix3d::toString)
    ;
    /*
     * sdf.h
     */
    class_<SignedDistanceFunction, SignedDistanceFunction*>("SignedDistanceFunction")
                .def("getDist",&SignedDistanceFunction::getDist)
                .def("writePVPScript", writePVPScript)
                .def("__str__",&SignedDistanceFunction::toString)
                ;
    class_<std::vector<SignedDistanceFunction*>>("std_vector_SDF_")
            .def(vector_indexing_suite<std::vector<SignedDistanceFunction*>>() )
            ;
    class_<SDF_PlantBox, bases<SignedDistanceFunction>>("SDF_PlantBox",init<double,double,double>())
                .def("getDist",&SDF_PlantBox::getDist)
                .def("__str__",&SDF_PlantBox::toString)
                ;
    class_<SDF_PlantContainer, bases<SignedDistanceFunction>>("SDF_PlantContainer",init<>())
                .def(init<double,double,double,double>())
                .def("getDist",&SDF_PlantContainer::getDist)
                .def("__str__",&SDF_PlantContainer::toString)
                ;
    class_<SDF_RotateTranslate, bases<SignedDistanceFunction>>("SDF_RotateTranslate",init<SignedDistanceFunction*,double,int,Vector3d&>())
                .def(init<SignedDistanceFunction*,Vector3d&>())
                .def("getDist",&SDF_RotateTranslate::getDist)
                .def("__str__",&SDF_RotateTranslate::toString)
                ;
    enum_<SDF_RotateTranslate::SDF_Axes>("SDF_Axis")
            .value("xaxis", SDF_RotateTranslate::SDF_Axes::xaxis)
            .value("yaxis", SDF_RotateTranslate::SDF_Axes::yaxis)
            .value("zaxis", SDF_RotateTranslate::SDF_Axes::zaxis)
            ;
    class_<SDF_Intersection, bases<SignedDistanceFunction>>("SDF_Intersection",init<std::vector<SignedDistanceFunction*>>())
                .def(init<SignedDistanceFunction*,SignedDistanceFunction*>())
                .def("getDist",&SDF_Intersection::getDist)
                .def("__str__",&SDF_Intersection::toString)
                ;
    class_<SDF_Union, bases<SDF_Intersection>>("SDF_Union",init<std::vector<SignedDistanceFunction*>>())
                .def(init<SignedDistanceFunction*,SignedDistanceFunction*>())
                .def("getDist",&SDF_Union::getDist)
                .def("__str__",&SDF_Union::toString)
                ;
    class_<SDF_Difference, bases<SDF_Intersection>>("SDF_Difference",init<std::vector<SignedDistanceFunction*>>())
        .def(init<SignedDistanceFunction*,SignedDistanceFunction*>())
        .def("getDist",&SDF_Difference::getDist)
        .def("__str__",&SDF_Difference::toString)
    ;
    class_<SDF_Complement, bases<SignedDistanceFunction>>("SDF_Complement",init<SignedDistanceFunction*>())
        .def("getDist",&SDF_Complement::getDist)
        .def("__str__",&SDF_Complement::toString)
    ;
    class_<SDF_HalfPlane, bases<SignedDistanceFunction>>("SDF_HalfPlane",init<Vector3d&,Vector3d&>())
        .def(init<Vector3d&,Vector3d&,Vector3d&>())
        .def("getDist",&SDF_HalfPlane::getDist)
        .def_readwrite("o", &SDF_HalfPlane::o)
        .def_readwrite("n", &SDF_HalfPlane::n)
        .def_readwrite("p1", &SDF_HalfPlane::p1)
        .def_readwrite("p2", &SDF_HalfPlane::p2)
        .def("__str__",&SDF_HalfPlane::toString)
    ;

    /*
     * OrganParameter.h
     */
    class_<OrganParameter, OrganParameter*>("OrganParameter")
                .def_readwrite("name",&OrganParameter::subType)
        ;
    class_<OrganTypeParameter>("OrganTypeParameter", init<PlantBase*>())
            .def("realize",&OrganTypeParameter::realize, return_value_policy<reference_existing_object>())
            .def("toString",&OrganTypeParameter::toString)
            .def_readwrite("name",&OrganTypeParameter::name)
            .def_readwrite("organType",&OrganTypeParameter::organType)
            .def_readwrite("subType",&OrganTypeParameter::subType)
    ;

//    /*
//     * sdf_rs.h
//     */
//    class_<SDF_RootSystem,SDF_RootSystem*, bases<SignedDistanceFunction>>("SDF_RootSystem", init<std::vector<Vector3d>, std::vector<Vector2i>, std::vector<double>, double>())
//        .def(init<Root&, double>())
//        .def(init<RootSystem&, double>())
//        .def("getDist",&SDF_RootSystem::getDist)
//        .def("__str__",&SDF_RootSystem::toString)
//    ;
//
//    /*
//     * soil.h
//     */
//    class_<SoilLookUp_Wrap, SoilLookUp_Wrap*, boost::noncopyable>("SoilLookUp",init<>())
//                .def("getValue",&SoilLookUp_Wrap::getValue)
//                .def("__str__",&SoilLookUp_Wrap::toString)
//                ;
//    class_<SoilLookUpSDF, SoilLookUpSDF*, bases<SoilLookUp>>("SoilLookUpSDF",init<>())
//                .def(init<SignedDistanceFunction*, double, double, double>())
//                .def_readwrite("sdf", &SoilLookUpSDF::sdf)
//                .def_readwrite("fmax", &SoilLookUpSDF::fmax)
//                .def_readwrite("fmin", &SoilLookUpSDF::fmin)
//                .def_readwrite("slope", &SoilLookUpSDF::slope)
//                .def("__str__",&SoilLookUpSDF::toString)
//                ;
//    class_<MultiplySoilLookUps, MultiplySoilLookUps*, bases<SoilLookUp>>("MultiplySoilLookUps",init<SoilLookUp*, SoilLookUp*>())
//                .def(init<std::vector<SoilLookUp*>>())
//                .def("getValue", &MultiplySoilLookUps::getValue, getValue_overloads())
//                .def("__str__",&MultiplySoilLookUps::toString)
//                ;
//    class_<ProportionalElongation, ProportionalElongation*, bases<SoilLookUp>>("ProportionalElongation",init<>())
//                .def("getValue", &ProportionalElongation::getValue, getValue_overloads())
//                .def("setScale", &ProportionalElongation::setScale)
//                .def("setBaseLookUp", &ProportionalElongation::setBaseLookUp)
//                .def("__str__",&ProportionalElongation::toString)
//                ;
//    class_<Grid1D, Grid1D*, bases<SoilLookUp>>("Grid1D",init<>())
//                .def(init<size_t, std::vector<double>, std::vector<double>>())
//                .def("getValue", &Grid1D::getValue, getValue_overloads())
//                .def("__str__",&Grid1D::toString)
//                .def("map",&Grid1D::map)
//                .def_readwrite("n", &Grid1D::n)
//                .def_readwrite("grid", &Grid1D::grid)
//                .def_readwrite("data", &Grid1D::data)
//                ;
//    class_<EquidistantGrid1D, EquidistantGrid1D*, bases<Grid1D>>("EquidistantGrid1D",init<double, double, size_t>())
//                .def(init<double, double, std::vector<double>>())
//                .def("getValue", &EquidistantGrid1D::getValue, getValue_overloads())
//                .def("__str__",&EquidistantGrid1D::toString)
//                .def("map",&EquidistantGrid1D::map)
//                .def_readwrite("n", &EquidistantGrid1D::n)
//                .def_readwrite("grid", &EquidistantGrid1D::grid)
//                .def_readwrite("data", &EquidistantGrid1D::data)
//                ;
//    /**
//     * tropism.h
//     */
//    class_<Tropism_Wrap, Tropism_Wrap*, boost::noncopyable>("Tropism",init<>())
//                .def("getHeading",&Tropism_Wrap::getHeading)
//                .def("tropismObjective",&Tropism_Wrap::tropismObjective, tropismObjective_overloads())
//                .def("copy",&Tropism_Wrap::copy, return_value_policy<reference_existing_object>())
//                .def("setTropismParameter",&Tropism_Wrap::setTropismParameter)
//                .def("setSeed",&Tropism_Wrap::setSeed)
//                .def("setGeometry",&Tropism_Wrap::setGeometry)
//                .def("rand",&Tropism_Wrap::rand)
//                .def("randn",&Tropism_Wrap::randn)
//                ;
//    class_<Tropism, Tropism*>("TropismBase",init<>()) // Base class for the following tropisms
//                .def("getHeading",&Tropism::getHeading)
//                .def("tropismObjective",&Tropism::tropismObjective, tropismObjective_overloads())
//                .def("copy",&Tropism::copy, return_value_policy<reference_existing_object>())
//                .def("setTropismParameter",&Tropism::setTropismParameter)
//                .def("setSeed",&Tropism::setSeed)
//                .def("setGeometry",&Tropism::setGeometry)
//                .def("rand",&Tropism::rand)
//                .def("randn",&Tropism::randn)
//                ;
//    class_<Gravitropism, Gravitropism*, bases<Tropism>>("Gravitropism",init<double, double>())
//        ;
//    class_<Plagiotropism, Plagiotropism*, bases<Tropism>>("Plagiotropism",init<double, double>())
//        ;
//    class_<Exotropism, Exotropism*, bases<Tropism>>("Exotropism",init<double, double>())
//        ;
//    class_<Hydrotropism, Hydrotropism*, bases<Tropism>>("Hydrotropism",init<double, double, SoilLookUp*>())
//        ;
//    //	class_<CombinedTropism, CombinedTropism*, bases<Tropism>>("CombinedTropism",init<>()) // Todo needs some extra work
//    //	;
//    /*
//     * ModelParameter.h
//     */
//    class_<RootTypeParameter, RootTypeParameter*, bases<OrganTypeParameter>>("RootTypeParameter", init<PlantBase*>())
//                .def("realize",&RootTypeParameter::realize, return_value_policy<reference_existing_object>())
//                .def("getLateralType",&RootTypeParameter::getLateralType)
//                .def("getK",&RootTypeParameter::getK)
//                .def_readwrite("type", &RootTypeParameter::type)
//                .def_readwrite("lb", &RootTypeParameter::lb)
//                .def_readwrite("lbs", &RootTypeParameter::lbs)
//                .def_readwrite("la", &RootTypeParameter::la)
//                .def_readwrite("las", &RootTypeParameter::las)
//                .def_readwrite("ln", &RootTypeParameter::ln)
//                .def_readwrite("lns", &RootTypeParameter::lns)
//                .def_readwrite("nob", &RootTypeParameter::nob)
//                .def_readwrite("nobs", &RootTypeParameter::nobs)
//                .def_readwrite("r", &RootTypeParameter::r)
//                .def_readwrite("rs", &RootTypeParameter::rs)
//                .def_readwrite("a", &RootTypeParameter::a)
//                .def_readwrite("a_s", &RootTypeParameter::as) // as is a keyword in python
//                .def_readwrite("colorR", &RootTypeParameter::colorR)
//                .def_readwrite("colorG", &RootTypeParameter::colorG)
//                .def_readwrite("colorB", &RootTypeParameter::colorB)
//                .def_readwrite("tropismT", &RootTypeParameter::tropismT)
//                .def_readwrite("tropismN", &RootTypeParameter::tropismN)
//                .def_readwrite("tropismS", &RootTypeParameter::tropismS)
//                .def_readwrite("dx", &RootTypeParameter::dx)
//                .def_readwrite("theta", &RootTypeParameter::theta)
//                .def_readwrite("thetas", &RootTypeParameter::thetas)
//                .def_readwrite("rlt", &RootTypeParameter::rlt)
//                .def_readwrite("rlts", &RootTypeParameter::rlts)
//                .def_readwrite("gf", &RootTypeParameter::gf)
//                .def_readwrite("name", &RootTypeParameter::name)
//                .def_readwrite("successor", &RootTypeParameter::successor)
//                .def_readwrite("successorP", &RootTypeParameter::successorP)
//                .def_readwrite("se", &RootTypeParameter::se)
//                .def_readwrite("sa", &RootTypeParameter::sa)
//                .def_readwrite("sbp", &RootTypeParameter::sbp)
//                .def("__str__",&RootTypeParameter::toString)
//                ;
//    class_<RootParameter, RootParameter*, bases<OrganParameter> >("RootParameter", init<>())
//                .def(init<int , double, double, const std::vector<double>&, int, double, double, double, double>())
//                .def("set",&RootParameter::set)
//                .def_readwrite("type", &RootParameter::type)
//                .def_readwrite("lb", &RootParameter::lb)
//                .def_readwrite("la", &RootParameter::la)
//                .def_readwrite("ln", &RootParameter::ln)
//                .def_readwrite("nob", &RootParameter::nob)
//                .def_readwrite("r", &RootParameter::r)
//                .def_readwrite("a", &RootParameter::a)
//                .def_readwrite("theta", &RootParameter::theta)
//                .def_readwrite("rlt", &RootParameter::rlt)
//                .def("getK",&RootParameter::toString)
//                .def("__str__",&RootParameter::toString)
//                ;
//    class_<RootSystemParameter>("RootSystemParameter", init<>())
//                .def("set",&RootSystemParameter::set)
//                .def_readwrite("seedPos", &RootSystemParameter::seedPos)
//                .def_readwrite("firstB", &RootSystemParameter::firstB)
//                .def_readwrite("delayB", &RootSystemParameter::delayB)
//                .def_readwrite("maxB", &RootSystemParameter::maxB)
//                .def_readwrite("nC", &RootSystemParameter::nC)
//                .def_readwrite("firstSB", &RootSystemParameter::firstSB)
//                .def_readwrite("delaySB", &RootSystemParameter::delaySB)
//                .def_readwrite("delayRC", &RootSystemParameter::delayRC)
//                .def_readwrite("nz", &RootSystemParameter::nz)
//                .def("__str__",&RootSystemParameter::toString)
//                ;
//    /**
//     * Root.h (no members, just data)
//     */
//    class_<Root, Root*, bases<Organ>>("Root", init<RootSystem*, int, Vector3d, double, Root*, double, int>())
//            .def(init<Root&>())
//            .def("__str__",&Root::toString)
//            .def_readwrite("rootsystem", &Root::rootsystem)
//            .def_readwrite("param", &Root::param)
//            .def_readwrite("id", &Root::id)
//            .def_readwrite("parent_base_length", &Root::parent_base_length)
//            .def_readwrite("parent_ni", &Root::parent_ni)
//            .def_readwrite("alive", &Root::alive)
//            .def_readwrite("active", &Root::active)
//            .def_readwrite("age", &Root::age)
//            .def_readwrite("length", &Root::length)
//            .def_readwrite("parent", &Root::parent)
//            .def_readwrite("laterals", &Root::laterals)
//            ;
//    class_<std::vector<Root*>>("std_vector_Root_")
//            .def(vector_indexing_suite<std::vector<Root*>>() )
//            ;
//    class_<std::vector<std::vector<Vector3d>>>("std_vector_vector_Vector3d_")
//            .def(vector_indexing_suite<std::vector<std::vector<Vector3d>>>() )
//            ;
//    class_<std::vector<std::vector<double>>>("std_vector_vector_double_")
//            .def(vector_indexing_suite<std::vector<std::vector<double>>>() )
//            ;
////    /*
//     * PlantBase.h
//     */
//    class_<PlantBase, PlantBase*>("PlantBase")
//        ;
//    /*
//     * RootSystem.h
//     */
//    class_<RootSystem, RootSystem*, bases<PlantBase>>("RootSystem", init<>()) //
//             .def(init<RootSystem&>())
//             .def("setRootTypeParameter", &RootSystem::setRootTypeParameter)
//             .def("getRootTypeParameter", &RootSystem::getRootTypeParameter, return_value_policy<reference_existing_object>())
//             .def("setRootSystemParameter", &RootSystem::setRootSystemParameter)
//             .def("getRootSystemParameter", &RootSystem::getRootSystemParameter, return_value_policy<reference_existing_object>()) // tutorial: "naive (dangerous) approach"
//             .def("openFile", &RootSystem::openFile, openFile_overloads())
//             .def("setGeometry", &RootSystem::setGeometry)
//             .def("setSoil", &RootSystem::setSoil)
//             .def("reset", &RootSystem::reset)
//             .def("initialize", initialize1)
//             .def("initialize", initialize2)
//             .def("setTropism", &RootSystem::setTropism)
//             .def("simulate",simulate1, simulate1_overloads())
//             .def("simulate",simulate2)
//             .def("simulate",simulate3, simulate3_overloads())
//             .def("getSimTime", &RootSystem::getSimTime)
//             .def("getNumberOfNodes", &RootSystem::getNumberOfNodes)
//             .def("getNumberOfSegments", &RootSystem::getNumberOfSegments)
//             .def("getNumberOfRoots", &RootSystem::getNumberOfRoots, getNumberOfRoots_overloads())
//             .def("getRoots", &RootSystem::getRoots)
//             .def("getBaseRoots", &RootSystem::getBaseRoots)
//             .def("getNodes", &RootSystem::getNodes)
//             .def("getPolylines", &RootSystem::getPolylines)
//             .def("getSegments", &RootSystem::getSegments)
//             .def("getShootSegments", &RootSystem::getShootSegments)
//             .def("getSegmentsOrigin", &RootSystem::getSegmentsOrigin)
//             .def("getNETimes",&RootSystem::getNETimes, getNETimes_overloads())
//             .def("getScalar", &RootSystem::getScalar)
//             .def("getRootTips", &RootSystem::getRootTips)
//             .def("getRootBases", &RootSystem::getRootBases)
//             .def("write", &RootSystem::write)
//             .def("setSeed",&RootSystem::setSeed)
//             .def("getNumberOfNewNodes",&RootSystem::getNumberOfNewNodes)
//             .def("getNumberOfNewRoots",&RootSystem::getNumberOfNewRoots)
//             .def("getUpdatedNodeIndices",&RootSystem::getUpdatedNodeIndices)
//             .def("getUpdatedNodes",&RootSystem::getUpdatedNodes)
//             .def("getNewNodes",&RootSystem::getNewNodes)
//             .def("getNewSegments",&RootSystem::getNewSegments)
//             .def("getNewSegmentsOrigin",&RootSystem::getNewSegmentsOrigin)
//             .def("getNewSegmentsTimes",&RootSystem::getNewSegmentsTimes)
//             .def("push",&RootSystem::push)
//             .def("pop",&RootSystem::pop)
//             .def("rand",&RootSystem::rand)
//             .def("randn",&RootSystem::randn)
//             .def("setPoreGeometry",&RootSystem::setPoreGeometry)
//             .def("setPoreConductivity",&RootSystem::setPoreConductivity)
//             .def("setPoreLocalAxes",&RootSystem::setPoreLocalAxes)
//             .def("applyPoreConductivities",&RootSystem::applyPoreConductivities)
//             ;
//    enum_<RootSystem::TropismTypes>("TropismType")
//            .value("plagio", RootSystem::TropismTypes::tt_plagio)
//            .value("gravi", RootSystem::TropismTypes::tt_gravi)
//            .value("exo", RootSystem::TropismTypes::tt_exo)
//            .value("hydro", RootSystem::TropismTypes::tt_hydro)
//            ;
//    enum_<RootSystem::GrowthFunctionTypes>("GrowthFunctionType")
//            .value("negexp", RootSystem::GrowthFunctionTypes::gft_negexp)
//            .value("linear", RootSystem::GrowthFunctionTypes::gft_linear)
//            ;
//    enum_<RootSystem::ScalarTypes>("ScalarType")
//            .value("type", RootSystem::ScalarTypes::st_type)
//            .value("radius", RootSystem::ScalarTypes::st_radius)
//            .value("order", RootSystem::ScalarTypes::st_order)
//            .value("time", RootSystem::ScalarTypes::st_time)
//            .value("length", RootSystem::ScalarTypes::st_length)
//            .value("volume", RootSystem::ScalarTypes::st_volume)
//            .value("surface", RootSystem::ScalarTypes::st_surface)
//            .value("one", RootSystem::ScalarTypes::st_one)
//            .value("userdata1", RootSystem::ScalarTypes::st_userdata1)
//            .value("userdata2", RootSystem::ScalarTypes::st_userdata2)
//            .value("userdata3", RootSystem::ScalarTypes::st_userdata3)
//            .value("parenttype", RootSystem::ScalarTypes::st_parenttype)
//            .value("lb", RootSystem::ScalarTypes::st_lb)
//            .value("la", RootSystem::ScalarTypes::st_la)
//            .value("nob", RootSystem::ScalarTypes::st_nob)
//            .value("r", RootSystem::ScalarTypes::st_r)
//            .value("theta", RootSystem::ScalarTypes::st_theta)
//            .value("rlt", RootSystem::ScalarTypes::st_rlt)
//            .value("meanln", RootSystem::ScalarTypes::st_meanln)
//            .value("sdln", RootSystem::ScalarTypes::st_sdln)
//            ;
//    /*
//     * analysis.h
//     */
//    class_<SegmentAnalyser, SegmentAnalyser*>("SegmentAnalyser")
//        .def(init<RootSystem&>())
//        .def(init<SegmentAnalyser&>())
//        .def("addSegments",addSegments1)
//        .def("addSegments",addSegments2)
//        .def("crop", &SegmentAnalyser::crop)
//        .def("filter", filter1)
//        .def("filter", filter2)
//        .def("pack", &SegmentAnalyser::pack)
//        .def("getScalar", &SegmentAnalyser::getScalar)
//        .def("getSegmentLength", &SegmentAnalyser::getSegmentLength)
//        .def("getSummed", getSummed1)
//        .def("getSummed", getSummed2)
//        .def("distribution", distribution_1)
//        .def("distribution", distribution_2)
//        .def("distribution2", distribution2_1)
//        .def("distribution2", distribution2_2)
//        .def("getRoots", &SegmentAnalyser::getRoots)
//        .def("getNumberOfRoots", &SegmentAnalyser::getNumberOfRoots)
//        .def("cut", cut1)
//        .def("addUserData", &SegmentAnalyser::addUserData)
//        .def("clearUserData", &SegmentAnalyser::clearUserData)
//        .def("write", &SegmentAnalyser::write)
//        .def_readwrite("nodes", &SegmentAnalyser::nodes)
//        .def_readwrite("segments", &SegmentAnalyser::segments)
//        .def_readwrite("ctimes", &SegmentAnalyser::ctimes)
//        // .def("cut", cut2) // not working, see top definition of cut2
//        ;
//    class_<std::vector<SegmentAnalyser>>("std_vector_SegmentAnalyser_")
//            .def(vector_indexing_suite<std::vector<SegmentAnalyser>>() )
//            ;
//    /*
//     * exudation.h
//     */
//    class_<ExudationModel, ExudationModel*>("ExudationModel", init<double, double, int, RootSystem&>())
//            .def(init<double, double, double, int, int, int, RootSystem&>())
//		    .def("calculate", &ExudationModel::calculate)
//		    .def_readwrite("Q", &ExudationModel::Q)
//		    .def_readwrite("Dl", &ExudationModel::Dl)
//		    .def_readwrite("theta", &ExudationModel::theta)
//		    .def_readwrite("R", &ExudationModel::R)
//		    .def_readwrite("k", &ExudationModel::k)
//		    .def_readwrite("l", &ExudationModel::l)
//		    .def_readwrite("type", &ExudationModel::type)
//		    .def_readwrite("n0", &ExudationModel::n0)
//		    .def_readwrite("thresh13", &ExudationModel::thresh13)
//            .def_readwrite("calc13", &ExudationModel::calc13)
//            .def_readwrite("observationRadius", &ExudationModel::observationRadius)
//            ;
//    enum_<ExudationModel::IntegrationType>("IntegrationType")
//            .value("mps_straight", ExudationModel::IntegrationType::mps_straight)
//            .value("mps", ExudationModel::IntegrationType::mps)
//            .value("mls", ExudationModel::IntegrationType::mls)
//            ;

}


} // end namespace CRootBox

#endif /* PYTHON_WRAPPER */

#endif /* PY_ROOTBOX_H_ */
