#ifndef PY_ROOTBOX_H_
#define PY_ROOTBOX_H_


/**
 *  A Python module for CRootbox based on boost.python
 *
 *  build shared library
 *  1. export LD_LIBRARY_PATH=~/boost_1_62_0/stage/lib
 *  2. g++ -std=c++11 -O3 -fpic -shared -o py_rootbox.so -Wl,-soname,"py_rootbox.so" PythonRootSystem.cpp -I/usr/include/python3.5 -L/home/daniel/boost_1_62_0/stage/lib -lboost_python Debug/ModelParameter.o Debug/Root.o Debug/RootSystem.o Debug/analysis.o Debug/sdf.o Debug/tropism.o
 *
 * sdf.h 		writePVPScript is not exposed, use RootSystem::write to write the geometry script
 * mymath.h		currently only Vector3d is exposed (lets see if we will need anything else)
 *
 */
// #define PYTHON_WRAPPER // UNCOMMENT TO BUILD SHARED LIBRARY

#ifdef PYTHON_WRAPPER

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "mymath.h"
#include "sdf.h"
#include "RootSystem.h"
#include "analysis.h"

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
double (AnalysisSDF::*getSummed1)(int st) const = &AnalysisSDF::getSummed;
double (AnalysisSDF::*getSummed2)(int st, SignedDistanceFunction* geometry) const = &AnalysisSDF::getSummed;

/**
 * Virtual functions (not sure if needed, or only if we derive classes from it in pyhton?), not working...
 *
 * it seems a bit tricky to make polymorphism work, we have to wrap the classes
 *
 * Tutorial example:
 * struct BaseWrap : Base, wrapper<Base>
 * {
 *    int f()
 *    {
 *        if (override f = this->get_override("f"))
 *            return f(); // *note*
 *        return Base::f();
 *    }
 *
 *    int default_f() { return this->Base::f(); }
 * };
 */
//struct SignedDistanceFunction_Wrap : SignedDistanceFunction, wrapper<SignedDistanceFunction>
//{
//	double getDist(const Vector3d& v) const {
//		if (override getDist = this->get_override("getDist"))
//			return getDist(v);
//		return SignedDistanceFunction::getDist(v);
//	}
//	double default_getDist(const Vector3d& v) { return this->SignedDistanceFunction::getDist(v); }
//};



/**
 * Expose classes to Pyhton module
 */
BOOST_PYTHON_MODULE(py_rootbox)
{
    /* general */
    class_<std::vector<double>>("std_vector_double_")
        .def(vector_indexing_suite<std::vector<double>>() );
    class_<std::vector<int>>("std_vector_int_")
        .def(vector_indexing_suite<std::vector<int>>() );
	/* mymath.h */
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
	/* sdf.h */
//	class_<SignedDistanceFunction_Wrap, boost::noncopyable>("SignedDistanceFunction")
//	    .def("getDist", &SignedDistanceFunction_Wrap::getDist, &SignedDistanceFunction_Wrap::default_getDist)
//	; // TODO how does polymorphism work...
	class_<SignedDistanceFunction>("SignedDistanceFunction")
			.def("getDist",&SignedDistanceFunction::getDist)
			.def("__str__",&SignedDistanceFunction::toString)
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
	// TODO SDF_Intersection, SDF_Union, SDF_Difference, SDF_Complement, SDF_HalfPlane

	/* ModelParameter */
	class_<RootTypeParameter>("RootTypeParameter", init<>())
			.def(init<RootTypeParameter&>())
			//.def("set",&RootTypeParameter::set) // why not?
			.def("realize",&RootTypeParameter::realize)
			.def("__str__",&RootTypeParameter::toString)
	;
	class_<RootParameter>("RootParameter", init<>())
			.def(init<int , double, double, const std::vector<double>&, int, double, double, double, double>())
			.def("set",&RootParameter::set)
			.def_readwrite("type", &RootParameter::type)
			.def_readwrite("lb", &RootParameter::lb)
			.def_readwrite("la", &RootParameter::la)
			.def_readwrite("ln", &RootParameter::ln)
			.def_readwrite("nob", &RootParameter::nob)
			.def_readwrite("r", &RootParameter::r)
			.def_readwrite("a", &RootParameter::a)
			.def_readwrite("theta", &RootParameter::theta)
			.def_readwrite("rlt", &RootParameter::rlt)
			.def("getK",&RootParameter::toString)
			.def("__str__",&RootParameter::toString)
	;
	/* RootSystem.h */
        class_<RootSystem>("RootSystem")
        .def("openFile", &RootSystem::openFile)
		.def("setGeometry", &RootSystem::setGeometry)
		.def("setSoil", &RootSystem::setSoil)
		.def("initialize", &RootSystem::initialize)
		.def("simulate",&RootSystem::simulate)
		.def("getNumberOfNodes", &RootSystem::getNumberOfNodes)
		.def("write",&RootSystem::write)
	;
	enum_<RootSystem::ScalarTypes>("ScalarType")
	    .value("type", RootSystem::ScalarTypes::st_type)
	    .value("radius", RootSystem::ScalarTypes::st_radius)
	    .value("order", RootSystem::ScalarTypes::st_order)
	    .value("red", RootSystem::ScalarTypes::st_red)
	    .value("green", RootSystem::ScalarTypes::st_green)
	    .value("blue", RootSystem::ScalarTypes::st_blue)
	    .value("time", RootSystem::ScalarTypes::st_time)
	    .value("length", RootSystem::ScalarTypes::st_length)
	    .value("surface", RootSystem::ScalarTypes::st_surface)
	;
	enum_<RootSystem::OutputTypes>("OutputType")
	    .value("segments", RootSystem::OutputTypes::ot_segments)
	    .value("polylines", RootSystem::OutputTypes::ot_polylines)
	;
    /* analysis.h */
    class_<AnalysisSDF>("AnalysisSDF",init<RootSystem&>()) //
    	.def(init<AnalysisSDF&>())
		.def("pack", &AnalysisSDF::pack)
		.def("getScalar", &AnalysisSDF::getScalar)
		.def("getSummed", getSummed1)
		.def("getSummed", getSummed2)
                .def("getNumberOfRoots", &AnalysisSDF::getNumberOfRoots)
		.def("write",&AnalysisSDF::write)
    ;

}

/*
 *  currently not exposed..
 *
  	class_<Vector2i>("Vector2i", init<>())
			.def(init<int,int>())
			.def(init<Vector2i&>())
			.def_readwrite("x",&Vector2i::x)
			.def_readwrite("y",&Vector2i::y)
			.def("__str__",&Vector2i::toString)
			;
	class_<Vector2d>("Vector2d", init<>())
			.def(init<double,double>())
			.def(init<Vector2d&>())
			.def_readwrite("x",&Vector2d::x)
			.def_readwrite("y",&Vector2d::y)
			.def("__str__",&Vector2d::toString)
	;
 */

/**
 * solution to wrap vectors from Stackoverflow
 *
// C++ code
typedef std::vector<std::string> MyList;
class MyClass {
  MyList myFuncGet();
  void myFuncSet(const Mylist& list);
  //       stuff
};

// Wrapper code

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;


BOOST_PYTHON_MODULE(mymodule)
{
    class_<MyList>("MyList")
        .def(vector_indexing_suite<MyList>() );

    class_<myClass>("MyClass")
        .def("myFuncGet", &MyClass::myFuncGet)
        .def("myFuncSet", &MyClass::myFuncSet)
        ;
}
*
 */

#endif /* PYTHON_WRAPPER */

#endif /* PY_ROOTBOX_H_ */
