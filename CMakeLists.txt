# build the CRootBox library using CMake
# Simply type cmake . && make
# although out-of-source-builds are encouraged

cmake_minimum_required (VERSION 2.8.12)
project (CRootBox)

# set compiler options
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(RELEASE_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare")
set(DEBUG_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align")
set(RELEASE_OPTS "-std=c++11 -fno-strict-aliasing -fstrict-overflow -fno-finite-math-only -O3 -march=native -funroll-loops -g0")
set(CMAKE_CXX_FLAGS_RELEASE "${RELEASE_WARNING_OPTS} ${RELEASE_OPTS}")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -ggdb -Wall")

# set default build type to release
set(CMAKE_BUILD_TYPE Release)

# add source directory to the include path
include_directories(${PROJECT_SOURCE_DIR})

#
# 1. Make RootBox library
#

# file(GLOB all_sources *.cpp) # i changed that to exclude PythonRootSystem.cpp
add_library(CRootBox
analysis.cpp
ModelParameter.cpp
Root.cpp
RootSystem.cpp
sdf.cpp
tropism.cpp
)

#
# 2. Make py_rootbox library
#

# Find necessary packages for py_rootbox
find_package(PythonLibs 3)
find_package(PythonInterp 3) # creates PYTHON_VERSION_MAJOR, PYTHON_VERSION_MINOR

# Only build the library if Python 3 was found
if (PythonLibs_FOUND)
  include_directories( ${PYTHON_INCLUDE_DIRS} )
  find_package( Boost COMPONENTS python-py${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED )
  include_directories( ${Boost_INCLUDE_DIR} )

  # make shared library for Pyhton binding
  add_library(py_rootbox SHARED
              PythonRootSystem.cpp
              analysis.cpp
              ModelParameter.cpp
              Root.cpp
              RootSystem.cpp
              sdf.cpp
              tropism.cpp
              examples/Exudation/gauss_legendre.cpp)
  set_property(TARGET py_rootbox PROPERTY POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(py_rootbox ${Boost_LIBRARIES} )
  set_target_properties(py_rootbox PROPERTIES PREFIX "" )
else ()
  message(WARNING "Python 3 was not found on your system! You can only use the C++ interface of CRootBox. \
                   Install Python 3 and rerun CMake to build the py_rootbox library.")
endif ()
