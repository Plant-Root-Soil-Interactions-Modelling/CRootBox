# build the CRootBox library using CMake
# Simply type cmake . && make
# although out-of-source-builds are encouraged

cmake_minimum_required (VERSION 2.8.12)
project (CRootBox)

# set compiler options
set(RELEASE_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align -Wno-sign-compare")
set(DEBUG_WARNING_OPTS "-Wall -Wunused -Wmissing-include-dirs -Wcast-align")
set(RELEASE_OPTS "-std=c++11 -fno-strict-aliasing -fstrict-overflow -fno-finite-math-only -O3 -march=native -funroll-loops -g0")
set(CMAKE_CXX_FLAGS_RELEASE "${RELEASE_WARNING_OPTS} ${RELEASE_OPTS}")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -O0 -ggdb -Wall")

# set default build type to release
set(CMAKE_BUILD_TYPE Release)

# add source directory to the include path
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/external)

# add subdirectories
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(python)
add_subdirectory(lib)
