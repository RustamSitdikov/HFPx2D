cmake_minimum_required(VERSION 3.6)


project(HFPx2DUnitTest CXX)

project(HFPx2D CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")

# If you want to run another compiler, such as clang
if (APPLE)
set(CMAKE_C_COMPILER /usr/bin/clang)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)
elseif (UNIX)
set(CMAKE_C_COMPILER /usr/bin/gcc)
set(CMAKE_CXX_COMPILER /usr/bin/g++)

    elseif(WIN32) # to fill up here
endif()

#  FOR HFPx2D
#  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

set(SOURCE_FILES
        src/core/Mesh.cpp
        src/core/Mesh.h
        src/core/Utilities.h
        src/core/Utilities.cpp
        src/core/SegmentData.h
        src/core/FluidProperties.h
        src/core/ElasticProperties.h
        src/core_dev/SolidEvolution.h
        src/core/Sources.h
        src/elasticity/AssemblyDDM.cpp
        src/elasticity/AssemblyDDM.h
        src/elasticity/Simplified3D.cpp
        src/elasticity/Simplified3D.h
        src/elasticity/PlaneStrainInfinite.cpp
        src/elasticity/PlaneStrainInfinite.h
        src/solvers/SimpleElastic.cpp
        src/solvers/SimpleElastic.h
        src/core/SimulationParameters.h
        src/input/LoadArguments.cpp
        src/input/LoadArguments.h
        src/input/loadInput.cpp
        src/input/loadInput.h
        src/input/Geometry/AutoLineMeshUtilities.cpp
        src/input/Geometry/AutoLineMeshUtilities.h
        src/input/Geometry/AutoLineMeshCreation.cpp
        src/input/Geometry/AutoLineMeshCreation.h
        src/input/Geometry/LoadGeometry.cpp
        src/input/Geometry/LoadGeometry.h
        src/input/Geometry/LoadMeshFile.cpp
        src/input/Geometry/LoadMeshFile.h
        src/input/Properties/LoadProperties.cpp
        src/input/Properties/LoadProperties.h
        src/input/Geometry/AutoLineMeshInfo.cpp
        src/input/Geometry/AutoLineMeshInfo.h
        src/input/Properties/LoadSingleMaterial.cpp
        src/input/Properties/LoadSingleMaterial.h
        src/input/FindUtilities.cpp
        src/input/FindUtilities.h
        il/io/toml/toml.h
        il/io/toml/toml.cpp
        src/input/Conditions/LoadConditions.cpp
        src/input/Conditions/LoadConditions.h
        src/core/Solution.h
        src/core/SolidProperties.h
        src/tip/TipAsymptote.cpp
        src/tip/TipAsymptote.h
        src/solvers/FluidInjFrictWeakDilatFault.cpp
        src/solvers/FluidInjFrictWeakDilatFault.h
        src/core_dev/FractureEvolution.h
        src/core/InSituStress.h
        src/devt/FromEdgeToCol.cpp
        src/devt/FromEdgeToCol.h
        src/elhsolvers/ReynoldsP1.cpp
        src/elhsolvers/ReynoldsP1.h
        src/devt/FiniteVolumeRoutines.cpp
        src/devt/FiniteVolumeRoutines.h
        src/util/json.hpp
        src/devt/IterativeSolvers.cpp
        src/devt/IterativeSolvers.h)

if(APPLE)
link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib )
elseif(UNIX)
    link_directories(/opt/intel/mkl/lib/intel64
                     /opt/intel/compilers_and_libraries/linux/lib/intel64)
endif()

add_executable(HFPx2D ${SOURCE_FILES}  main.cpp)

target_compile_definitions(HFPx2D PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)

target_link_libraries(HFPx2D mkl_intel_lp64 mkl_core mkl_intel_thread dl iomp5 pthread m)

if(APPLE)
    target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
elseif(UNIX)
    target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} /opt/intel/mkl/include)
endif()

#if(APPLE)
# below for some under mac os sierra ;(
#add_custom_command(TARGET HFPx2D
        #POST_BUILD COMMAND
        #/usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
        #$<TARGET_FILE:HFPx2D>)
#endif()


################################################################################
# For Unit Tests  HFPx2DUnitTest
################################################################################

set(SR_FOLDER "${CMAKE_SOURCE_DIR}/src")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_FOLDER='\"${SR_FOLDER}\"'")

set(UNIT_TEST_FILES
        src/_test/elastic_test.cpp
        gtest/src/gtest-all.cc
        src/_test/mesh_test.cpp
        src/tip/_test/tip_test.cpp
        src/_test/FV_test.cpp
        src/elasticity/_test/elasticity_tests.cpp)

add_executable(HFPx2DUnitTest ${SOURCE_FILES} ${UNIT_TEST_FILES} test.cpp)

target_compile_definitions(HFPx2DUnitTest PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)

target_link_libraries(HFPx2DUnitTest mkl_intel_lp64 mkl_core mkl_intel_thread dl iomp5 pthread m)

# should be the same -> make sure the MKLROOT is defined under linux
if(APPLE)
    target_include_directories(HFPx2DUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/gtest $ENV{MKLROOT}/include)
elseif(UNIX)
    target_include_directories(HFPx2DUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/gtest  /opt/intel/mkl/include)
endif()


# For unit tests: The precondition of our fonctions are checked with assert
# macros that terminate the program in debug mode. In order to test those macros
# in our unit tests, the flag IL_UNIT_TEST turns them into exceptions that can
# be caught and checked.
target_compile_options(HFPx2DUnitTest PRIVATE "-DIL_UNIT_TEST")

#if(APPLE)
    # below for some under mac os sierra ;(
#    add_custom_command(TARGET HFPx2DUnitTest
#        POST_BUILD COMMAND
#        /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
#        $<TARGET_FILE:HFPx2DUnitTest>)
#endif()
