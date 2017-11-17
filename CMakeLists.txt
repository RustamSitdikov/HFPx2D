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
        src/Core/Mesh.cpp
        src/Core/Mesh.h
        src/Core/Utilities.h
        src/Core/Utilities.cpp
        src/Core/SegmentData.h
        src/Core/FluidProperties.h
        src/Core/ElasticProperties.h
        src/Core_dev/SolidEvolution.h
        src/Core/Sources.h
        src/Elasticity/AssemblyDDM.cpp
        src/Elasticity/AssemblyDDM.h
        src/Elasticity/Simplified3D.cpp
        src/Elasticity/Simplified3D.h
        src/Elasticity/PlaneStrainInfinite.cpp
        src/Elasticity/PlaneStrainInfinite.h
        src/Solvers/SimpleElastic.cpp
        src/Solvers/SimpleElastic.h
        src/Core/SimulationParameters.h
        src/Input/LoadArguments.cpp
        src/Input/LoadArguments.h
        src/Input/loadInput.cpp
        src/Input/loadInput.h
        src/Input/Geometry/AutoLineMeshUtilities.cpp
        src/Input/Geometry/AutoLineMeshUtilities.h
        src/Input/Geometry/AutoLineMeshCreation.cpp
        src/Input/Geometry/AutoLineMeshCreation.h
        src/Input/Geometry/LoadGeometry.cpp
        src/Input/Geometry/LoadGeometry.h
        src/Input/Geometry/LoadMeshFile.cpp
        src/Input/Geometry/LoadMeshFile.h
        src/Input/Properties/LoadProperties.cpp
        src/Input/Properties/LoadProperties.h
        src/Input/Geometry/AutoLineMeshInfo.cpp
        src/Input/Geometry/AutoLineMeshInfo.h
        src/Input/Properties/LoadSingleMaterial.cpp
        src/Input/Properties/LoadSingleMaterial.h
        src/Input/FindUtilities.cpp
        src/Input/FindUtilities.h
        il/io/toml/toml.h
        il/io/toml/toml.cpp
        src/Input/Conditions/LoadConditions.cpp
        src/Input/Conditions/LoadConditions.h
        src/Input/Sources/LoadSources.cpp
        src/Input/Sources/LoadSources.h
        #src/FractureFluidFlow/ReynoldsP0.cpp
        #src/FractureFluidFlow/ReynoldsP0.h
        #src/Solvers/DevLHFMSolver.cpp
        #src/Solvers/DevLHFMSolver.h
        src/Core/Solution.h
        src/Core/SolidProperties.h
        src/Tip/TipAsymptote.cpp
        src/Tip/TipAsymptote.h
        src/Solvers/FluidInjFrictWeakDilatFault.cpp
        src/Solvers/FluidInjFrictWeakDilatFault.h
        src/Core_dev/FractureEvolution.h
        src/Core/InSituStress.h
        src/Devt/FromEdgeToCol.cpp
        src/Devt/FromEdgeToCol.h
        src/FractureFluidFlow/ReynoldsP1.cpp
        src/FractureFluidFlow/ReynoldsP1.h
        src/Devt/FiniteVolumeRoutines.cpp
        src/Devt/FiniteVolumeRoutines.h)

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

set(UNIT_TEST_FILES src/_test/elastic_test.cpp
        gtest/src/gtest-all.cc
        src/_test/mesh_test.cpp
        src/Tip/_test/tip_test.cpp)

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

if(APPLE)
    # below for some under mac os sierra ;(
    add_custom_command(TARGET HFPx2DUnitTest
        POST_BUILD COMMAND
        /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
        $<TARGET_FILE:HFPx2DUnitTest>)
endif()
