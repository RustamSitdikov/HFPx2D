cmake_minimum_required(VERSION 3.6)

project(HFPx2D CXX)
project(HFPx2DUnitTest CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")


IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

set(SOURCE_FILES
        src/Elasticity/AssemblyDDM.cpp
        src/Elasticity/AssemblyDDM.h
        src/core/Mesh.cpp
        src/core/Mesh.h
        src/core/DOF_Handles.h
        src/core/DOF_Handles.cpp
        src/Elasticity/Simplified3D.cpp
        src/Elasticity/Simplified3D.h
        il/math.h
        src/Elasticity/PlaneStrainInfinite.cpp
        src/Elasticity/PlaneStrainInfinite.h
        src/Solvers/SimpleElastic.cpp src/Solvers/SimpleElastic.h src/core/SolutionClass.h src/core/SimulationParameters.h src/core/loadArguments.cpp src/core/loadArguments.h src/core/loadInput.cpp src/core/loadInput.h)

link_directories($ENV{MKLROOT}/lib/intel64 /opt/intel/compilers_and_libraries/linux/lib/intel64)

add_executable(HFPx2D ${SOURCE_FILES}  main.cpp src/core/ElasticProperties.h)

target_compile_definitions(HFPx2D PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)

target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)

target_link_libraries(HFPx2D mkl_intel_lp64 mkl_core mkl_intel_thread dl iomp5 pthread m)

#add_custom_command(TARGET HFPx2D
#        POST_BUILD COMMAND
#        /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
#        $<TARGET_FILE:HFPx2D>)



################################################################################
# For Unit Tests
################################################################################

set(SR_FOLDER "${CMAKE_SOURCE_DIR}/src")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_FOLDER='\"${SR_FOLDER}\"'")

set(UNIT_TEST_FILES
        src/_test/elastic_test.cpp
        gtest/src/gtest-all.cc
        )

add_executable(HFPx2DUnitTest ${SOURCE_FILES} ${UNIT_TEST_FILES} test.cpp)

target_compile_definitions(HFPx2DUnitTest PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)

target_include_directories(HFPx2DUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/il ${CMAKE_SOURCE_DIR}/gtest $ENV{MKLROOT}/include)

target_link_libraries(HFPx2DUnitTest mkl_intel_lp64 mkl_core mkl_intel_thread dl iomp5 pthread m)


# For unit tests: The precondition of our fonctions are checked with assert
# macros that terminate the program in debug mode. In order to test those macros
# in our unit tests, the flag IL_UNIT_TEST turns them into exceptions that can
# be caught and checked.
target_compile_options(HFPx2DUnitTest PRIVATE "-DIL_UNIT_TEST")

#add_custom_command(TARGET HFPx2DUnitTest
#        POST_BUILD COMMAND
#        /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
#        $<TARGET_FILE:HFPx2DUnitTest>)

