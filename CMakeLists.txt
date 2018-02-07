cmake_minimum_required(VERSION 3.6)

project(HFPx2DUnitTest CXX)

project(HFPx2D CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11")

################################################################################
# choose math & parallel framework
################################################################################

set(IL_MKL 0)
set(IL_OPENBLAS 1)
set(IL_OPENMP 1)

################################################################################
# custom command(s) for some under mac os sierra
################################################################################

set(CUSTOM 1)

################################################################################
# choose compiler: clang, gcc, ???
################################################################################

if(APPLE)
    set(CMAKE_C_COMPILER /usr/bin/clang)
    set(CMAKE_CXX_COMPILER /usr/bin/clang++)
elseif(UNIX)
    set(CMAKE_C_COMPILER /usr/bin/gcc)
    set(CMAKE_CXX_COMPILER /usr/bin/g++)
elseif(WIN32)
    if(MINGW)
        set(CMAKE_C_COMPILER C:/msys64/mingw64/bin/gcc)
        set(CMAKE_CXX_COMPILER C:/msys64/mingw64/bin/g++)
    endif()
endif()

################################################################################
#  treat warnings as errors in Debug mode
################################################################################

IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        #        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        #        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -pedantic")
    endif()
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

################################################################################
# parallel framework
################################################################################

# For OpenMP
if(IL_OPENMP)
    if(UNIX)
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            # link_directories(/opt/clang-3.9.0/lib)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        endif()
    elseif(APPLE)
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
            # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            # link_directories(/opt/clang-3.9.0/lib)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif()
    elseif(WIN32)
        if(MINGW)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif()
    endif()
endif()

# For TBB
if(IL_TBB)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_TBB")
    # target_link_libraries(InsideLoop tbb)
endif()

# For Cilk
if(IL_CILK)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_CILK")
endif()

################################################################################
# math framework
################################################################################

# For Intel MKL
if(IL_MKL)
    if(UNIX)
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
            link_directories(/opt/intel/mkl/lib/intel64/opt/intel/compilers_and_libraries/linux/lib/intel64)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
            link_directories($ENV{MKLROOT}/lib/intel64_lin)
            link_directories(/opt/intel/mkl/lib/intel64/opt/intel/compilers_and_libraries/linux/lib/intel64)
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel -DIL_MKL -DIL_BLAS")
            link_directories(/opt/intel/mkl/lib/intel64/opt/intel/compilers_and_libraries/linux/lib/intel64)
        endif()
    elseif(APPLE)
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
            link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib )
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -DIL_MKL -DIL_BLAS")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")
            set(CMAKE_MKL_LIBRARIES mkl_intel_lp64 mkl_core mkl_intel_thread iomp5 dl pthread m)
            link_directories($ENV{MKLROOT}/lib/intel64_lin)
            link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib )
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel -DIL_MKL -DIL_BLAS")
        endif()
    elseif(WIN32) # to fill up here ...
        if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
            # to fill up here ...
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            # to fill up here ...
        elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            # to fill up here ...
        endif()
    endif()
# For OpenBlas
elseif(IL_OPENBLAS)
    if(UNIX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
        set(CMAKE_OPENBLAS_LIBRARIES "openblas")
    elseif(WIN32)
        if(MINGW)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_OPENBLAS -DIL_BLAS")
            set(CMAKE_OPENBLAS_LIBRARIES "openblas")
        endif()
    elseif(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_MKL -DIL_OPENBLAS -DIL_BLAS")
        set(CMAKE_OPENBLAS_LIBRARIES "openblas")
    endif()
endif()

################################################################################
#  FOR HFPx2D
################################################################################

# set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

set(SOURCE_FILES
        src/core/Mesh.cpp src/core/Mesh.h
        src/core/SegmentData.h
        src/core/Fluid.h
        src/core/ElasticProperties.h
        src/core_dev/SolidEvolution.h
        src/core/Sources.h
        src/core/DomainMesh.h
        src/core/SimulationParameters.h
        src/core/Solution.h
        src/core/SolidProperties.h
        src/core/InSituConditions.h
        src/elasticity/AssemblyDDM.cpp src/elasticity/AssemblyDDM.h
        src/elasticity/Simplified3D.cpp src/elasticity/Simplified3D.h
        src/elasticity/PlaneStrainInfinite.cpp src/elasticity/PlaneStrainInfinite.h
        src/solvers/SimpleElasticBenchmarks.cpp src/solvers/SimpleElasticBenchmarks.h
        #        src/input/loadArguments.cpp
        #        src/input/loadArguments.h
        #        src/input/loadInput.cpp
        #        src/input/loadInput.h
        #        src/input/loadInput.cpp
        #        src/input/loadInput.h
        #        src/input/geometry/autoLineMeshUtilities.cpp
        #        src/input/geometry/autoLineMeshUtilities.h
        #        src/input/geometry/autoLineMeshCreation.cpp
        #        src/input/geometry/autoLineMeshCreation.h
        #        src/input/geometry/loadGeometry.cpp
        #        src/input/geometry/loadGeometry.h
        #        src/input/geometry/loadMeshFile.cpp
        #        src/input/geometry/loadMeshFile.h
        ##        src/input/properties/loadProperties.cpp
        #        src/input/properties/loadProperties.h
        #src/input/geometry/autoLineMeshInfo.cpp
        #src/input/geometry/autoLineMeshInfo.h
        #        src/input/properties/loadSingleMaterial.cpp
        #        src/input/properties/loadSingleMaterial.h
        src/input/findUtilities.cpp src/input/findUtilities.h
        il/io/toml/toml.h il/io/toml/toml.cpp
        #        src/input/conditions/loadConditions.cpp
        #        src/input/conditions/loadConditions.h
        #        src/input/Sources/loadSources.cpp
        #        src/input/Sources/loadSources.h
        src/ehlsolvers/ReynoldsP0.cpp src/ehlsolvers/ReynoldsP0.h
        src/tip/tipAsymptote.cpp src/tip/tipAsymptote.h
        src/util/json.hpp
        src/input/json/loadJsonMesh.cpp src/input/json/loadJsonMesh.h
        src/wellbore/WellInjection.cpp src/wellbore/WellInjection.h
        src/wellbore/WellMesh.h
        src/wellbore/WellFlowP0.cpp src/wellbore/WellFlowP0.h
        src/wellbore/WellSolution.h
        src/wellbore/SimpleWellFlowSolverBenchmark.cpp
        src/wellbore/SimpleWellFlowSolverBenchmark.h
        src/solvers/MultiFracsSolution.h
        src/input/json/LoadInputMultistage.cpp src/input/json/LoadInputMultistage.h
        src/solvers/HFPropagationP0.cpp src/solvers/HFPropagationP0.h
        src/solvers/MultiFracsSolver.cpp src/solvers/MultiFracsSolver.h
        src/util/RootFinder.cpp src/util/RootFinder.h)

add_executable(HFPx2D ${SOURCE_FILES}  main.cpp)
target_link_libraries(HFPx2D ${CMAKE_MKL_LIBRARIES} ${CMAKE_OPENBLAS_LIBRARIES})

if(IL_MKL)
    if(APPLE)
        target_compile_definitions(HFPx2D PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
    elseif(UNIX)
        target_compile_definitions(HFPx2D PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} /opt/intel/mkl/include)
    elseif(WIN32)
        target_compile_definitions(HFPx2D PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)
    endif()
elseif(IL_OPENBLAS)
    target_compile_definitions(HFPx2D PRIVATE IL_OPENBLAS=1 IL_BLAS=1 IL_OPENMP=1)
    target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR})
endif()

################################################################################
# below for some under mac os sierra ;(
################################################################################

if(APPLE AND CUSTOM)
    add_custom_command(TARGET HFPx2D
        POST_BUILD COMMAND
        /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
        $<TARGET_FILE:HFPx2D>)
endif()


################################################################################
# For Unit Tests  HFPx2DUnitTest
################################################################################

set(SR_FOLDER "${CMAKE_SOURCE_DIR}/src")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_FOLDER='\"${SR_FOLDER}\"'")

set(UNIT_TEST_FILES src/_test/elastic_benchmarks.cpp
        gtest/src/gtest-all.cc
        src/_test/mesh_test.cpp
        src/tip/_test/tip_test.cpp
        src/_test/domain_mesh_test.cpp
        src/wellbore/SimpleWellFlowSolverBenchmark.cpp
        src/wellbore/SimpleWellFlowSolverBenchmark.h
        src/elasticity/_test/elasticity_tests.cpp)

add_executable(HFPx2DUnitTest ${SOURCE_FILES} ${UNIT_TEST_FILES} test.cpp)

target_link_libraries(HFPx2DUnitTest ${CMAKE_MKL_LIBRARIES} ${CMAKE_OPENBLAS_LIBRARIES})

if(IL_MKL)
    if(APPLE)
        target_compile_definitions(HFPx2DUnitTest PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(HFPx2DUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/gtest $ENV{MKLROOT}/include)
    elseif(UNIX)
        target_compile_definitions(HFPx2DUnitTest PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(HFPx2DUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/gtest /opt/intel/mkl/include)
    elseif(WIN32)
        target_compile_definitions(HFPx2DUnitTest PRIVATE IL_MKL=1 IL_BLAS=1 IL_OPENMP=1)
        target_include_directories(HFPx2DUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/gtest $ENV{MKLROOT}/include)
    endif()
elseif(IL_OPENBLAS)
    target_compile_definitions(HFPx2DUnitTest PRIVATE IL_OPENBLAS=1 IL_BLAS=1 IL_OPENMP=1)
    target_include_directories(HFPx2DUnitTest PUBLIC ${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/gtest)
endif()

# For unit tests: The precondition of our fonctions are checked with assert
# macros that terminate the program in debug mode. In order to test those macros
# in our unit tests, the flag IL_UNIT_TEST turns them into exceptions that can
# be caught and checked.
target_compile_options(HFPx2DUnitTest PRIVATE "-DIL_UNIT_TEST")

################################################################################
# below for some under mac os sierra ;(
################################################################################

if(APPLE AND CUSTOM)
    add_custom_command(TARGET HFPx2DUnitTest
        POST_BUILD COMMAND
        /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
        $<TARGET_FILE:HFPx2DUnitTest>)
endif()