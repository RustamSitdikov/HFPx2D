cmake_minimum_required(VERSION 3.6)
project(HFPx2D)
set(CMAKE_CXX_STANDARD 11)

IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra -pedantic")
ENDIF(CMAKE_BUILD_TYPE MATCHES DEBUG)

set(SOURCE_FILES_DIR ${PROJECT_SOURCE_DIR}/src)

set(SOURCE_FILES
        main.cpp
        il/io/toml.cpp
        ${SOURCE_FILES_DIR}/Elasticity2D.cpp
        ${SOURCE_FILES_DIR}/Elasticity2D.h
        ${SOURCE_FILES_DIR}/Mesh.cpp
        ${SOURCE_FILES_DIR}/Mesh.h
        ${SOURCE_FILES_DIR}/DOF_Handles.h
        ${SOURCE_FILES_DIR}/DOF_Handles.cpp
        ${SOURCE_FILES_DIR}/AssemblyDDM.cpp
        ${SOURCE_FILES_DIR}/AssemblyDDM.h
        ${SOURCE_FILES_DIR}/Friction.cpp
        ${SOURCE_FILES_DIR}/Friction.h
        ${SOURCE_FILES_DIR}/Dilatancy.cpp
        ${SOURCE_FILES_DIR}/Dilatancy.h
        ${SOURCE_FILES_DIR}/FromEdgeToCol.cpp
        ${SOURCE_FILES_DIR}/FromEdgeToCol.h
        ${SOURCE_FILES_DIR}/FVM.cpp
        ${SOURCE_FILES_DIR}/FVM.h
        ${SOURCE_FILES_DIR}/ConductivitiesNewtonian.cpp
        ${SOURCE_FILES_DIR}/ConductivitiesNewtonian.h
        ${SOURCE_FILES_DIR}/TimeIncr.cpp
        ${SOURCE_FILES_DIR}/TimeIncr.h
        ${SOURCE_FILES_DIR}/Output_results.h
        ${SOURCE_FILES_DIR}/MC_criterion.cpp
        ${SOURCE_FILES_DIR}/MC_criterion.h
        ${SOURCE_FILES_DIR}/EHLDs.cpp
        ${SOURCE_FILES_DIR}/EHLDs.h
        ${SOURCE_FILES_DIR}/Permeability.cpp
        ${SOURCE_FILES_DIR}/Permeability.h
        InputData.toml
        )

link_directories($ENV{MKLROOT}/lib $ENV{MKLROOT}/../compiler/lib)

add_executable(HFPx2D ${SOURCE_FILES})

target_compile_definitions(HFPx2D PRIVATE IL_MKL=1 IL_BLAS=1)

target_include_directories(HFPx2D PRIVATE ${CMAKE_SOURCE_DIR} $ENV{MKLROOT}/include)

target_link_libraries(HFPx2D mkl_intel_lp64 mkl_core mkl_intel_thread dl iomp5 pthread m)
