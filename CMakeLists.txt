cmake_minimum_required(VERSION 3.6)
project(HFPx2D)
link_directories($ENV{MKLROOT}/lib)
set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES main.cpp Elasticity2D.cpp Elasticity2D.h Mesh.cpp Mesh.h DOF_Handles.h DOF_Handles.cpp AssemblyDDM.cpp AssemblyDDM.h)
add_executable(HFPx2D ${SOURCE_FILES})
target_compile_definitions(HFPx2D PRIVATE IL_MKL=1)
target_include_directories(HFPx2D PRIVATE $ENV{MKLROOT}/include)
target_link_libraries(HFPx2D mkl_intel_lp64 mkl_core mkl_intel_thread dl iomp5 pthread m)
add_custom_command(TARGET HFPx2D
        POST_BUILD COMMAND
        /usr/bin/install_name_tool -add_rpath /opt/intel/compilers_and_libraries_2017.0.102/mac/mkl/lib
        $<TARGET_FILE:HFPx2D>)
