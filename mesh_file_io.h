//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 1/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef INC_3D_BEM_LOAD_MESH_FROM_FILE_H
#define INC_3D_BEM_LOAD_MESH_FROM_FILE_H

#include <cstdio>
#include <complex>
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

namespace hfp3d {

    template <typename T>
    void save_data_to_csv(
            const il::Array<T> &vector,
            const std::string &trg_dir,
            const std::string &of_name) {
        std::string f_path = trg_dir + of_name;
        const char* format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi\n";
        } else {
            format = "%.16g\n";
        }
        FILE* of=std::fopen(f_path.c_str(),"w");
        for (int j=0; j < vector.size(); ++j){
            std::fprintf(of, format, vector[j]);
        }
        std::fclose(of);
    }

    template <typename T, il::int_t m>
    void save_data_to_csv(
            const il::StaticArray<T, m> &vector,
            const std::string &trg_dir,
            const std::string &of_name) {
        std::string f_path = trg_dir + of_name;
        const char* format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi\n";
        } else {
            format = "%.16g\n";
        }
        FILE* of=std::fopen(f_path.c_str(),"w");
        for (int j=0; j < vector.size(); ++j){
            T out = vector[j];
            std::fprintf(of, format, out);
        }
        std::fclose(of);
    }

    template <typename T>
    void save_data_to_csv(
            const il::Array2D<T> &matrix,
            const std::string &trg_dir,
            const std::string &of_name) {
        std::string f_path = trg_dir + of_name;
        const char* format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi";
        } else {
            format = "%.16g";
        }
        FILE* of=std::fopen(f_path.c_str(),"w");
        for (int j=0; j < matrix.size(0); ++j) {
            for (int k=0; k < matrix.size(1); ++k) {
                T out = matrix(j, k);
                std::fprintf(of, format, out);
                if (k < matrix.size(1)-1) std::fprintf(of, ",");
            }
            std::fprintf(of, "\n");
        }
        std::fclose(of);
    }

    template <class T, il::int_t m, il::int_t n>
    void save_data_to_csv(
            const il::StaticArray2D<T, m, n> &matrix,
            const std::string &trg_dir,
            const std::string &of_name) {
        std::string f_path = trg_dir + of_name;
        const char* format;
        if (sizeof(T) == sizeof(std::complex<double>)) {
            format = "%.16g%+.16gi";
            } else {
            format = "%.16g";
        }
        FILE* of=std::fopen(f_path.c_str(),"w");
        for (int j=0; j < matrix.size(0); ++j) {
            for (int k=0; k < matrix.size(1); ++k) {
                T out = matrix(j, k);
                std::fprintf(of, format, out);
                if (k < matrix.size(1)-1) std::fprintf(of, ",");
            }
            std::fprintf(of, "\n");
        }
        std::fclose(of);
    }

}

#endif //INC_3D_BEM_LOAD_MESH_FROM_FILE_H
