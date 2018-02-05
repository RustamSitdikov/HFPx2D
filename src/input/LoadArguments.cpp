//
// HFPx2D project.
//
// Created by Federico Ciardo on 02.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

// Inclusion from the project
#include "LoadArguments.h"

namespace hfp2d {

void loadArguments(int argc, char const *argv[], il::io_t,
                   il::String &config_filename,
                   il::String &path_output_directory) {
  // loadArguments will receive the arguments that are passed to the program
  // from the command line as input

  // Initial check
  if (argc != 3) {
    std::cerr << "ERROR -> Usage of the program is " << std::endl;
    std::cerr << "  HFPx2D input_file path_output_directory " << std::endl;
    std::cerr << "-- Press ENTER to exit...";
    std::cin.get();
    exit(EXIT_FAILURE);
  }

  // Instantiate an object from stringstream class
  std::stringstream arg_stream;

  // Insert arguments in the object
  for (il::int_t i = 0; i < argc; ++i) arg_stream << argv[i] << " ";

  // Extract arguments from arg_stream object
  std::string first_arg;
  std::string second_arg;
  std::string third_arg;
  arg_stream >> first_arg;
  arg_stream >> second_arg;
  arg_stream >> third_arg;

  il::int_t status_of_access;
  config_filename =
      il::String(il::StringType::Ascii, second_arg.c_str(), second_arg.size());

  status_of_access = access(config_filename.asCString(), F_OK | R_OK);

  if (status_of_access != 0) {
    std::cerr << "Error: impossible to access the input file "
              << config_filename << std::endl;
    std::cerr << strerror(errno) << std::endl;
    exit(EXIT_FAILURE);
  }

  path_output_directory =
      il::String(il::StringType::Ascii, third_arg.c_str(), third_arg.size());

  const int output_dir_status = mkdir(path_output_directory.asCString(), 0777);

  if (output_dir_status != 0) {
    status_of_access =
        access(path_output_directory.asCString(), F_OK | W_OK | R_OK);
  }

  if (status_of_access != 0) {
    std::cerr << "Error: impossible to create or access the directory "
              << path_output_directory << std::endl;
    std::cerr << strerror(errno) << std::endl;
    exit(EXIT_FAILURE);

  } else {
    cleanOutputDir(path_output_directory.asCString());
  }
}

///// This script removes the files in *path that have only some particular
/// extensions
void cleanOutputDir(const char *path) {
  // open the stream of files in the directory *path
  DIR *dir_file = opendir(path);

  // if we have a successfull stream
  if (dir_file != nullptr) {
    // create the pointer to the list of files
    struct dirent *h_file;

    // for each file in the folder
    while ((h_file = readdir(dir_file)) != nullptr) {
      // Avoid dealing with current and parent folders
      if (!(std::string(h_file->d_name).compare("."))) continue;
      if (!(std::string(h_file->d_name).compare(".."))) continue;

      // in linux hidden files all start with '.' --> Ignore them or use a
      // boolean variable such as gIgnoreHidden
      // which has to be set from the outside.
      // if ( gIgnoreHidden && ( h_file->d_name[0] == '.' )) continue;
      if (h_file->d_name[0] == '.') continue;

      // now we compare the extensions
      // save the file as string
      std::string filename = h_file->d_name;
      // find location of the dot for the extension
      size_t dot = filename.find_last_of('.');  // '.' is a char, "." is a
      // string (that is a sequence of char, more expensive)

      // if the dot is found before the extension
      if (dot != std::string::npos)  // npos is returned by find_last_of if the
                                     // dot is not found
      {
        // compare the file extension with the ones to be removed
        std::string ext = filename.substr(dot, filename.size() - dot);
        if ((ext == ".hf2o") || (ext == ".txt")) {
          // Recover the full path for the file to delete
          const std::string full_file_name = std::string(path) + "/" + filename;

          // Delete the file
          int removeStatus = remove(full_file_name.c_str());

          // Check successful delete
          if (removeStatus != 0) {
            std::cerr << "Error while trying to delete " << full_file_name
                      << std::endl;
            std::cerr << strerror(errno) << std::endl;
            exit(1);
          }
        }

      }  // else it is a folder or a file without extension and we do not
         // touch it
    }
    // Close the stream
    closedir(dir_file);
  } else {
    std::cerr << "Error while trying to open directory " << path
              << " for cleaning." << std::endl;
    std::cerr << strerror(errno) << std::endl;
    exit(1);
  }
}
}
