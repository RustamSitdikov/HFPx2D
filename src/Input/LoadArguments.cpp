//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 02.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include "LoadArguments.h"

namespace hfp2d {

void loadArguments(int argc, char const *argv[], il::io_t, bool check_input,
                   bool check_output, bool check_restart,
                   il::String &input_filename, il::String &restart_filename,
                   il::String &path_output_directory) {
  // loadArguments will receive the arguments that were passed to the program
  // from the command line as Input.
  // they are organized as following:
  // .. argc = argument count, the number of arguments;
  // .. argv = a set of characters which are the Input arguments.
  //
  // To determine which argument is which, we will use the switches -i, -o, -r

  if (argc < 5) {
    std::cerr << "ERROR -> Usage of the program is " << std::endl;
    std::cerr << "  HFPx2D -i input_file -o path_output_directory     for new "
                 "analyses;"
              << std::endl;
    std::cerr << "-- Press ENTER to exit...";
    std::cin.get();
    exit(EXIT_FAILURE);
  }

  // create an object from stringstream class
  std::stringstream arg_stream;

  // Insert arguments in the object
  for (int i = 0; i < argc; ++i) arg_stream << argv[i] << " ";

  std::string first_arg, second_arg, third_arg, fourth_arg, fifth_arg;
  arg_stream >> first_arg;
  arg_stream >> second_arg;
  arg_stream >> third_arg;
  arg_stream >> fourth_arg;
  arg_stream >> fifth_arg;

  int status_of_access;
  if (second_arg == "-i") {
    input_filename =
        il::String(il::StringType::Bytes, third_arg.c_str(), third_arg.size());

    status_of_access = access(input_filename.asCString(), F_OK | R_OK);

    if (status_of_access != 0) {
      std::cerr << "Error: impossible to access the Input file "
                << input_filename << std::endl;
      std::cerr << strerror(errno) << std::endl;
      exit(EXIT_FAILURE);
    } else {
      check_input = true;
    }

  } else if (second_arg == "-r") {
    restart_filename =
        il::String(il::StringType::Bytes, third_arg.c_str(), third_arg.size());
    status_of_access = access(restart_filename.asCString(), F_OK | R_OK);

    if (status_of_access != 0) {
      std::cerr << "Error: impossible to access the Input file "
                << restart_filename << std::endl;
      std::cerr << strerror(errno) << std::endl;
      exit(EXIT_FAILURE);
    } else {
      check_restart = true;
    }

  } else {
    std::cout << "Wrong argument for Input or restart analysis. Please change "
                 "second argument."
              << std::endl;
    std::cout << "-- Press ENTER to exit...";
    std::cin.get();
    exit(1);
  }

  if (fourth_arg == "-o") {
    path_output_directory =
        il::String(il::StringType::Bytes, fifth_arg.c_str(), fifth_arg.size());
    check_output = true;
  } else {
    std::cout << "Wrong fourth argument for output directory. Please change "
                 "fourth argument."
              << std::endl;
    std::cout << "-- Press ENTER to exit...";
    std::cin.get();
    exit(1);
  }

  /// Creating the OUTPUT directory

  if (check_output) {
    int output_dir_status = mkdir(path_output_directory.asCString(), 0777);
    if (output_dir_status != 0) {
      status_of_access =
          access(path_output_directory.asCString(), F_OK | W_OK | R_OK);
    }

    if (status_of_access !=
        0) {  // if there is an error on access as well, no option left:
      // we have to stop because there is no way to work on the folder
      std::cerr << "Error: impossible to create or access the directory "
                << path_output_directory << std::endl;
      std::cerr << strerror(errno) << std::endl;
      exit(EXIT_FAILURE);

    } else {
      if (!check_restart) {
        // The directory is accessible and we can work in it.
        // Let us clean the directory ONLY IF THE ANALYSIS IS NEW !!!
        cleanOutputDir(path_output_directory.asCString());
      }
    }
  }
}

///// This script removes the files in *path that have only some particular
/// extensions
// TODO: clean up, change to the il::String form and check that is still
// working properly
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

}  // End of script
}
