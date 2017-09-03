//
// HFPx2D project.
//
// Created by Lorenzo Benedetti on 02.09.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#include <iostream>
#include <cstring>
#include <fstream>
//#include <cstdio>
#include <sys/stat.h>
//#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include "loadArguments.h"

namespace hfp2d {

void loadArguments(int argc, char *argv[],
                   bool &checkInput, std::string &inputFileName,
                   bool &checkRestart, std::string &restartFileName,
                   bool &checkOutput, std::string &outputDirectory) {

// loadArguments will receive the arguments that were passed to the program from the command line as input.
// they are organized as following:
// .. argc = argument count, the number of arguments;
// .. argv = a set of characters which are the input arguments.
//
// To determine which argument is which, we will use the switches -i, -o, -r

/// GENERAL SCHEME OF THE SCRIPT //////////////////////////////////////
// First load arguments
// Check if input file exists, output directory exists or restart file exists (see stat() from POSIX).
// If not stop the program.
// Check if restart and input are given together or they are both missing.
// If it is so, stop the program.
// Check if directory exists.
// If not, create directory and check status of creation. On failure, stop the program.
// If it exists, clean just the output files from directory. Check status of cleaning.
// TODO: decide which extension can be given to the output files (.hf2o ??)
// Check possibility of writing in folder.
///////////////////////////////////////////////////////////////////////


// Now let us recover the arguments
// First, we check that we have at least 3 arguments of which the first is
// the name of the program (by default of C++), whereas the second and the
// third are (at least) a switch (input or reset) and the passed string.
// If no output directory is not given, the directory of the input file is used.

// The script will return which of the input arguments where passed and their value back to the main.

  // check the number of passed arguments
  if (argc < 3) {

    // If we do not have enough arguments, let us throw an error and close the program
    std::cout << "Usage of the program is " << std::endl;
    std::cout << "  HFPx2D -i inputfile -o outputdirectory     for new analyses;" << std::endl;
    std::cout << "  HFPx2D -r restartfile -o outputdirectory   for restarted analyses." << std::endl;
    std::cout << "-- Press ENTER to exit...";
    std::cin.get();
    exit(1);

  } else {

    // Initialize status variable to evaluate the existence of files
    int statusOfAccess;

    // if we got enough parameters, load them
    for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
                                      * Note that we're starting on 1 because we don't need to know the
                                      * path of the program, which is stored in argv[0] */

      if (std::string(argv[i]) == "-i") {

        // We know the next argument *should* be the input filename:
        inputFileName = argv[i + 1];

        // We add (to be sure) the initial "./" to the address
        inputFileName = "./" + inputFileName;

        // Check for existence and readability of the file
        statusOfAccess = access(inputFileName.c_str(), F_OK | R_OK);

        if (statusOfAccess != 0) {
          std::cerr << "Error: impossible to access the input file " << inputFileName << std::endl;
          std::cerr << strerror(errno) << std::endl;
          exit(1);
        } else {
          checkInput = true;
          i++;
        }

      } else if (std::string(argv[i]) == "-o") {

        // We know the next argument *should* be the output directory:
        outputDirectory = argv[i + 1];

        // We will check for the properties of the output folder later
        checkOutput = true;
        i++;

      } else if (std::string(argv[i]) == "-r") {

        // We know the next argument *should* be the restart filename:
        restartFileName = argv[i + 1];

        // We add (to be sure) the initial "./" to the address
        // inputFileName = "./" + inputFileName;
        inputFileName.insert(0,"./");

        // Check for existence and readability of the file
        statusOfAccess = access(restartFileName.c_str(), F_OK | R_OK);

        if (statusOfAccess != 0) {

          std::cerr << "Error: impossible to access the restart file " << inputFileName << std::endl;
          std::cerr << strerror(errno) << std::endl;
          exit(1);

        } else {

          // All good with restart
          checkRestart = true;
          i++;

        }

      } else {

        // There is an invalid argument, so send an error.
        std::cerr << "Invalid argument " << std::string(argv[i]) << std::endl;
        std::cerr << "-- Press ENTER to exit...";
        std::cin.get();
        exit(1);

      }

    }

    // now we check if both input and restart are inserted as argument
    if (!checkInput && !checkRestart) {

      std::cout << "No input or restart file set. Please provide an input file." << std::endl;
      std::cout << "-- Press ENTER to exit...";
      std::cin.get();
      exit(1);

    }

    // we are checking if checkInput and checkRestart are true, since we already checked that one of them is true
    if (checkInput && checkRestart) {

      std::cout << "Don't be schizophrenic... do you want to do a new analysis or restart one?" << std::endl;
      std::cout << "-- Press ENTER to exit...";
      std::cin.get();
      exit(1);

    }

    // check if the user gave an output folder or not, otherwise set the current folder as output
    if (!checkOutput) {

      std::cout << "No output directory set. Present directory will be taken as output directory." << std::endl;
      std::cout << "-- Press ENTER to continue...";
      std::cin.get();

    }

    /// Creating the OUTPUT directory
    // Construct the directory path
    outputDirectory = "./" + outputDirectory;

    if (checkOutput) { // if the output directory parameter has been passed

      // CASE 1: Directory is not there and we need to create one as required.
      // Try to create a directory
      int outputDirStatus = mkdir(outputDirectory.c_str(), 0777);

      // Check if it was possible to create the directory
      if (outputDirStatus != 0) {
        // If there is an error in creating the directory, then:
        // CASE 2: the directory is already there and we would like to clean it.
        // Check for existence, readability and possibility to write in the folder
        statusOfAccess = access(outputDirectory.c_str(), F_OK | W_OK | R_OK);

        if (statusOfAccess != 0) { // if there is an error on access as well, no option left:
                                   // we have to stop because there is no way to work on the folder
          std::cerr << "Error: impossible to create or access the directory " << outputDirectory << std::endl;
          std::cerr << strerror(errno) << std::endl;
          exit(1);

        } else {
          if (!checkRestart) {
            // The directory is accessible and we can work in it.
            // Let us clean the directory ONLY IF THE ANALYSIS IS NEW !!!
            cleanOutputDir(outputDirectory.c_str());
          }
        }
      }  // End "IF" for check directory status
    } // End "IF" for creation of directory
  } // End "IF" passed arguments number

} // End of script



///// This script removes the files in *path that have only some particular extensions
void cleanOutputDir(const char *path) {

  // open the stream of files in the directory *path
  DIR *dirFile = opendir(path);

  // if we have a successfull stream
  if (dirFile != nullptr) {

    // create the pointer to the list of files
    struct dirent *hFile;

    // for each file in the folder
    while ((hFile = readdir(dirFile)) != nullptr) {

      // Avoid dealing with current and parent folders
      if (!(std::string(hFile->d_name).compare(".")))
        continue;
      if (!(std::string(hFile->d_name).compare( "..")))
        continue;

      // in linux hidden files all start with '.' --> Ignore them or use a boolean variable such as gIgnoreHidden
      // which has to be set from the outside.
      //if ( gIgnoreHidden && ( hFile->d_name[0] == '.' )) continue;
      if (hFile->d_name[0] == '.')
        continue;

      // now we compare the extensions
      // save the file as string
      std::string fileName = hFile->d_name;
      // find location of the dot for the extension
      size_t dot = fileName.find_last_of('.'); // '.' is a char, "." is a string (that is a sequence of char, more expensive)

      // if the dot is found before the extension
      if (dot != std::string::npos) // npos is returned by find_last_of if the dot is not found
      {
        // compare the file extension with the ones to be removed
        std::string ext = fileName.substr(dot, fileName.size() - dot);
        if ((ext == "hf2o") || (ext == "txt")) { // TODO: ADD CORRECT FILE EXTENSIONS

          // Recover the full path for the file to delete
          const std::string full_file_name = std::string(path) + "/" + fileName;

          // Delete the file
          int removeStatus = remove(full_file_name.c_str());

          // Check successful delete
          if (removeStatus != 0) {
            std::cerr << "Error while trying to delete " << full_file_name << std::endl;
            std::cerr << strerror(errno) << std::endl;
            exit(1);
          }

        }

      } // else it is a folder or a file without extension and we do not touch it

    }
    // Close the stream
    closedir(dirFile);
  } else {
    std::cerr << "Error while trying to open directory " << path << " for cleaning." << std::endl;
    std::cerr << strerror(errno) << std::endl;
    exit(1);
  }
} // End of script

}
