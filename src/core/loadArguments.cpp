//
// Created by lorenzo on 9/2/17.
//

#include <iostream>
#include <fstream>
#include <sys/stat.h>
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

    // if we got enough parameters, load them
    for (int i = 1; i < argc; i++) { /* We will iterate over argv[] to get the parameters stored inside.
                                      * Note that we're starting on 1 because we don't need to know the
                                      * path of the program, which is stored in argv[0] */

      if (std::string(argv[i]) == "-i") {

        // We know the next argument *should* be the input filename:
        inputFileName = argv[i + 1];
        checkInput = true;
        i++;
        std::cout << " Loaded input " << std::endl;

      } else if (std::string(argv[i]) == "-o") {

        // We know the next argument *should* be the output directory:
        outputDirectory = argv[i + 1];
        checkOutput = true;
        i++;
        std::cout << " Loaded output " << std::endl;

      } else if (std::string(argv[i]) == "-r") {

        // We know the next argument *should* be the restart filename:
        restartFileName = argv[i + 1];
        checkRestart = true;
        i++;

      } else {

        // There is an invalid argument, so send an error.
        std::cout << "Invalid argument " << std::string(argv[i]) << std::endl;
        std::cout << "-- Press ENTER to exit...";
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

    if (checkOutput) { // if the output directory parameter is not empty

      // Try to create a directory
      int outputDirStatus = mkdir(outputDirectory.c_str(), 0777);

      // Check if it was possible to create the directory
      if (outputDirStatus == -1) { // If there is an error in creating the directory

        // If the directory is already there, then try to remove the directory and its content
        system(("exec rm -r " + outputDirectory).c_str());

        // Try to create new empty folder
        int outputDirStatus = mkdir(outputDirectory.c_str(), 0777);

        // If also this time it is not possible, stop the program
        if (outputDirStatus != 0) {
          std::cerr << "Error " << outputDirStatus
                    << " in creating the output directory " << outputDirectory
                    << std::endl;
          std::cout << "-- Press ENTER to exit...";
          std::cin.get();
          exit(1);
        }
      }
    } // End "IF" directory creation
  } // End "IF" passed arguments number

} // End of script

}