/**
 * General DREAM C++ interface.
 */


#include <cmath>
#include <iostream>
#include <H5Cpp.h>
#include <string>
#include <unistd.h>

// If "not debugging" is defined, then we're in
// debug mode and would like to active floating-point
// exceptions
#ifndef NDEBUG
#   include <fenv.h>
#endif

#include <softlib/SOFTLibException.h>

#include "DREAM/config.h"
#include "DREAM/Init.h"
#include "DREAM/IO.hpp"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SFile.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Simulation.hpp"
#include "FVM/FVMException.hpp"


using namespace std;

struct cmd_args {
    bool display_settings=false;
    bool splash=true;
    string
        input_filename,
        output_filename;
};

void display_settings(DREAM::Settings *s=nullptr) {
    if (s == nullptr)
        s = DREAM::SimulationGenerator::CreateSettings();

    cout << endl << "LIST OF DREAM SETTINGS" << endl;
    cout <<         "----------------------" << endl;
    s->DisplaySettings();
}

/**
 * Print the DREAMi command-line argument help.
 */
void print_help() {
    cout << "Syntax: dreami INPUT [OPTIONS...]" << endl;
    cout << "   Run a thermal quench simulation according to the specifications made" << endl;
    cout << "   in the file 'INPUT'." << endl << endl;

    cout << "OPTIONS" << endl;
    cout << "  -h           Print this help." << endl;
    cout << "  -l           List all available settings in DREAM." << endl;
    cout << "  -o           Specify the name of the output file." << endl;
    cout << "  -s           Do not show the splash screen." << endl;
}

/**
 * Parse command-line arguments.
 *
 * argc: Number of command-line arguments (including program name).
 * argv: List of command-line arguments.
 */
struct cmd_args *parse_args(int argc, char *argv[]) {
    char c;

    struct cmd_args *a = new struct cmd_args;
    a->output_filename = "dream_output.h5";
    a->display_settings = false;

    while ((c = getopt(argc, argv, "hlo:s")) != -1) {
        switch (c) {
            case 'h':
                print_help();
                break;
            case 'l':
                display_settings();
                break;
            case 'o':
                a->output_filename = string(optarg);
                break;
            case 's':
                a->splash = false;
                break;
            case '?':
                if (optopt == 'o') {
                    cout << "Option -o requires an argument." << endl;
                    return nullptr;
                } else {
                    cout << "Unrecognized option: " << optopt << endl;
                    return nullptr;
                }
        }
    }

    if (optind+1 < argc) {
        cout << "Too many trailing input arguments." << endl;
        return nullptr;
    } else if (optind == argc) {
        cout << "No input file specified." << endl;
        return nullptr;
    }

    a->input_filename = string(argv[optind]);
    return a;
}

void splash() {
    cout << endl;
    cout << R"( It's time to...)" << endl;
    cout << endl;
    cout << R"( * * ________   _____  ______ ___* *  ___  ___       )" << endl;
    cout << R"(  * \\   __   \/ __  \/  ____//   \ //   \/   \      )" << endl;
    cout << R"(    //  / //  / /_// /  /__ // /\  \/  / / /  /   *  )" << endl;
    cout << R"(   //  / //  /      /   __/// /_/  /  / / /  /   * * )" << endl;
    cout << R"( _//  /_//  /  /\  \   /__//  __  /  / / /  /     *  )" << endl;
    cout << R"( \\________/__/ \\__\_____//_///_/__//_//__/     *   )" << endl;
    cout << R"(               * *          * *      *         *     )" << endl;
    cout << R"(      * *       *       * *      * *             *   )" << endl;
    cout << R"(     * * *                                 ...baby...)" << endl;
    cout << endl;
}


/**
 * Program entry point.
 *
 * argc: Number of command-line arguments (including program name).
 * argv: List of command-line arguments.
 */
int main(int argc, char *argv[]) {
    // Initialize the DREAM library
    //dream_initialize(&argc, &argv);
    dream_initialize();

    // Parse command-line arguments
    struct cmd_args *a = parse_args(argc, argv);
    if (a == nullptr)
        return -1;

    // Make sure that anything written to stdio/stderr is
    // written immediately
    cout.setf(ios_base::unitbuf);

    if (a->splash)
        splash();

    cout << "alpha version (commit " << DREAM_GIT_SHA1 << ")" << endl;

    // Except on NaN (but only in debug mode)
#ifndef NDEBUG
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif

    try {
        DREAM::Settings *settings = DREAM::SimulationGenerator::CreateSettings();
        
        DREAM::SettingsSFile::LoadSettings(settings, a->input_filename);

        DREAM::Simulation *sim = DREAM::SimulationGenerator::ProcessSettings(settings);
        sim->Run();
        sim->Save("output.h5");

        // TODO Generate output
        
    } catch (DREAM::FVM::FVMException &ex) {
        DREAM::IO::PrintError(ex.what());
    } catch (SOFTLibException &ex) {
        DREAM::IO::PrintError(ex.what());
    } catch (H5::FileIException &ex) {
        DREAM::IO::PrintError(ex.getDetailMsg());
    }

    // De-initialize the DREAM library
    dream_finalize();

    return 0;
}

