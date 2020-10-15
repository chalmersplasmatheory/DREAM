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
#   include <csignal>
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
    bool print_adas=false;
    bool splash=true;
    bool save_initial = false;
    string
        input_filename,
        initial_filename,
        output_filename;
};

void display_settings(DREAM::Settings *s=nullptr) {
    if (s == nullptr)
        s = DREAM::SimulationGenerator::CreateSettings();

    cout << endl << "LIST OF DREAM SETTINGS" << endl;
    cout <<         "----------------------" << endl;
    s->DisplaySettings();
}

void display_adas(DREAM::Simulation *sim) {
    sim->GetADAS()->PrintElements();
    cout << endl;
}

/**
 * Print the DREAMi command-line argument help.
 */
void print_help() {
    cout << "Syntax: dreami INPUT [OPTIONS...]" << endl;
    cout << "   Run a thermal quench simulation according to the specifications made" << endl;
    cout << "   in the file 'INPUT'." << endl << endl;

    cout << "OPTIONS" << endl;
    cout << "  -a           Print list of elements in ADAS database." << endl;
    cout << "  -h           Print this help." << endl;
    cout << "  -i           Save the initial simulation state to the named file." << endl;
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
    a->initial_filename = "";
    a->output_filename = "";
    a->display_settings = false;

    while ((c = getopt(argc, argv, "ahi:lo:s")) != -1) {
        switch (c) {
            case 'a':
                a->print_adas = true;
                break;
            case 'h':
                print_help();
                break;
            case 'i':
                a->save_initial = true;
                a->initial_filename = string(optarg);
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
 * Handle floating point exceptions.
 */
void sig_fpe(int) {
    throw DREAM::FVM::FVMException("Floating-point error.");
}

/**
 * Construct fake command-line arguments.
 */
char ***construct_fake_args(vector<string> &args, int &argc) {
    argc = args.size();

    char ***argv = new char**;
    // +1: PETSc is apparently buggy...
    *argv = new char*[argc+1];

    for (int i = 0; i < argc; i++) {
        size_t l = args[i].size();
        (*argv)[i] = new char[l+1];
        args[i].copy((*argv)[i], l);
        (*argv)[i][l] = 0;
    }

    (*argv)[argc] = nullptr;

    return argv;
}


/**
 * Program entry point.
 *
 * argc: Number of command-line arguments (including program name).
 * argv: List of command-line arguments.
 */
int main(int argc, char *argv[]) {
    int exit_code = 0;

    // The code below can be used to make PETSc print a list of
    // citations to cite based on the current simulation
    /*int argc2;
    vector<string> args({"dreami", "-citations", "petsc-citations.txt"});
    char ***argv2 = construct_fake_args(args, argc2);
    dream_initialize(&argc2, argv2);*/

    // Initialize the DREAM library
    dream_initialize();
#ifndef NDEBUG
    PetscPopSignalHandler();
    std::signal(SIGFPE, sig_fpe);
#endif

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

    DREAM::Simulation *sim = nullptr;
    try {
        DREAM::Settings *settings = DREAM::SimulationGenerator::CreateSettings();
        
        DREAM::SettingsSFile::LoadSettings(settings, a->input_filename);

        sim = DREAM::SimulationGenerator::ProcessSettings(settings);

        if (a->print_adas)
            display_adas(sim);

        if (a->save_initial)
            sim->Save(a->initial_filename);

        sim->Run();
    } catch (DREAM::FVM::FVMException &ex) {
        DREAM::IO::PrintError(ex.what());
        exit_code = 1;
    } catch (SOFTLibException &ex) {
        DREAM::IO::PrintError(ex.what());
        exit_code = 2;
    } catch (H5::FileIException &ex) {
        DREAM::IO::PrintError(ex.getDetailMsg().c_str());
        exit_code = 3;
    }

    if (sim != nullptr) {
        try {
            if (a->output_filename != "")
                sim->Save(a->output_filename);
            else
                sim->Save();
        } catch (H5::FileIException &ex) {
            DREAM::IO::PrintError(ex.getDetailMsg().c_str());
            exit_code = 4;
        }
    }

    // De-initialize the DREAM library
    dream_finalize();

    return exit_code;
}

