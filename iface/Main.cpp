/**
 * General DREAM C++ interface.
 */


#include <iostream>
#include <H5Cpp.h>
#include <string>
#include <unistd.h>

#include <softlib/SOFTLibException.h>

#include "DREAM/config.h"
#include "DREAM/Init.h"
#include "DREAM/Settings/Settings.hpp"
#include "DREAM/Settings/SFile.hpp"
#include "DREAM/Settings/SimulationGenerator.hpp"
#include "DREAM/Simulation.hpp"


using namespace std;

struct cmd_args {
    bool splash=true;
    string
        input_filename,
        output_filename;
};

/**
 * Print the DREAMi command-line argument help.
 */
void print_help() {
    cout << "Syntax: dreami INPUT [OPTIONS...]" << endl;
    cout << "   Run a thermal quench simulation according to the specifications made" << endl;
    cout << "   in the file 'INPUT'." << endl << endl;

    cout << "OPTIONS" << endl;
    cout << "  -h           Print this help." << endl;
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

    while ((c = getopt(argc, argv, "ho:s")) != -1) {
        switch (c) {
            case 'h':
                print_help();
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

void display_settings(DREAM::Settings *s) {
    cout << endl << "LIST OF DREAM SETTINGS" << endl;
    cout <<         "----------------------" << endl;
    s->DisplaySettings();
}

void splash() {
    cout << endl;
    cout << R"( It's time to...)" << endl;
    cout << endl;
    cout << R"( * *  _______   _____  ______  __* *  ___  ___       )" << endl;
    cout << R"(  *  \   __   \/ __  \/  ____//   \ */   \/   \      )" << endl;
    cout << R"(     /  /  /  / /__| /  /__  / /\  \/  / / /  /   *  )" << endl;
    cout << R"(    /  /  /  /      /   __/ / /_/  /  / / /  /   * * )" << endl;
    cout << R"(  _/  /__/  /  /\  \   /___/  __  /  / / /  /     *  )" << endl;
    cout << R"(  \________/__/  \__\_____//_/ /_/__/___/__/     *   )" << endl;
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
    dream_initialize(&argc, &argv);

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

    try {
        DREAM::Settings *settings = DREAM::SimulationGenerator::CreateSettings();
        display_settings(settings);
        /*DREAM::SettingsSFile::LoadSettings(settings, a->input_filename);

        DREAM::Simulation *sim = DREAM::SimulationGenerator::ProcessSettings(settings);
        sim->Run();*/

        // TODO Generate output
        
    } catch (SOFTLibException &ex) {
        cout << ex.what() << endl;
    } catch (H5::FileIException &ex) {
        cout << ex.getDetailMsg() << endl;
    }

    // De-initialize the DREAM library
    dream_finalize();

    return 0;
}

