/**
 * General TQS C++ interface.
 */


#include <iostream>
#include <H5Cpp.h>
#include <string>
#include <unistd.h>

#include <softlib/SOFTLibException.h>

#include "TQS/config.h"
#include "TQS/Init.h"


using namespace std;

struct cmd_args {
    string
        input_filename,
        output_filename;
};

/**
 * Print the TQSi command-line argument help.
 */
void print_help() {
    cout << "Syntax: tqsi INPUT [OPTIONS...]" << endl;
    cout << "   Run a thermal quench simulation according to the specifications made" << endl;
    cout << "   in the file 'INPUT'." << endl << endl;

    cout << "OPTIONS" << endl;
    cout << "  -h           Print this help." << endl;
    cout << "  -o           Specify the name of the output file." << endl;
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
    a->output_filename = "tqs_output.h5";

    while ((c = getopt(argc, argv, "ho:")) != -1) {
        switch (c) {
            case 'h':
                print_help();
                break;
            case 'o':
                a->output_filename = string(optarg);
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


/**
 * Program entry point.
 *
 * argc: Number of command-line arguments (including program name).
 * argv: List of command-line arguments.
 */
int main(int argc, char *argv[]) {
    // Initialize the TQS library
    tqs_initialize(&argc, &argv);

    // Parse command-line arguments
    struct cmd_args *a = parse_args(argc, argv);
    if (a == nullptr)
        return -1;

    // Make sure that anything written to stdio/stderr is
    // written immediately
    cout.setf(ios_base::unitbuf);

    cout << "Welcome to TQS alpha (commit " << TQS_GIT_SHA1 << ")" << endl;

    try {
        // DO WORK
    } catch (SOFTLibException &ex) {
        cout << ex.what() << endl;
    } catch (H5::FileIException &ex) {
        cout << ex.getDetailMsg() << endl;
    }

    // De-initialize the TQS library
    tqs_finalize();

    return 0;
}

