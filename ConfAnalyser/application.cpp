#include "application.h"
#include "plane_3D.h"
#include "cyclohexane.h"
#include "cyclopentane.h"
#include "benzene.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <getopt.h>

#define EMPTY        -1
#define CYCLOHEXANE   0
#define CYCLOPENTANE  1
#define BENZENE       2

using namespace std;


map<string, vector<vector<string>>> Molecule::atom_names;


Application::Application(int _argc, char **_argv)
{
        argc = _argc;
        argv = _argv;
        print_summary = true;
        print_list = true;
        analysis_type = EMPTY;
        string input_file_list = string();
}


Application::~Application()
{
        for (vector<Molecule*>::iterator itr = molecules.begin();
                                        itr != molecules.end(); ++itr) {
                delete(*itr);
        }
        molecules.clear();
}


bool Application::read_PDB(string &file_name, vector<Atom*> &molecule)
{
        string line;
        string record_name;

        ifstream ifile;
        ifile.open(file_name);
	if (ifile.fail()) {
                cerr << "Could not open file " << file_name << "..." << endl;
                return false; 
        }

	size_t line_number = 1; /* keep line number for case of error */
       	while (!ifile.eof()) {
	        getline(ifile, line);
                if (line.length() >= 6) {
                        record_name = line.substr(0, 6); 
                        if (record_name == "ATOM  " || record_name == "HETATM") {
                                Atom *atom = new Atom();
                                atom->set_line_number(line_number);
                                atom->read_entry(line);
                                molecule.push_back(atom);
                        }
                }
                line_number++;
        }

        ifile.close();

        return true;
}


bool Application::read_atom_names()
{
        size_t atoms_count = 0;
        switch (analysis_type) {
                case CYCLOPENTANE:
                        atoms_count = 5;
                        break;
                case CYCLOHEXANE:
                        atoms_count = 6;
                        break;
                case BENZENE:
                        atoms_count = 6;
                        break;
                default:
                        cerr << "Can`t deduce number of atoms from given analysis type!" << endl;
                        return false;
        }

        ifstream ifile;
        ifile.open(atom_names_list);
	if (ifile.fail()) {
                cerr << "Could not open file " << atom_names_list << "..." << endl;
                return false; 
        }

        string line;
	size_t line_number = 1; /* keep line number for case of error */
        stringstream ss;
       	while (getline(ifile, line)) {
	        ss.clear();
	        ss.str(line);

                /* Get ligand name */
                string ligand_name;
                ss >> ligand_name;

                if (ss.fail()) {
                        ss.clear();
                        cerr << atom_names_list << ": Wrong syntax on line nr. "
                             << line_number << "..." << endl;
                        line_number++;
                        continue;
                }

                /* Read atom names to temporary container */
                vector<string> tmp_vec;
                string tmp_str;
                while (ss >> tmp_str) {
                        ss.clear();
                        tmp_vec.push_back(tmp_str);
                }

                if (tmp_vec.size() != atoms_count) {
                        cerr << atom_names_list << ": Wrong number of atom names on line nr. "
                             << line_number << " (expected " << atoms_count << ", was "
                             << tmp_vec.size() << "), entry ommited..." << endl;
                        line_number++;
                        continue;
                }

                /* Create new lingand entry, if this one doesn`t exist */
                if (Molecule::atom_names.count(ligand_name) == 0) {
                        vector<vector<string>> all_names_list;
                        Molecule::atom_names.insert(map<string, vector<vector<string>>>::value_type(ligand_name, all_names_list));
                }

                /* fill the atom names */
                for (size_t i = 0; i < atoms_count; i++) {
                        if (Molecule::atom_names.at(ligand_name).size() <= i) {
                                vector<string> current_names_list;
                                Molecule::atom_names.at(ligand_name).push_back(current_names_list);
                        }
                        Molecule::atom_names.at(ligand_name).at(i).push_back(tmp_vec[i]);
                }
                
                line_number++;
        }

        ifile.close();

        return true;

}

bool Application::process_file(Molecule* mol, string file_name)
{
        vector<Atom*> molecule;
        
        if (!read_PDB(file_name, molecule))
        {
                cerr << "Problem while reading PDB " << file_name << "...\n";
                goto ERROR;
        }

        if (!mol->initialize(molecule))
        {
                cerr << "Problem while initializing molecule " << file_name << "...\n";
                goto ERROR;
        }

        if (!mol->analyse()) {
                cerr << "Problem while analysing " << file_name << "...\n";
                goto ERROR;
        }

        return true;
 ERROR:
        cout << file_name << ": ommited\n";
        return false;
 }


void Application::help() const
{
        cout << "Usage:" << endl;
        cout << "   " << argv[0]
             << " [-h] -i file_list.txt -n name_list.txt --(ring_type) [-l | -s | -a]"
             << endl << endl ;
        cout << "Required:" << endl;
        cout << "   -i --input_list=FILE" << endl
             << "      read list of molecules to process from FILE - each line is treated as path to single PDB file" << endl;
        cout << "   -n --name_list=FILE" << endl
             << "      read list of names of atom ring from FILE. Each line represents one ligand, first word on the" << endl
             << "      line is treated as ligand name, all the following words are treated as atom names (if ligand is" << endl
             << "      not known or if name of the atom is not found in this list, processed atom will be ommited)." << endl
             << "      In case of multiple name variations, more lines with the same ligand name has to be present." << endl
             << "      Atom order matters!" << endl;
        cout << "   --(ring_type)" << endl
             << "      perform analysis of this type of molecule ring" << endl
             << "      currently supported:" << endl
             << "         --cyclohexane" << endl
             << "         --cyclopentane" << endl
             << "         --benzene" << endl;
        cout << "Optional:" << endl;
        cout << "   -h --help" << endl
             << "      display this help" << endl;
        cout << "   -l --list" << endl
             << "      display results only as a list of molecules and their conformations" << endl;
        cout << "   -s --summary" << endl
             << "      display results only as a short summary of relative occurances of conformations among tested molecules" << endl;
        cout << "   -a --all" << endl
             << "      display both list and summary (turned on by default, unless one of -l/-s options is detected)" << endl;
}


void Application::parse_options()
{
        /* option returned by getopt_long */
        int opt = 0;
        /* getopt_long stores the option index here. */
        int index = 0;
        /* watch that only one of s/l/a options will by set */
        bool display_option_set = false;
        bool analysis_type_set = false;
        /* long options */
        static struct option long_opt[] =
        {
                {"help",         no_argument,       nullptr,        'h'},
                {"cyclohexane",  no_argument,       &analysis_type, 0  },
                {"cyclopentane", no_argument,       &analysis_type, 1  },
                {"benzene",      no_argument,       &analysis_type, 2  },
                {"list",         no_argument,       nullptr,        'l'},
                {"summary",      no_argument,       nullptr,        's'},
                {"all",          no_argument,       nullptr,        'a'},
                {"input_list",   required_argument, nullptr,        'i'},
                {"name_list",    required_argument, nullptr,        'n'},
                {0, 0, 0, 0}
        };
        /* short options */
        static const char *short_opt = "hlsai:n:";

        /* Proces all of the arguments */
        while(true) {
                opt = getopt_long(argc, argv, short_opt, long_opt, &index);
                if (opt == EMPTY) {
                        break;
                }

                switch (opt) {
                        case 0:
                                if (*(long_opt[index].flag) != -1 && !analysis_type_set) {
                                        analysis_type_set = true;
                                        break;
                                }
                                cout << "Single and valid molecule type is required!";
                                goto END;
                        case 'h':
                                help();
                                exit(0);
                        case 'l':
                                if (display_option_set) {
                                        cout << "Only one of -l/-s/-a options can be specified!";
                                        goto END;
                                }
                                display_option_set = true;
                                print_summary = false;
                                break;
                        case 's':
                                if (display_option_set) {
                                        cout << "Only one of -l/-s/-a options can be specified!";
                                        goto END;
                                }
                                display_option_set = true;
                                print_list = false;
                                break;
                        case 'a':
                                if (display_option_set) {
                                        cout << "Only one of -l/-s/-a options can be specified!";
                                        goto END;
                                }
                                display_option_set = true;
                                break;
                        case 'i':
                                if (optarg == nullptr) {
                                        cout << "Input file was not specified!";
                                        goto END;
                                }
                                input_file_list = optarg;
                                break;
                        case 'n':
                                if (optarg == nullptr) {
                                        cout << "List of atom names was not provided!";
                                        goto END;
                                }
                                atom_names_list = optarg;
                                break;
                        case '?':
                                cout << "Terminating...";
                                goto END;
                        default:
                                cout << "Error: getopt_long() returned unexpected char: '" << opt << "'" << endl;
                                abort ();
                }
        }

        /* check that required arguments were found */
        if (input_file_list.empty() || analysis_type == EMPTY) {
                cout << "Some required arguments are missing!";
                goto END;
        }

        return;

END:
        cout << endl << endl;
        help();
        exit(0);
}


void Application::results(vector<Molecule*> molecules){
        if (print_list) {
                for (auto x : molecules) {
                        cout << *x;
                }
        }

        if (print_list && print_summary) {
                cout << endl;
        }

        if (print_summary) {
                if (!molecules.empty()) {
                        cout << "SUMMARY" << endl << "-------" << endl;
                        Molecule::statistics(molecules);
                } else {
                        cout << "No molecules detected!" << endl;
                }
        }
}


int Application::run()
{
        /* Parse command line arguments */
        parse_options();

        /* Read list of atom names */
        if (!read_atom_names()) {
                return EXIT_FAILURE;
        }

        /* Open list of molecules */
        string line;
        ifstream f(input_file_list);
        if (!f.is_open()) {
                cerr << "Error while opening file " << input_file_list << endl;
                return EXIT_FAILURE;
        }

        /* Read molecules from list of molecules and proccess them */ 
        Molecule* tmp;
        while(getline(f, line)) {
                switch (analysis_type) {
                        case CYCLOHEXANE:
                                tmp = new Cyclohexane(line);
                                break;
                        case CYCLOPENTANE:
                                tmp = new Cyclopentane(line);
                                break;
                        case BENZENE:
                                tmp = new Benzene(line);
                                break;
                        default:
                                cerr << "Unknown type of analysis!" << endl;
                                return EXIT_FAILURE;
                }

                if (!process_file(tmp, line)) {
                        delete(tmp);
                        continue;
                }

                molecules.push_back(tmp);
        }

        /* Print results */
        results(molecules);

        return EXIT_SUCCESS; 
}
