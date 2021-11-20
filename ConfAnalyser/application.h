#ifndef APPLICATION_H
#define APPLICATION_H

#include "molecule.h"
#include <vector>
#include <map>
#include <string>

class Application
{
        public:
                Application(int _argc = 0, char **_argv = nullptr);
                ~Application();
                int run();
        private:
                bool read_PDB(std::string &file_name,
                              std::vector<Atom*> &molecule);
                bool read_atom_names();
                bool process_file(Molecule* mol, std::string file_name);
                void help() const;
                void parse_options();
                void results(std::vector<Molecule*> molecules);
                std::vector<Molecule*> molecules;
                int argc;
                char ** argv;
                bool print_summary;
                bool print_list;
                int analysis_type;
                std::string input_file_list;
                std::string atom_names_list;
};

#endif
