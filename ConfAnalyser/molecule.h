#ifndef MOLECULE_H
#define MOLECULE_H

#include "atom.h"
#include <iostream>
#include <vector>
#include <map>
#include <string>

class Molecule
{
        public:
                Molecule() = delete;
                Molecule(std::string _structure);
                virtual ~Molecule();
                short get_conformation() const;
                virtual std::string translate_conformation() const;
                virtual std::ostream& print(std::ostream& out);
                virtual bool initialize(const std::vector<Atom*> &atoms) = 0;
                virtual bool analyse() = 0;
                static void statistics(const std::vector<Molecule*> vec);
                friend std::ostream& operator<<(std::ostream& out,
                                                        Molecule &mol);
                /* List of names of ring atoms in given ligand */
                static std::map<std::string,
                        std::vector<std::vector<std::string>>> atom_names;
        protected:
                /* Possible conformations, molecule-specific */
                static std::map<std::string, short> conformations;
                /* Data members */
                std::string structure;
		std::string ligand;
                short conformation;
                bool filled;
                bool analysed;
};


#endif
