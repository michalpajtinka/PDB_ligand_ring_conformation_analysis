/*
 *	file: atom.h
 */

#ifndef ATOM_H
#define ATOM_H

#include "point_3D.h"
#include <iostream>
#include <string>

/* Class representing single atom from PDB structure */
class Atom : public Point_3D
{
        public:
                /* values for setting recordType */
                static const int RECORD_UNKNOWN = 0;
                static const int RECORD_ATOM = 1;
                static const int RECORD_HETATM = 2;
 
                /* constructor */
                Atom();

                /* reading one line of PDB file */
                void read_entry(const std::string &line);

                /* writing one line to PDB file */
                void write_entry(std::ofstream &ofile);

                /* print testing info about atom */
                void print() const;

                /* set line number */
                void set_line_number(size_t num);

                /* get atom name */
                std::string get_atom_name() const;

		/* get residue name */
                std::string get_residue_name() const;

		/* get residue number */
		int get_residue_number() const;

		/* get element name */
                std::string get_element_name() const;
        private:
                /* atom attributes */
                size_t line_number; /* keep line number in case of error */
                int record_type;
                int atom_number;
                std::string atom_name;
                char alternate_location;
                std::string residue_name;
                char chain_id;
                int residue_number;
                char i_code;
                double occupancy;
                double temp_factor;
                std::string element_name;
                std::string formal_charge;
                bool is_occupancy;
                bool is_temp_factor;
};


#endif
