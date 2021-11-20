/*
 *	file: atom.cpp
 */

#include "atom.h"
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

Atom::Atom()
{
        record_type =  RECORD_UNKNOWN;
        atom_number = 0;
        atom_name = "    ";
        alternate_location = ' ';
        residue_name = "   ";
        chain_id = ' ';
        residue_number = 0;
        i_code = ' ';
        X = 0;
        Y = 0;
        Z = 0;
        occupancy = 0;
        temp_factor = 0;
        element_name = "  ";
        formal_charge = "  ";
        is_occupancy = false;
        is_temp_factor = false;
}

void Atom::read_entry(const string &line)
{
        string record_name;
        string s; /* string temporary storing parsed out information */
        istringstream sstream;
        /* check minimal line size */
        if (line.length() < 53) {
                cout<<"The line number " << line_number << " is too short!" << endl;
                return;
        }

        /* record type */
        record_name = line.substr(0, 6);
        if (record_name == "ATOM  ") {
                record_type = RECORD_ATOM;
        } else if (record_name == "HETATM") {
                record_type = RECORD_HETATM;
        } else {
                return;
        }

        /* atom number */
        s = line.substr(6, 5);
        sstream.str(s);
        sstream.clear();
        sstream >> atom_number;

        if (sstream.fail()) {
                cout << "Line " << line_number
                     << ": Error reading atom number!" << endl;
                return;
        }

        /* atom name */
        atom_name = line.substr(12, 4);

        /* location indicator */
        alternate_location = line[16];

        /* residue name */
        residue_name = line.substr(17, 3);

        /* chain ID */
        chain_id = line[21];

        /* residue number */
        s = line.substr(22, 4);
        sstream.str(s);
        sstream.clear();
        sstream >> residue_number;

        /* i code */
        i_code = line[26];

        /* reading X coordinate */
        s = line.substr(30, 8);
        sstream.str(s);
        sstream.clear();
        sstream >> X;
        if (sstream.fail()) {
                cout << "Line " << line_number
                     << ": Error while reading coordinates!" << endl;
                return;
        }

        /* reading Y coordinate */
        s = line.substr(38, 8);
        sstream.str(s);
        sstream.clear();
        sstream >> Y;
        if (sstream.fail()) {
                cout << "Line " << line_number
                     << ": Error while reading coordinates!" << endl;
                return;
        }

        /* reading Z coordinate */
        s = line.substr(46, 8);
        sstream.str(s);
        sstream.clear();
        sstream >> Z;
        if (sstream.fail()) {
                cout << "Line " << line_number
                     << ": Error while reading coordinates!" << endl;
                return;
        }

        /* read occupancy, if the line is long enought */
        if (line.length() >= 60) {
                s = line.substr(54, 6);
                sstream.str(s);
                sstream.clear();
                sstream >> occupancy;
                if (!sstream.fail()) {
                        is_occupancy = true;
                }
        }

        /* tempFactor */
        if (line.length() >= 66) {
                s = line.substr(60, 6);
                sstream.str(s);
                sstream.clear();
                sstream >> temp_factor;
                if (!sstream.fail()) {
                        is_temp_factor = true;
                }
        }

        /* elementName */
        if (line.length() >= 78) {
                element_name = line.substr(76, 2);
        }

        /* formalCharge */
        if (line.length() >= 80) {
                formal_charge = line.substr(78, 2);
        }
}


void Atom::write_entry(ofstream &ofile)
{
        /* print record name */
        if (record_type == RECORD_ATOM) {
                ofile << "ATOM  ";
        } else if (record_type == RECORD_HETATM) {
                ofile << "HETATM";
        } else {
                return;
        }

        /* atom number, 5 chars aligned to right */
        ofile << right << setw(5) << atom_number;

        /* obligatory white space */
        ofile << ' ';

        /* atom name, 4 chars aligned to left */
        ofile << left << setw(4) << atom_name;

        /*alternate location */
        ofile << alternate_location;

        /* residue name, 3 chars aligned to left */
        ofile << left << setw(3) << residue_name;

        /* obligatory space */
        ofile << ' ';

        /* chain ID */
        ofile << chain_id;

        /* residue number, 4 chars aligned to right */
        ofile << right << setw(4) << residue_number;

        /* insert code, 1 char */
        ofile << i_code;

        /* 3 obligatory spaces */
        ofile << "   ";

        /* Cartesian coordinates, manipulator fixed, 8 chars, 3 decimals */
        ofile << right << fixed << setprecision(3);
        ofile << setw(8) << X;
        ofile << setw(8) << Y;
        ofile << setw(8) << Z;

        /* occupancy or spaces */
        if (is_occupancy) {
                ofile << right << fixed << setprecision(2)
                                                << setw(6)<< occupancy;
        } else {
                ofile << "      ";
        }

        /* temp factor or spaces */
        if (is_temp_factor) {
                ofile << right << fixed << setprecision(2)
                                                << setw(6) << temp_factor;
        } else {
                ofile << "      ";
        }

        /* 10 obligator spaces */
        ofile << "          ";

        /* element name */
        ofile << right << setw(2) << element_name;

        /* formal charge */
        ofile << left << setw(2) <<formal_charge;

        /* new line character */
        ofile << endl;
}


void Atom::print() const
{
        cout << atom_number << " " << atom_name << " [" << X << ", " <<
                        Y << ", " << Z << "] " <<  element_name << endl;
}


void Atom::set_line_number(size_t num)
{
        line_number = num;
}


string Atom::get_atom_name() const
{
        return atom_name;
}


string Atom::get_residue_name() const
{
	return residue_name;
}


int Atom::get_residue_number() const
{
	return residue_number;
}
