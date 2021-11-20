#include "molecule.h"
#include <iomanip>

using namespace std;

map<string, short> Molecule::conformations = {
        {"UNANALYSED", 0},
        {"UNDEFINIED", 1}
};


Molecule::Molecule(string _structure)
{
        structure = _structure;
	ligand = "";
        conformation = conformations["UNANALYSED"];
        filled = false;
        analysed = false;
}


Molecule::~Molecule() {}


short Molecule::get_conformation() const
{
        return conformation;
}


string Molecule::translate_conformation() const
{
        for (auto conf : conformations) {
                if(conf.second == conformation) {
                        return conf.first;
                }
        }

        return "";
}


ostream& Molecule::print(ostream& out)
{
        size_t sep = structure.find_last_of("/");
        string tmp = (sep == string::npos) ? structure :
                structure.substr(sep + 1, structure.size() - sep - 1);
        return out << tmp << ": " << translate_conformation() << endl;
}


void Molecule::statistics(const std::vector<Molecule*> vec)
{
        size_t sum = 0;
        size_t *conf_num = new size_t[conformations.size()]();

	for (auto x : vec) {
	        sum++;
                conf_num[x->get_conformation()]++;
        }

        for (auto conf : conformations) {
                cout << setw(14) << left << string(conf.first)+": "
                     << conf_num[conf.second]
                     << " ("
                     << conf_num[conf.second] / (float)sum * 100
                     << "%)"
                     << endl;
        }
        cout << setw(14) << left << "TOTAL: " << sum << endl;

        delete[] conf_num;
}


ostream& operator<<(ostream& out, Molecule &mol)
{
        return mol.print(out);
}
