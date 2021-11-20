#!/usr/bin/python
import sys
import os
import platform
import subprocess
import argparse
import re

# Object representing single conformation
class Conformation:
    total_num = 0

    # Constructor
    def __init__(self, _name, _template):
        self.name = _name.strip()
        self.template = _template
        self.num = 0

    # Calculating RMSD of given file when compared to the template
    # Expects executable SiteBinderCMD.exe stored in ./SBCMD/ directory
    # Expects pdb files to be examined stored in ./structures/ directory
    def get_RMSD(self, _program_executable, use_mono, _file):
        RMSD = sys.float_info.max
        for _template in self.template:
            mono = "wine" if use_mono else ""
            args_SiteBinder = [mono, _program_executable, _file, _template]
            try:
                process_SiteBinder = subprocess.Popen(args_SiteBinder,
                                            stdout=subprocess.PIPE, shell=False)
                output = process_SiteBinder.communicate()[0]
                process_SiteBinder.stdout.close()
            except BaseException as error:
                print "An exception occurred while trying to execute following command:"
                print "   "+" ".join(args_SiteBinder)
                print "Terminating program..."
                sys.exit()
            #Verify that SiteBinder returned the expected result
            pattern = re.compile("^RMSD = \d+[.]\d* Atoms = \d.*$")
            tmp = sys.float_info.max
            if pattern.match(output):
                try:
                    tmp = float(output.split()[2])
                except ValueError:
                        print "Data doesn`t match the expected result format:"
                        print "   " + output
            else:
                print "---------"
                print "WARNING: ANALYSIS MAY BE INACCURATE!"
                print "Command \""+" ".join(args_SiteBinder) + "\" has returned the following message:"
                print output
                print "---------" 
            if tmp < RMSD:
                RMSD = tmp
        return RMSD

    # Find relative occurance of this conformation among all conformations
    def get_percentage(self):
        if Conformation.total_num == 0:
            return 0
        else:
            return (float(self.num) / Conformation.total_num) * 100

class Molecule:
    # Constructor
    def __init__(self, _sitebinder, _name, _location, conformation_list, use_mono, tolerance):
        self.name = _name
        self.location = _location
        self.RMSD_list = []
        self.__find_RMSD_values(_sitebinder, use_mono, conformation_list)
        if tolerance is None:
            self.conformation = self.__determine_conformation(conformation_list)
        else:
            self.conformation = self.__determine_conformation(conformation_list, tolerance)
        self.__add_to_statistics()

    # Function to find RMSD with every conformation
    def __find_RMSD_values(self, _sitebinder, use_mono, conformation_list):
        for conf in conformation_list[:-1]:
            self.RMSD_list.append(conf.get_RMSD(_sitebinder, use_mono, self.get_path()))

    # Function to find conformation of given molecule
    def __determine_conformation(self, conformation_list, tolerance = sys.float_info.max):
        conformation = conformation_list[-1]
        lowest_RMSD = sys.float_info.max
        for i, conf in enumerate(conformation_list[:-1]):
            if self.RMSD_list[i] < min(lowest_RMSD, tolerance):
                lowest_RMSD = self.RMSD_list[i]
                conformation = conf
        return conformation

    # Function to add molecule to statistics
    def __add_to_statistics(self):
        self.conformation.num += 1;
        Conformation.total_num += 1;

    # Function to generate path to file
    def get_path(self):
        delim = '\\' if platform.system() == "Windows" else '/'
        return self.location + ("" if self.location[-1] == delim else delim) + self.name


# Function to generate CSV with all RMSDs of given molecule
def generate_RMSD_chart(conformation_list, molecules):
    sys.stdout.write(';')
    for conf in conformation_list[:-1]:
        sys.stdout.write(conf.name+";")
    sys.stdout.write('\n')
    for mol in molecules:
        sys.stdout.write(mol.name+";")
        for RMSD in mol.RMSD_list:
                sys.stdout.write(str(RMSD)+";")
        sys.stdout.write('\n')


# Function to generate list of files and their conformations
def generate_result_list(molecules):
    for mol in molecules:
        sys.stdout.write(mol.name+": "+mol.conformation.name+"\n")


# Function to generate overall summary
def generate_summary(conformations):
    for conf in conformations:
        sys.stdout.write("%-13s %3d (%.2f%%)\n" % (conf.name+":", conf.num,
                                                    conf.get_percentage()))
    sys.stdout.write("%-13s %3d\n" % ("TOTAL:", Conformation.total_num))


# Function to process every file in given directory
def process_files(_sitebinder, _dir, conf_list, use_mono, tolerance):
    molecules = []
    for _file in sorted(os.listdir(_dir)):
        if _file.endswith(".pdb"):
            tmp = Molecule(_sitebinder, _file, _dir, conf_list, use_mono, tolerance)
            molecules.append(tmp)
    return molecules


# Function to fill list of conformations:
def fill_conformation_list(template_location):
    conformation_list = []
    for _conf in os.listdir(template_location):
        file_list = []
        delim = '\\' if platform.system() == "Windows" else '/'
        for _file in os.listdir(template_location+delim+_conf):
            if _file.endswith(".pdb"):
                tmp = template_location + ("" if template_location[-1] == delim else delim) + _conf + delim + _file
                file_list.append(tmp)
        conformation_list.append(Conformation(_conf, file_list))
    conformation_list.append(Conformation("UNKNOWN", []))
    return conformation_list


def main():
    parser = argparse.ArgumentParser(description="Determine conformation of carbon rings.")
    parser.add_argument('-t', '--template', required=True, type=str,
                        help='path to direcotory containing tempaltes')
    parser.add_argument('-i', '--input', required=True, type=str,
                        help='path to direcotory containing molecules to be tested')
    parser.add_argument('-e', '--executable', required=True, type=str,
                        help='path to SiteBinder executable')
    parser.add_argument('-m', '--use_mono', action='store_true',
                        help='use mono - use this flag under Linux')
    parser.add_argument('-s', '--tolerance', type=float, default=None,
                        help='highest accaptable RMSD value')
    parser.add_argument('-A', '--all', action='store_true',
                        help='display all variants of output')
    parser.add_argument('-L', '--list', action='store_true',
                        help='display list of molecules')
    parser.add_argument('-S', '--summary', action='store_true',
                        help='display overall statistics')
    parser.add_argument('-R', '--RMSD_chart', action='store_true',
                        help='display RMSD chart in CSV format')
    args = parser.parse_args()
    template_location = args.template
    molecule_location = args.input
    sitebinder_location = args.executable
    use_mono = args.use_mono
    tolerance = args.tolerance
    display_RMSD_chart = args.RMSD_chart
    display_list = args.list
    display_summary = args.summary
    if args.all:
        display_RMSD_chart = True
        display_list = True
        display_summary = True

    conformation_list = fill_conformation_list(template_location)
    molecules = process_files(sitebinder_location, molecule_location,
                    conformation_list, use_mono, tolerance)
    if display_RMSD_chart:
        generate_RMSD_chart(conformation_list, molecules)
        sys.stdout.write('\n')
    if display_list:
        generate_result_list(molecules)
        sys.stdout.write('\n')
    if display_summary:
        generate_summary(conformation_list)
        sys.stdout.write('\n')

if __name__ == "__main__":
    main()
