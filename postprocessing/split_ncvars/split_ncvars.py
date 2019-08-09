#!/usr/local/python/2.7.3/bin/python -u
#
# Program that will split out variables
# fron netCDF files.
#
#----------------------------------------

import argparse
import netCDF4
import os
import re
import sys
import time
import shutil
import subprocess
from subprocess import PIPE

try:
    from nco import Nco
    nco = Nco()
except:
    raise Exception("ERROR: Unable to load pynco, please use python 2.7.3 on PPAN")

class OpenNetCDF():

    def __init__(self,input_file,output_dir,variables,static=False):
        """ Initialize the Open NetCDF file for processing """
        self.variables = variables.split(",")
        self.output_dir = output_dir
        self.filename = input_file
        if static:
            self.static = True
        try:
            self.input_file = netCDF4.Dataset(input_file, "r")
        except:
            raise Exception("ERROR: Unable to open file %s" % self.filename)
        self.found_variables = self.input_file.variables

    def extract_variables(self):
        """ Gets the list of available variables and pull them out of
        the netCDF file that is being worked on."""
        for variable in self.variables:
            cell_methods_dict = {}
            extract_variables = []
            extract_variables.append(variable)
            if variable not in self.found_variables:
                raise Exception("ERROR: Requested variable %s was not found in input file" % variable)
            metadata = self.input_file.variables["%s" % variable]
            if "time" not in metadata.dimensions:
                """ We have a static variable """
                if len(metadata.dimensions) == 1:
                    self.run_split_command(extract_variables,variable)
                    self.move_output_file(variable)
                else:
                    cell_measures_dict = self.get_cell_measures_dict(metadata, variable)
                    if cell_measures_dict:
                        [extract_variables.append(value) for value in cell_measures_dict.values() if self.double_check_var_presence(value)]
                    variable_coords = self.get_variable_coords(metadata, variable)
                    if variable_coords:
                        [extract_variables.append(coord) for coord in variable_coords if self.double_check_var_presence(coord)]
                    extra_vars = self.check_for_extra_vars(metadata.dimensions,variable)
                    if extra_vars:
                        [extract_variables.append(extra_var) for extra_var in extra_vars if self.double_check_var_presence(extra_var)]
                    self.run_split_command(extract_variables,variable)
                    self.move_output_file(variable)
            else:
                """ This variable is time varying """
                cell_measures_dict = self.get_cell_measures_dict(metadata, variable)
                if cell_measures_dict:
                    [extract_variables.append(value) for value in cell_measures_dict.values() if self.double_check_var_presence(value)]
                variable_coords = self.get_variable_coords(metadata, variable)
                if variable_coords:
                    [extract_variables.append(coord) for coord in variable_coords if self.double_check_var_presence(coord)]
                time_average_info = self.get_time_average_info(metadata, variable)
                if time_average_info:
                    [extract_variables.append(time_avg) for time_avg in time_average_info if self.double_check_var_presence(time_avg)]
                try:
                    time_data = self.input_file.variables["time"]
                    extract_variables.append(time_data.bounds)
                    extra_vars = self.check_for_extra_vars(metadata.dimensions,variable)
                    if extra_vars:
                        [extract_variables.append(extra_var) for extra_var in extra_vars if self.double_check_var_presence(extra_var)]
                except:
                    """ No time bounds found in the file """
                self.run_split_command(extract_variables,variable)
                self.move_output_file(variable)

    def move_output_file(self,variable):
        """ move the temporary file and append to file if it currently exists """
        outfile = os.path.join(self.output_dir,".".join([".temp",variable,"nc"]))
        destination_outfile = os.path.join(self.output_dir,".".join([variable,"nc"]))
        tmpfile = os.path.join(self.output_dir,".".join(["","tmp","nc"]))
        if os.path.isfile(destination_outfile):
            """Assuming that this it the files are sequential at this point"""
            shutil.move(destination_outfile,tmpfile)
            nco.ncrcat(input=[tmpfile,outfile], output=destination_outfile, options='-h')
        else:
            """This file obviously doesn't exist, so we need to simply move the file"""
            shutil.move(outfile,destination_outfile)
            nco.ncatted(input=destination_outfile, options=" ".join(["-h","-O","-a","\"filename,global,m,c,"+variable+".nc\""]))
        if os.path.isfile(tmpfile):
            os.remove(tmpfile)
        if os.path.isfile(outfile):
            os.remove(outfile)

    def run_split_command(self,extract_variables,variable):
        """ Portion of the program that runs the split_commands """
        ncks_options = " ".join(["-h","--64bit","--header_pad","16384"])
        nco.ncks(variable=extract_variables, input=self.filename, 
            output=os.path.join(self.output_dir,".".join([".temp",variable,"nc"])), options=ncks_options)

    def double_check_var_presence(self,varname):
        """ Routine to double check that we have the variables
        that we are attempting to split out before running the split command """
        if varname in self.found_variables:
            return True
        else:
            return False

    def check_for_extra_vars(self,dimensions,variable):
        """ Routine to check for third dimension edges """
        extra_vars = []
        for dimension in dimensions:
            dim_metadata = self.input_file.variables["%s" % dimension]
            if re.search("(?i)^p(lev|half|full)", dimension):
                [extra_vars.append(variable) for variable in ["bk","pk","ps"]]
            try:
                extra_vars.append(dim_metadata.edges)
            except:
                """ No edge information associated with the variable """
                
        return extra_vars
                             
    def get_time_average_info(self,metadata,variable):
        """ Routine to extract time average information associate with 
        a given variable, returned list associated will be empty if 
        no information is found."""
        try:
            variable_time_avg_str = metadata.time_avg_info
            return variable_time_avg_str.split(",")
        except:
            """print("INFO: No time average information found in association with %s" % variable)"""
        return 

    def get_variable_coords(self,metadata,variable):
        """ Routine to extract variable coordinates from the 
        netCDF file associated with the given variable."""
        variable_coords_list = []
        try:
            variable_coords_str = metadata.coordinates
            return variable_coords_str.split(" ")
        except:
            """print("INFO: Requested variable %s has no associated coordinates")"""
            return 

    def get_cell_measures_dict(self,metadata,variable):
        """ Routine to get the cell_measures_dictionary,
        or return empty dicitonary if there are no 
        cell_measrues associated with the requested variable."""
        cell_measures_dict = {}
        try:
            cell_measures_dict = self.format_cell_methods_and_measures(metadata.cell_measures)
        except:
            """print("INFO: No cell measures attribute associated with %s" % variable)"""
        return cell_measures_dict

    def format_cell_methods_and_measures(self,variable_string):
        """ Routine to change string from netCDF variable 
        cell_measures and cell_methods to a python dict."""
        attribute_dictionary = {}
        variable_string_bits = variable_string.split(" ")
        for i in range(0, len(variable_string_bits), 2):
            key = re.sub(":$", "", variable_string_bits[i])
            value = variable_string_bits[i+1]
            attribute_dictionary[key] = value
        return attribute_dictionary

    def __del__(self):
        self.input_file.close()


if __name__ == "__main__":

    print("DEPRECATION WARNING. UNSUPPORTED TOOL, INVOKING split_ncvars.pl")
    os.execvp('split_ncvars.pl', sys.argv)

    parser = argparse.ArgumentParser(description="Program to split variables and related fields from netCDF files")
    parser.add_argument("-i", "--input", help="directory to look for netCDF files", required=True)
    parser.add_argument("-o", "--output", help="directory to save newly created netCDF files", required=True)
    parser.add_argument("-v", "--variables", help="comma separated list of variables to split out of a given netCDF files", required=True)
    parser.add_argument("-f", "--file", help="filename to save all variables in a single netCDF file", required=False)
    parser.add_argument("-d", "--debug", action="store_true", help="more verbose output from the tool", required=False)
    parser.add_argument("-s", "--static", action="store_true", help="let the tool know you're extract static variables and exclude check for third dimension", required=False, default=False)
    parser.add_argument("-t", "--tiled", action="store_true", help="flag to decide if files on commandline are tiled or timevarying", required=False)
    parser.add_argument("files", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if len(args.files) == 0:
        raise Exception("ERROR: no files supplied on the commandline")

    for filename in args.files:
        NetCDFDoc = OpenNetCDF(os.path.join(args.input,filename),args.output,args.variables,args.static)
        NetCDFDoc.extract_variables()
    exit()
