"""
    BIANA: Biologic Interactions and Network Analysis
    Copyright (C) 2009  Javier Garcia-Garcia, Emre Guney, Baldo Oliva

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""


"""
File        : bianaParser.py
Author      : Javier Garcia
Creation    : October 2007
Contents    : General parser to introduce information into biana
Called from : command line

=======================================================================================================

This file implements a program that fills up tables in database biana with information from distinct databases

"""

## STEP 1: IMPORT NECESSARY MODULES

import sys
import getopt
import re
import time
import gzip
import traceback
import os
#import tarfile


from biana.BianaDB import BianaDBaccess
from biana.BianaObjects import *


class BianaParser(object):
    """
    General Parser Class to biana
    """

    def __init__(self, default_db_description = None,
                 default_script_name = "bianaParser.py",
                 default_script_description = "This file implements a program that fills up tables in database biana with information from distinct databases",
                 #content_type_list = [],
                 additional_compulsory_arguments = [],
                 additional_optional_arguments = []):
        
        """
        Starts the bianaParser Object
        """

        print "Parser object started"

        self.compulsory_arguments = [ ("input-identifier=",None,"path or file name of input file(s) containing database data. Path names must end with \"/\"."),
                                      ("biana-dbname=",None,"name of database biana to be used"),
                                      ("biana-dbhost=",None,"name of host where database biana to be used is placed"),
                                      ("database-name=",None,"internal identifier name to this database (it must be unique in the database)"),
                                      ("database-version=",None,"version of the database to be inserted") ]

        self.compulsory_arguments.extend(additional_compulsory_arguments)
        

        self.optional_arguments = [ ("biana-dbuser=",None,"username accessing the database (not required in most systems)"),
                                    ("biana-dbpass=",None,"password of username accessing the database (not required in most systems"),
                                    ("help",None,"prints this message and exits"),
                                    ("verbose",0,"prints process info to stdout"),
                                    ("log-file=",None,"Prints a log file of the parsing result (number of inserted proteins, references...)"),
                                    ("time-control",None,"prints to stderr a control of the timing of the parser"),
                                    ("database-description=",default_db_description,"Description of the database to be inserted."),
                                    ("optimize-for-parsing",None,"Optimizes database for parsing"),
				    ("promiscuous",False,"sets the database to be parsed as promiscuous (whose entities can be included in multi user entities)") ]
                                    #("mode=","scratch","sets mode to be used by parser. Valid modes are: \"scratch\" (biana database is empty, create it from scratch) or \"tables\" (fill only tables indicated in tables_to_fill (see code)")]   
                                           
        self.optional_arguments.extend(additional_optional_arguments)

        self.script_name = default_script_name
        self.script_description = default_script_description

        #Parse general methods
        self.arguments_dic = self.parseArguments()
        self.input_file = self.arguments_dic["input-identifier"]
        self.biana_dbname = self.arguments_dic["biana-dbname"]
        self.biana_dbhost = self.arguments_dic["biana-dbhost"]
        self.sourcedb_name = self.arguments_dic["database-name"]
        self.sourcedb_version = self.arguments_dic["database-version"]
        self.biana_dbuser = self.arguments_dic["biana-dbuser"]
        self.biana_dbpass = self.arguments_dic["biana-dbpass"]
        self.help = self.arguments_dic["help"]
        self.verbose = self.arguments_dic["verbose"]
        self.time_control = self.arguments_dic["time-control"]
        self.log_file = self.arguments_dic["log-file"]
        self.optimize_for_parsing = self.arguments_dic["optimize-for-parsing"]
        #self.mode = self.arguments_dic["mode"]
	self.is_promiscuous = self.arguments_dic["promiscuous"] # Flag deciding whether database gives information that is going to be added to more than one user entiries

        self.database = None
        if self.arguments_dic.has_key("default-attribute"):
            self.default_eE_attribute = self.arguments_dic["default-attribute"] # default externalEntityAttribute specified by the particular database parser (it will be overwritten in the parser if not given as argument)
        else:
            self.default_eE_attribute = ""
        
        #self.content_type_list = content_types


    def start(self):


        print "Parser started"
        if isinstance(self.sourcedb_name,int) or isinstance(self.sourcedb_version,int):
            sys.stderr.write("You must insert correctly the database name and database version\n")
            sys.exit(1)
            
        #if( self.mode=="scratch" ):
        self.database_description = self.arguments_dic["database-description"]

        # Log dictionary where all log information will be stored
        self.log = {}
        if self.log_file:
            self.log_file_fd = file(self.log_file, 'w')

        self.biana_access = BianaDBaccess(dbname=self.biana_dbname, dbhost=self.biana_dbhost, dbuser=self.biana_dbuser, use_buffer=True, dbpassword=self.biana_dbpass, lock_tables=True, check_integrity=True )


        # check data consistency

        # Time related
        self.initial_time = time.time()

        # Insert the information associated to the parsed database
        # Introduce database info into biana database
        #if( self.mode=="scratch" ):

        self.database = ExternalDatabase( databaseName = self.sourcedb_name,
                                          databaseVersion = self.sourcedb_version,
                                          databaseFile = self.input_file.split(os.sep)[-1],
                                          databaseDescription = self.database_description,
                                          defaultExternalEntityAttribute = self.default_eE_attribute,
					  isPromiscuous = self.is_promiscuous ) 
                                          #content_type_list = self.content_type_list)

        self.biana_access.insert_new_external_database( externalDatabase = self.database )
                                                               
        # Open the input file descriptor
        # This is a responsability of subclasses method
            
        try:
            if self.optimize_for_parsing:
            	self.biana_access.optimize_database_for(mode="parsing")

            self.parse_database()
            
            # set the parsing time
            self.database.set_parsing_time( int(time.time() - self.initial_time) )

            # Updates the information that this external database has inserted
            self.biana_access.update_external_database_external_entity_attributes( self.database )

            self.close()


        except:
            traceback.print_exc()
            sys.stderr.write("ERROR WHILE PARSING. ALL MODIFICATIONS ARE GOING TO BE DELETED\n")
            self.biana_access._rollback()
            sys.exit(1)
        

    # METHODS

    def close(self):
        ## LAST STEP: CLOSE DATABASE CONNECTION    IMPORTANT !!!!
        ## As bianaDBaccess uses an internal buffer, it is necessary to close the connection to sure that all inserts are correctly done, as well as unlock tables


        self.biana_access.close()
        
        if self.time_control:
            sys.stderr.write("Total time: %s seconds\n" %(time.time()-self.initial_time))

        if self.log_file:
            self.log_file_fd.write(self.get_log_string())
            self.log_file_fd.close()

        if self.verbose:
            sys.stderr.write("\n Total time: %s \n" %(time.time()-self.initial_time) )
            sys.stderr.write(self.get_log_string())



    ## GENERAL PARSER METHODS ##
    def parseArguments(self):
        """
        Method that returns a dictionary with the values of the arguments
        
        """
        
        arguments = self.compulsory_arguments+self.optional_arguments
        
        # Set all default values
        #return_values = [i[1] for i in arguments]
        return_dict = {}
        for i in arguments:
            return_dict[i[0].replace("=","")] = i[1]
        
        # Obtain a list with the names of all arguments
        list_arguments = [argument[0] for argument in arguments]
        # It can be of the following way because it contains "=" digit
        #list_arguments = return_dict.keys()
        
        
        # Parse arguments
        try:
            opts, args = getopt.getopt(sys.argv[2:], "", list_arguments)
            
        except getopt.GetoptError, bad_opt:
            # return error in parsing parameters, and return void list
            raise ValueError("%s\n" %(bad_opt.__str__()) )

        # If there is no error, continue with the parsing
        for option,value in opts:
            if option=="--help":
                self.print_help()
                sys.exit(2)
            for actual_argument  in list_arguments:
                # Delete the "=" value if it has
                temp_arg = actual_argument.replace("=","")
                if option=="--"+temp_arg:
                    if value=="":
                        return_dict[temp_arg]=1
                    else:
                        return_dict[temp_arg]=value

        
        # Check for all compulsory arguments:
        for comp_arg in self.compulsory_arguments:
            if return_dict[comp_arg[0].replace("=","")] is None:
                sys.stderr.write("%s argument is not defined!\n" %(comp_arg[0].replace("=","")))
                self.print_help()
                sys.exit(2)

        return return_dict


    ## LOG RELATED METHODS ##
    def add_to_log(self,key):
        """
        Increment the counter of a log dictionary for a given key
        
        Used in parsers
        
        """

        try:
            self.log[key] += 1
        except KeyError:
            self.log[key] = 1

            
    def get_log_string(self):
        """
        Returns a string with the content of the log dictionary
        
        Format: key: value
        """
        
        string_list = []

        for log_element in self.log.keys():
            string_list.append("%s: %s" %(log_element,self.log[log_element]))

        return "\n".join(string_list)
    

    def print_help(self):

        print "--------------------------------------------------------------------------------------------------------------"
        print "DESCRIPTION:"
        print "\t"+self.script_description

        usage = "\tpython %s " %(self.script_name)

        for argument in self.compulsory_arguments:
            if re.search("=",argument[0]):
                usage = usage + "--%s%s " %(argument[0],argument[0].rstrip("="))
            else:
                usage = usage + "--%s " %(argument[0])

        for argument in self.optional_arguments:
            if re.search("=",argument[0]):
                usage = usage + "[--%s=%s] " %(argument[0],argument[0].rstrip("="))
            else:
                usage = usage + "[--%s] " %(argument[0])

        print "\n"
        print "USAGE:"
        print usage

        print "\nWHERE:\n"

        if len(self.compulsory_arguments)>0:
            print "COMPULSORY ARGUMENTS:"

            for argument in self.compulsory_arguments:
                sys.stdout.write("\t%s:" %(argument[0].rstrip("=")))
                sys.stdout.write("%s" %(self._indent(3,len(argument[0])-1)))
                argument_description = self._splitsize(string=argument[2],size=80)
                if len(argument_description)==1:
                    sys.stdout.write("%s\n" %(argument[2]))
                else:
                    sys.stdout.write("%s\n" %(argument_description[0]))
                    for i in range(1,len(argument_description)):
                        sys.stdout.write("\t\t\t\t%s\n" %(argument_description[i]))

        if len(self.optional_arguments)>0:
            print
            print "OPTIONAL ARGUMENTS:"

            for argument in self.optional_arguments:
                sys.stdout.write("\t%s:" %(argument[0].rstrip("=")))
                sys.stdout.write("%s" %(self._indent(3,len(argument[0])-1)))
                argument_description = self._splitsize(string=argument[2]+" [default: %s]" %(argument[1]),size=80)
                if len(argument_description)==1:
                    sys.stdout.write("%s [default: %s]\n" %(argument[2],argument[1]))
                else:
                    sys.stdout.write("%s\n" %(argument_description[0]))
                    for i in range(1,len(argument_description)):
                        sys.stdout.write("\t\t\t\t%s\n" %(argument_description[i]))
                
        print "--------------------------------------------------------------------------------------------------------------"



    def _indent(self,max_num_tabulators, initial_length):

        num_tabulators = max_num_tabulators - (initial_length+1)/8

        values_to_return = []

        #print "num tabulators: %s" %(num_tabulators)

        for i in xrange(num_tabulators):
            values_to_return.append("\t")

        return "".join(values_to_return)

    def _splitsize(self, string, size):
        """
        Split a string in substrings with a determined size
        """

        list_return = []

        final_position=0

        if len(string)<=size:
            list_return = [string]
        else:
            for i in xrange(len(string)/size):
                initial_position = i*size + final_position - i*size
                final_position = (i+1)*size
                while final_position<len(string) and string[final_position] != " " and  string[final_position] != "\t":
                    final_position += 1
                list_return.append(string[initial_position:final_position])

        return list_return


    def initialize_input_file_descriptor(self):
        """
        Create the input file descriptor given input database file name. Handles gzipped data as well.
        """

        self.input_file_fd = None

        if self.input_file != "None":

	    if not os.path.exists(self.input_file):
		raise ValueError("File or directory %s is missing" % (self.input_file))

            if( os.path.isfile(self.input_file) ):
                #if( self.input_file.endswith("tar.gz") ):
                #    self.input_file_fd = tarfile.open(self.input_file,'r')
                if( self.input_file.endswith(".gz") ):
                    self.input_file_fd = gzip.open(self.input_file,'r')
                else:
                    self.input_file_fd = file(self.input_file, 'r')
            elif( os.path.isdir(self.input_file) ):
                self.input_file_fd = None


    def parse_database(self):
        """
        Method to be overwritten by specific parsers

        The method must include the calls to control lock and unlock database procedures
        """
        return


