import sys
import os


# A very dirty way of doing, but... if anyone knows how to do it correctly, tell me

biana_path = None
site_packages_path = None

# get where the source code, looking subdirectories in the python system path
for current_path_index in xrange(1,len(sys.path)):
    if not os.path.exists(sys.path[current_path_index]):
	continue
    if not os.path.isdir(sys.path[current_path_index]):
        continue
    for current_dir in os.listdir(sys.path[current_path_index]):
	if not os.path.isdir(sys.path[current_path_index]+os.sep+current_dir):
	    continue
	if current_dir.lower() == "biana":
	    biana_path = sys.path[current_path_index]
	    break
        
if biana_path is None:
    biana_path = site_packages_path

biana_path = biana_path+os.sep+"biana"+os.sep+"BianaParser" #+os.sep

sys.path.append(biana_path)

list_all_files = os.listdir(biana_path)

not_allowed_import_files = ["bianaParser.py","database2biana.py","__init__.py", "psi_MiXMLParser.py"]

list_python_files = []
for current_file in list_all_files:
    if current_file.endswith(".py") and current_file not in not_allowed_import_files and not current_file.startswith("."):
        list_python_files.append(current_file)

parserObjects = {}

# Import all the parsers
imports = map(__import__, [ x[0:-3] for x in list_python_files] )

# Get the BianaParser class
for current_parser in imports[0].BianaParser.__subclasses__():
    parserObjects[current_parser.name.lower()] = current_parser


def run():
    """
    Runs the script to insert the parser to the database
    """
    parser = ""
    try:
        parser = sys.argv[1]
    except:
        sys.stderr.write("First argument must be the name of the parser to execute.\n")

    if parserObjects.has_key(parser.lower()):
        parserObjects[parser.lower()]().start()
    else:
        sys.stderr.write("Database parser not recognized\n")
        parsers_list = [x[0] for x in get_available_parsers_info()]
        parsers_list.sort()
        sys.stderr.write("Accepted parsers are: \n\t%s\n" %"\n\t".join(parsers_list))
        sys.exit(2)


def get_available_parsers_info():
    """
    Returns a list of tuples with the information of the parsers. (parser_name, parser_description, external_entity_definition)
    """

    return [ (current_parser.name, current_parser.description, current_parser.external_entity_definition) for current_parser in parserObjects.values() ]
        


