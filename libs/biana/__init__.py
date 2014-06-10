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

# Graphical interface communication
from OutBianaInterface import OutBianaInterface

# Administration commands
from biana_commands import administration

# Session commands
from biana_commands import create_new_session
from biana_commands import load_session
from biana_commands import save_session
from biana_commands import remove_session
from biana_commands import available_sessions #

from biana_commands import close
from biana_commands import ping

# Utility commands
#import biana_globals
import utilities

# Delete some commands to not appear when importing biana and making dir()
del biana_commands


credits = "Structural Bioinformatics Lab"

copyright = "Copylefted :P"

license = "GPL"

help = "go to SBI web page http://sbi.imim.es/"



## Add the external packages to the path ##

import sys
import os

biana_path = None
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

# Add the new biana path to the begining of the path, to have preference over other versions of mysqldb or networkx
if biana_path is not None:
	biana_path = biana_path+os.sep+"biana"+os.sep+"ext"
	new_sys_path = [biana_path]
	new_sys_path.extend(sys.path)
	sys.path = new_sys_path 

del biana_path
del sys
del os
