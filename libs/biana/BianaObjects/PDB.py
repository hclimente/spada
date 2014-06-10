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

from Sequence import ProteinSequence, RNASequence, DNASequence
import biana.biana_globals as biana_globals

import traceback
import math
import gzip
import sys
import os
import re


temp_acc = 0
class PDB(object):

    def __init__(self, name, resolution=None):

        self.name = name
        self.resolution = resolution

        self.num_atoms = 0

        self.chains = {}                            # A dictionary of chain objects. Key: chain_name
                                                    #                                Value: Dictionary:
                                                    #                                                 Key: residue_num
                                                    #                                                 Value: Residue objects list

        self.sorted_residues = {}                   # A list with the residues ordered by numeration. Key: chain_name
        self.sorted_chains = None

        self.C_structure = None

        self.executed_dssp = False

        #self.dssp_accessibility = None
        #self.rel_dssp_accessibility = None
        #self.dssp_ss = None


    # METHOD TO CHECK
    # BETTER IF IT WAS NOT DEPENDANT ON BIOPYTHON... TOO MANY DEPENDENCES
    def read_pdb_file(pdbFile, fragments=[], merge_fragments=False, chain_value='A'):
        """
        merge_fragments is used to reenumerate residues in fragments and consider all the fragments as a single chain

        "chain_value" is used only when merging fragments
        """

        #print "Reading %s" %pdbFile

        pdbObject = PDB(name=pdbFile)

        # specify the chains to obtain
        requested_chains = [ x.get_chain() for x in fragments ]
        
        # ATOM   2331  CE1 HIS B 154      15.127  96.397  74.300  1.00100.00           C
        # ATOM   3395  O   GLY Y 221A     22.992  34.279 -11.344  1.00 26.95      2PKA35
        # ATOM     18  CE1 TYR     2      12.696 -11.556  -8.351  1.00 45.77      1PPI 2
        # ATOM   2084  NH2 ARG A 256      -5.028  36.209  68.974  1.00 47.90           N
        # ATOM   2714  OD1 ASN B1076      59.163  -6.192  27.994  1.00 94.43           O
        # ATOM   3229  CE  LYS H 217      20.184  -2.298-100.965  1.00 28.60      1GIG33
        # ATOM   1975  O   HIS B 246      45.930  10.812  85.198
        # ATOM      4  O   GLY A  -1      57.715  -4.038 -11.555  1.00 39.51           O
        #atom_regex = re.compile("ATOM\s+(\d+)\s+(\w{1,4})\s*(\w+)\s+(\w{0,1})\s*(\-{0,1}\d+)\w*\s+(\-*\d+\.\d{3})\s*(\-*\d+\.\d{3})\s*(\-*\d+\.\d{3})\s+(\d+\.\d{2})\s*(\d*\.\d+)\s+")

        #atom_regex = re.compile("ATOM\s+(\d+)\s+(\w{1,4})\s*(\w+)\s+(\w{0,1})\s*(\-{0,1}\d+)\w*\s+(\-*\d*\.\d{3})\s*(\-*\d*\.\d{3})\s*(\-*\d*\.\d{3})\s*")
        atom_regex = re.compile("ATOM\s+(\d+)\s+(\w{1,4})\s*(\w+)\s+(\w{0,1})\s+(\-{0,1}\d+)\w*\s+(\-*\d*\.\d{3})\s*(\-*\d*\.\d{3})\s*(\-*\d*\.\d{3})\s*")
        
        
        if pdbFile.endswith(".gz"):
            pdb_fd = gzip.open(pdbFile)
        else:
            pdb_fd = open(pdbFile)


        control_num = 0

        residue_num_value = 0
        
        current_residue_num = None

        for line in pdb_fd:

            if line.startswith("ATOM"):
                m = atom_regex.match(line)

                if m:

                    chain_id = m.group(4)

                    if chain_id.strip() == '':
                        chain_id = " "

                    residue_num = int(m.group(5))

                    if len(fragments)>0:
                        if True in [ current_fragment.includes( chain=chain_id, res_num = residue_num ) for current_fragment in fragments ]:
                            selected_residue = True
                        else:
                            continue
                    if merge_fragments:
                        if residue_num != current_residue_num:
                            residue_num_value += 1
                            current_residue_num = residue_num
                        chain_id = chain_value
                    else:
                        residue_num_value = residue_num

                    
                    current_atom = PDBAtom( atom_num = int(m.group(1)),
                                            atom_type = m.group(2).strip(),  # TO CHECK!!! type==name??
                                            atom_name = m.group(2).strip(),
                                            x = float(m.group(6)),
                                            y = float(m.group(7)),
                                            z = float(m.group(8)) )
                    control_num += 1
                    pdbObject.add_atom( atomObject = current_atom, chain_name = chain_id, residue_num = residue_num_value, residue_type = m.group(3) )

                else:
                    sys.stderr.write("Alert. Following atom does not match with regular expression\n%s" %line)

        if control_num == 0:
            raise ValueError("No residues were added")
        
        return pdbObject

        # Using BioPython strategy:
        # Commented because it needed many libraries and crashed in some ocasions when running pprc

        from Bio.PDB.PDBParser import PDBParser

        parser = PDBParser()

        temp_struct = parser.get_structure(pdbFile,pdbFile)

        if len(temp_struct.get_list())>1:
                print "ALERT. PDB %s has more than a single model" %pdbFile

        pdbObject = PDB(name=pdbFile)

        # specify the chains to obtain
        requested_chains = [ x.get_chain() for x in fragments ]
        
        #if chain=='-':
        #    chain_list = temp_struct[0].get_list()
        #else:
        #    chain_list = [temp_struct[0][chain]]
        chain_list = temp_struct[0].get_list()

        residue_num_value = 0

        control_num = 0

        for current_chain in chain_list:
            
            #if current_chain.get_id() not in requested_chains:
            #    continue
            for current_residue in current_chain.get_list():
                if current_residue.get_resname() == 'HOH' or current_residue.get_resname() == 'GDP' or current_residue.get_resname()==" MG":
                    continue
                
                #Determine if residue is requested or not
                residue_num = int(current_residue.get_id()[1])

                if merge_fragments:
                    chain_id = 'A'
                else:
                    chain_id = current_chain.get_id()

                if len(fragments)>0:
                    if True in [ current_fragment.includes( chain=current_chain.get_id(), res_num = residue_num ) for current_fragment in fragments ]:
                        selected_residue = True
                    else:
                        continue

                #print "Adding residue %s from chain %s" %(current_residue.get_resname(),current_chain.get_id())

                #print "Residue selected!"      
                control_num = control_num + 1

                if merge_fragments:
                    residue_num_value += 1
                else:
                    residue_num_value = residue_num

                                
                for current_atom in current_residue.get_list():
                    current_atom = PDBAtom( atom_num = current_atom.get_serial_number(),
                                            atom_type = current_atom.get_name(),
                                            atom_name = current_atom.get_name(),
                                            x = current_atom.get_coord()[0],
                                            y = current_atom.get_coord()[1],
                                            z = current_atom.get_coord()[2])
                    pdbObject.add_atom( atomObject = current_atom, chain_name = chain_id, residue_num = residue_num_value, residue_type = current_residue.get_resname())

        if control_num == 0:
            raise ValueError("No residues were added")

        #print "Added %s residues" %control_num
        return pdbObject
    
    read_pdb_file = staticmethod(read_pdb_file)


    def get_rel_dssp_accessibility(self, chain_name):

        if self.executed_dssp is False:
            self._get_dssp_results()

        return [ self.chains[chain_name][x].dssp_rel_accessibility for x in self._get_sorted_residues(chain_name=chain_name) ]
            
        #if self.rel_dssp_accessibility is None:
            
        #    dssp_results = self._get_dssp_results()
            
        #    self.dssp_accessibility = dssp_results[0]
        #    self.rel_dssp_accessibility = dssp_results[1]
        #    self.dssp_ss = dssp_results[2]

        #return self.rel_dssp_accessibility
        
    def get_dssp_accessibility(self, chain_name):

        if self.executed_dssp is False:
            self._get_dssp_results()

        return [ self.chains[chain_name][x].dssp_accessibility for x in self._get_sorted_residues(chain_name=chain_name) ]

        #if self.dssp_accessibility is None:
            
        #    dssp_results = self._get_dssp_results()
            
        #    self.dssp_accessibility = dssp_results[0]
        #    self.rel_dssp_accessibility = dssp_results[1]
        #    self.dssp_ss = dssp_results[2]

        #return self.dssp_accessibility

    def get_dssp_ss(self, chain_name):

        if self.executed_dssp is False:
            self._get_dssp_results()

        return [ self.chains[chain_name][x].dssp_ss for x in self._get_sorted_residues(chain_name=chain_name) ]

        #if self.dssp_ss is None:
            
        #    dssp_results = self._get_dssp_results()
            
        #    self.dssp_accessibility = dssp_results[0]
        #    self.rel_dssp_accessibility = dssp_results[1]
        #    self.dssp_ss = dssp_results[2]

        #return self.dssp_ss

    def set_resolution(self, resolution):
        self.resolution = float(resolution)

    def _get_sorted_residues(self, chain_name):
        if self.sorted_residues[chain_name] is None:
            self.sorted_residues[chain_name] = self.chains[chain_name].keys()
            self.sorted_residues[chain_name].sort()
        return self.sorted_residues[chain_name]

    def _get_sorted_chains(self):
        """
        Returns a list with the chain names sorted
        """
        if self.sorted_chains is None:
            self.sorted_chains = self.chains.keys()
            self.sorted_chains.sort()
            
        return self.sorted_chains
        
    def get_chain_names(self):
        return self.chains.keys()

    def get_num_atoms(self):
        return self.num_atoms

    def get_chain_num_atoms(self, chain_name):
        temp = [ x.get_num_atoms() for x in self.chains[chain_name].values() ]
        return sum(temp)
    
    def get_residues(self,chain_name):
        return self.chains[chain_name].values()

    def get_num_residues(self,chain_name=None):
        if chain_name is None:
            return sum([len(x) for x in self.chains.values()])
        else:
            return len(self.chains[chain_name])

    def get_name(self):
        return self.name

    def add_residue(self, chain_name, residue_num=0, atoms_initial_list=[], residue_type=None,hssp_conservation=None, hssp_entropy=None, hssp_exposure=None, hssp_norm_entropy=None, hssp_variability=None, dssp=None, residue_object = None):
        """
        To add residue. Chain_name is mandatory always. residue_num and atoms_initial_list are mandatory only if residue_object is not specified

        if residue_object is specified, the rest of parameters except chain_name are not used
        """

        residue_num = int(residue_num)

        if residue_object is None:
            residue_object = PDBResidue(residue_num = residue_num, residue_type=residue_type, atoms_initial_list = atoms_initial_list,
                                        hssp_conservation=hssp_conservation, hssp_entropy=hssp_entropy, hssp_exposure=hssp_exposure,
                                        hssp_norm_entropy=hssp_norm_entropy, hssp_variability=hssp_variability, dssp=dssp)
        else:
            residue_num = residue_object.residue_num

        try:
            if self.chains[chain_name].has_key(residue_num):
                print "Trying to insert twice the same residue... error?"
                #[ self.chains[chain_name][residue_num].add_atom(atomObject) for atomObject in atoms_list ]
            else:
                self.chains[chain_name][residue_num] = residue_object
                self.sorted_residues[chain_name]=None
        except:
            self.chains[chain_name] = {residue_num: residue_object}
            self.sorted_residues[chain_name]=None
        

    def add_atom(self, chain_name, residue_num, residue_type=None, atomObject=None):
        """
        Atom object can be None in order to indicate that it must be added a residue but no atoms information
        """
        #self.atoms.append(atomObject)
        self.num_atoms += 1
        residue_num = int(residue_num)
        try:
            if self.chains[chain_name].has_key(residue_num):
                self.chains[chain_name][residue_num].add_atom(atomObject)
            else:
                self.chains[chain_name][residue_num] = PDBResidue(residue_num = residue_num, residue_type=residue_type, atoms_initial_list = [atomObject])
                self.sorted_residues[chain_name]=None
        except:
            self.chains[chain_name] = {residue_num: PDBResidue(residue_num = residue_num, residue_type=residue_type, atoms_initial_list = [atomObject])}
            self.sorted_residues[chain_name]=None


    def get_pdb_summary(self):

        return "PDB with id %s, with %s chains, %s residues and %s atoms" %(self.get_name(),
                                                                            len(self.chains),
                                                                            sum([len(x) for x in self.chains.values()]),
                                                                            self.num_atoms)
    def get_in_pdb_format(self):
        """
        Gets a string with the pdb atoms in PDB format
        """

        lines = []

        lines.append("HEADER")
        lines.append("COMPND compount")
        lines.append("SOURCE")
        lines.append("AUTHOR")

        for actual_chain in self.chains.keys():
            sorted_residues = self.chains[actual_chain].keys()
            sorted_residues.sort()
            for actual_residue in sorted_residues:
                residue_object = self.chains[actual_chain][actual_residue]
                lines.extend( ["ATOM%s  %s%s%s%s%s%s%s" %(str(actual_atom.get_num()).rjust(7),
                                                          str(actual_atom.get_name()).ljust(3),
                                                          residue_object.get_residue_type().rjust(4),
                                                          str(actual_chain).rjust(2),
                                                          str(residue_object.get_residue_num()).rjust(4),
                                                          ("%.3f" %float(actual_atom.get_x())).rjust(12),
                                                          ("%.3f" %float(actual_atom.get_y())).rjust(8),
                                                          ("%.3f" %float(actual_atom.get_z())).rjust(8)) for actual_atom in residue_object.get_atoms() ] )

        return "\n".join(lines)+"\nTER"


    def get_sequence(self, chain_name=None):
        """
        Returns a dictionary with chains as keys and sequence string as values
        """

        if chain_name is not None:
            return {chain_name: Sequence.ProteinSequence("".join( [ ProteinSequence.get_aminoacid_code_3to1(self.chains[chain_name][x].get_residue_type()) for x in self._get_sorted_residues(chain_name=chain_name)] ), sequenceID=self.name+"_"+chain_name)}
        else:
            if len(self.chains)>1:
                print "ALERT: Trying to get the sequence from a PDB with more than a single chain."

            chains = self.chains.keys()

            return dict( [(chain_name, Sequence.ProteinSequence("".join( [ ProteinSequence.get_aminoacid_code_3to1(self.chains[chain_name][x].get_residue_type()) for x in self._get_sorted_residues(chain_name=chain_name)] ), sequenceID=self.name+"_"+chain_name)) for chain_name in chains ] )
    

    def get_conservation(self, chain_name):

        return [ self.chains[chain_name][x].get_hssp_conservation() for x in self._get_sorted_residues(chain_name=chain_name) ]

    def get_entropy(self, chain_name):

        return [ self.chains[chain_name][x].get_hssp_entropy() for x in self._get_sorted_residues(chain_name=chain_name) ]

    def get_accessibility(self, chain_name):

        return [ self.chains[chain_name][x].get_hssp_accessibility() for x in self._get_sorted_residues(chain_name=chain_name) ]

    def get_norm_entropy(self, chain_name):

        return [ self.chains[chain_name][x].get_hssp_norm_entropy() for x in self._get_sorted_residues(chain_name=chain_name) ]

    def get_variability(self, chain_name):

        return [ self.chains[chain_name][x].get_hssp_variability() for x in self._get_sorted_residues(chain_name=chain_name) ]

    def get_dssp(self, chain_name):

        return "".join([ self.chains[chain_name][x].dssp for x in self._get_sorted_residues(chain_name=chain_name) ])

    def get_dssp_as_horiz(self, chain_name):
        """
        Returns the dssp result in psi-pred nomenclature
        """

        return self.get_dssp(chain_name=chain_name).replace('T','C').replace('G','C').replace('S','C').replace(' ','C').replace('B','E')



    def filter_residues(self, binary_list):
        """
        Returns a new PDB object with the residues filtered by binary_list
        
        The PDB must contain a single chain
        """
        
        if len(self.chains) != 1:
            raise ValueError("Trying to filter a PDB with a more than a single chain")

        chain_name = self.chains.keys()[0]

        new_pdb = PDB(name="filtered_%s" %self.name)
        
        sorted_residues = self._get_sorted_residues(chain_name)

        if len(sorted_residues) != len(binary_list):
            raise ValueError("PDB has a distinct number of residues (%s) than binary list (%s)" %(self.get_num_residues(chain_name),len(binary_list)))

        for actual_residue_x in xrange(len(sorted_residues)):
            if binary_list[actual_residue_x]==1:
                residue_object = self.chains[chain_name][sorted_residues[actual_residue_x]]
                new_pdb.add_residue( chain_name = chain_name, residue_object = residue_object )

        return new_pdb


    def get_contacting_pairs_indices(self, pdb2, chains1, chains2, type = "min", cutoff=5.0):
        """
        Returns the indices!!! Not the residue nums!!!!
        """
        # TODO!!!! to check because copied directly from get_distances...

        if type.lower() == "min":
            method = PDBResidue.get_min_distance
        elif type.lower() == "mean":
            method = PDBResidue.get_mean_distance
        elif type.lower() == "ca":
            method = PDBResidue.get_ca_distance
        elif type.lower() == "cb":
            method = PDBResidue.get_cb_distance
        else:
            raise "Trying to calculate distances with a method unknown (%s)" %type

        distances_list = []
        
        distances = {}
        
        for actual_chain1 in chains1:
            residues1 = self._get_sorted_residues(chain_name = actual_chain1)
            distances[actual_chain1] = {}
            distances_list = []
            for actual_chain2 in chains2:
                residues2 = pdb2._get_sorted_residues(chain_name = actual_chain2)
                contacting_pairs = []
                res_index1 = 0
                for actual_residue1 in residues1:
                    res_index2 = 0
                    for actual_residue2 in residues2:
                        if method(self.chains[actual_chain1][actual_residue1],residue2=pdb2.chains[actual_chain2][actual_residue2]) <= cutoff:
                            contacting_pairs.append((res_index1,res_index2))
                        res_index2 += 1
                    res_index1 += 1

        return contacting_pairs

    def get_distances(self, pdb2, chains1, chains2, type="min", fragments1=None, fragments2=None):
        """
        "type" can be "min", "mean", "ca" or "cb"

        "fragments" is used to determine which fragments are going to be used to calculate distances
                    TODO!!!

        For the moment, it returns a dictionary with the wcm of the different chain combinations
        """

        #chains1 = self._get_sorted_chains()
        #chains2 = pdb2._get_sorted_chains()

        global temp_acc
        temp_acc = 0
        
        try:
            import pprc
        except:
            raise ValueError("It is necessary pprc module to execute this method")

        if type.lower() == "min":
            method = PDBResidue.get_min_distance
        elif type.lower() == "mean":
            method = PDBResidue.get_mean_distance
        elif type.lower() == "ca":
            method = PDBResidue.get_ca_distance
        elif type.lower() == "cb":
            method = PDBResidue.get_cb_distance
        else:
            raise "Trying to calculate distances with a method unknown (%s)" %type

        distances_list = []
        
        distances = {}
        
        for actual_chain1 in chains1:
            residues1 = self._get_sorted_residues(chain_name = actual_chain1)
            distances[actual_chain1] = {}
            distances_list = []
            for actual_chain2 in chains2:
                #a=1
                #print "CHAIN %s - %s" %(actual_chain1, actual_chain2)
                residues2 = pdb2._get_sorted_residues(chain_name = actual_chain2)
                #counter+=len(self.chains[actual_chain1])*len(pdb2.chains[actual_chain2].keys())
                distances_list.extend( [ method(self.chains[actual_chain1][actual_residue1],residue2=pdb2.chains[actual_chain2][actual_residue2])
                                         for actual_residue1 in residues1 for actual_residue2 in residues2 ] )

                #if actual_chain1 != actual_chain2:
                #contacting_pairs = []
                #for actual_residue1 in residues1:
                #    for actual_residue2 in residues2:
                #        if method(self.chains[actual_chain1][actual_residue1],residue2=pdb2.chains[actual_chain2][actual_residue2]) < 5.0:
                #            #print "Contact between residue %s[%s] - %s[%s]: %s" %(self.chains[actual_chain1][actual_residue1].residue_num,actual_chain1,
                #            #                                                      pdb2.chains[actual_chain2][actual_residue2].residue_num,actual_chain2,
                #            #                                                      method(self.chains[actual_chain1][actual_residue1],residue2=pdb2.chains[actual_chain2][actual_residue2]))

                distances[actual_chain1][actual_chain2] = pprc.wcm(rows = self.get_num_residues(chain_name = actual_chain1),
                                                                   columns = pdb2.get_num_residues(chain_name = actual_chain2),
                                                                   values = distances_list,
                                                                   description = "Distances between cb atoms between structure %s chain %s and %s chain %s" %(self.name,
                                                                                                                                                              actual_chain1,
                                                                                                                                                              pdb2.name,
                                                                                                                                                              actual_chain2 ),
                                                                   type="DISTANCE_TYPE")
        
        return distances



    def _parse_dssp_results(self, fp):

        # relative_accessibilities = [ None for x in self.]
        # accessibilities = []
        # ss = []

        surface = {'A': 115,
                   'C': 149,
                   'D': 170,
                   'E': 207,
                   'F': 230,
                   'G': 86,
                   'H': 206,
                   'I': 187,
                   'K': 222,
                   'L': 192,
                   'M': 210,
                   'N': 184,
                   'P': 140,
                   'Q': 208,
                   'R': 263,
                   'S': 140,
                   'T': 164,
                   'V': 161,
                   'W': 269,
                   'Y': 257}

        
        start_regex = re.compile("  #  RESIDUE AA STRUCTURE BP1 BP2")
        in_results = None

        temp_seq = []
        
        for line in fp:

            if in_results is None:
                if re.search(start_regex,line):
                    in_results = 1
                continue
            
            #Skipping chain breaks?
            if line[13]=="!":
                continue

            temp_seq.append(line[13])
            acc = int(line[35:38])
            #accessibilities.append(acc)
            chain = line[11]
            res_num = int(line[5:10].strip())
            res_type = line[13]

            residue_object = self.chains[chain][res_num]

            residue_object.dssp_accessibility = acc
            if ProteinSequence.get_aminoacid_code_3to1(residue_object.get_residue_type())!=res_type and res_type!="X":
                print ProteinSequence.get_aminoacid_code_3to1(residue_object.get_residue_type())
                print res_type
                sys.stderr.write("DSSP LINE: %s" %line)
                raise ValueError("Discordance between PDB and DSSP")
            #if( float(acc)*100/surface[line[13]] > 100 ):
            #    print "How can be a percentage greater than 100? residue: %s. Accessible area: %s, total surface: %s" %(line[13],acc,surface[line[13]])
            try:
                residue_object.dssp_rel_accessibility = float(acc)*100/surface[line[13]]
                #relative_accessibilities.append(float(acc)*100/surface[line[13]])
            except:
                print "The code arrives here? if yes... why??"
                pass
                #print line
                #traceback.print_exc()
                #relative_accessibilities.append(0) #Added becuase if not it produces an error
                
            if line[16]==' ':
                residue_object.dssp_ss = 'N'
                #ss.append('N')
            else:
                residue_object.dssp_ss = line[16]
                #ss.append(line[16])

        #print "DSSP SEQ: %s" %"".join(temp_seq)
        #return (accessibilities,relative_accessibilities,ss)

    def safe_traceback(self):
	# Child processes catch exceptions so that they can exit using
        # os._exit() without fanfare.  They use this function to print
        # the traceback to stderr before dying.
        #import traceback
        sys.stderr.write("Error in child process, pid %d.\n" %os.getpid())
        sys.stderr.flush()
        traceback.print_exc()
        sys.stderr.flush()

    def _get_dssp_results(self):
        """
        Executes the program dssp to get the accessibility
        """

        import popen2
        import readline
        
        dssp_program = biana_globals.DSSP_EXEC
        args = ["--"]

        print dssp_program

        # Prepare pipes
        p_readfd, c_writefd = os.pipe()
        c_readfd, p_writefd = os.pipe()

        if os.fork():
            # Parent
            for fd in (c_readfd, c_writefd, p_writefd):
                os.close(fd)
            # Convert the pipe fd to a file object, so we can use its
            # read() method to read all data.
            fp = os.fdopen(p_readfd, 'r')
            
            #results =
            self._parse_dssp_results(fp)
            fp.close()                      # Will close p_readfd.

            #return results
            
        else:
            # Child
            try:
                if os.fork():
                    # Still the same child                    
                    os.write(p_writefd,self.get_in_pdb_format())
                else:
                    # Grandchild
                    try:
                        # Redirect the pipe to stdin.
                        os.close(0)
                        os.dup(c_readfd)
                        # Redirect stdout to the pipe.
                        os.close(1)
                        os.dup(c_writefd)

                        # Trying to close stderr
                        os.close(2)

                        # Now close unneeded descriptors.
                        for fd in (c_readfd, c_writefd, p_readfd, p_writefd):
                            os.close(fd)
                        # Finally, execute the external command.
                        print dssp_program
                        print args
                        os.execv(dssp_program, [dssp_program] + args)
                    except:
                        self.safe_traceback()
                        os._exit(127)
            except:
                self.safe_traceback()
                os._exit(127)
            else:
                os._exit(0)

        self.executed_dssp = True


class PDBResidue(object):

    def __init__(self, residue_num, residue_type, atoms_initial_list=[], hssp_conservation=None, hssp_entropy=None, hssp_exposure=None, hssp_norm_entropy=None, hssp_variability=None, dssp=None):
        self.residue_num = int(residue_num)
        self.residue_type = residue_type
        
        self.atoms = []

        for current_atom in atoms_initial_list:
            self.add_atom(current_atom)

        self.ca_index = None
        self.cb_index = None

        #hssp related
        self.hssp_conservation = hssp_conservation
        self.hssp_entropy = hssp_entropy
        self.hssp_exposure = hssp_exposure
        self.hssp_norm_entropy = hssp_norm_entropy
        self.hssp_variability = hssp_variability
        self.dssp = dssp

        self.dssp_accessibility = None
        self.dssp_rel_accessibility = None
        self.dssp_ss = None

    def get_dssp(self):
        return self.dssp

    def get_hssp_conservation(self):
        return self.hssp_conservation

    def get_hssp_entropy(self):
        return self.hssp_entropy

    def get_hssp_accessibility(self):
        return self.hssp_exposure

    def get_hssp_norm_entropy(self):
        return self.hssp_norm_entropy

    def get_hssp_variability(self):
        return self.hssp_variability

    def add_atom(self, atomObject):
        if atomObject.atom_name == "CA":
            self.ca_index = len(self.atoms)
        elif atomObject.atom_name == "CB":
            self.cb_index = len(self.atoms)
        self.atoms.append(atomObject)

    def get_number_atoms(self):
        return len(self.atoms)

    def get_residue_num(self):
        return self.residue_num

    def get_atoms(self):
        return self.atoms

    def get_residue_type(self):
        return self.residue_type

    def get_num_atoms(self):
        return len(self.atoms)

    def get_min_distance(self, residue2):

        return min([x.get_distance(y) for x in self.atoms for y in residue2.atoms])

    def get_mean_distance(self, residue2):

        raise "TODO"
        return sum([ x.get_distance(atom2=y) for x in self.atoms for y in residue2.atoms ])/(len(self.atoms)*len(residue2.atoms))

    def get_ca_distance(self, residue2):
        if self.ca_index and residue2.ca_index:
            return self.atoms[self.ca_index].get_distance(atom2 = residue2.atoms[residue2.ca_index])
        else:
            print "Trying to calculate distances between alpha atoms and coordinates are not known"
            return None

    def get_cb_distance(self, residue2):

        if self.cb_index and residue2.cb_index:
            distance= self.atoms[self.cb_index].get_distance(atom2 = residue2.atoms[residue2.cb_index])
        elif self.ca_index and residue2.cb_index:
            distance= self.atoms[self.ca_index].get_distance(atom2 = residue2.atoms[residue2.cb_index])
        elif self.cb_index and residue2.ca_index:
            distance= self.atoms[self.cb_index].get_distance(atom2 = residue2.atoms[residue2.ca_index])
        elif self.ca_index and residue2.ca_index:
            distance= self.atoms[self.ca_index].get_distance(atom2 = residue2.atoms[residue2.ca_index])
        else:
            print "Trying to calculate distances between beta atoms and coordinates are not known for cb and ca: Residue %s %s" %(self.residue_num, self.residue_type)
            distance= None

        global temp_acc

        #if distance <= 8:
        #print "COntact %s(%s) - %s(%s): %s" %(self.residue_num,self.residue_type,residue2.residue_num, residue2.residue_type, distance)

        if( distance <= 8 ):
            #print "Contact %s(%s) - %s(%s): %s" %(self.residue_num,self.residue_type,residue2.residue_num, residue2.residue_type, distance)
            temp_acc = temp_acc+1

        return distance
        

class PDBAtom(object):

    def __init__(self, atom_num, atom_type, atom_name, x, y, z):

        self.atom_num = int(atom_num)
        self.atom_type = atom_type
        self.atom_name = atom_name
        self.x = x
        self.y = y
        self.z = z

    def get_coordinates(self):
        return (self.x,self.y,self.z)

    def get_x(self):
        return self.x

    def get_y(self):
        return self.y

    def get_z(self):
        return self.z

    def get_name(self):
        return self.atom_name

    def get_type(self):
        return self.atom_type

    def get_num(self):
        return self.atom_num

    def get_distance(self, atom2):
        dx = self.x - atom2.x
        dy = self.y - atom2.y
        dz = self.z - atom2.z
        return math.sqrt(dx*dx+dy*dy+dz*dz)


class PDBFragment(object):
    """
    Class to represent a pdb fragment
    """
    
    def __init__(self, start_residue=None, end_residue=None, chain=None):


        self.start_residue = start_residue
        self.end_residue = end_residue
        self.chain = chain

    def __str__(self):
        return "PDB Fragment. Chain: %s. Start: %s. End: %s" %(self.chain, self.start_residue, self.end_residue)

    def __repr__(self):
        return self.__str__()

    def get_start_residue(self):
        return self.start_residue

    def get_end_residue(self):
        return self.end_residue

    def get_chain(self):
        return self.chain

    def includes(self, chain, res_num):
        """
        Determines if the residue belongs to this fragment or not
        """
        
        if self.chain is None or self.chain==chain:
            if self.start_residue is None and self.end_residue is None:
                return True
            elif self.start_residue is None and self.end_residue>=res_num:
                return True
            elif self.end_residue is None and self.start_residue<=res_num:
                return True

        #print "Not included"
        #print "\tChain:\t%s (%s)" %(self.chain,chain)
        #print "\Range: %s-%s (%s)" %(self.start_residue,self.end_residue,res_num)
        return False


    def fragment_parser(fragment_str, separator=';'):
    
        frag_regex = re.compile("(\w*)\:*(\-*\d*)\-*(\-*\d*)")


        splitted = fragment_str.split(separator)


        fragments = []

        for current_splitted in splitted:

            if current_splitted.strip()=='':
                continue

            if current_splitted=="-":
                current_splitted=''
                fragments.append(PDBFragment( chain = None, start_residue = None, end_residue = None ))
                continue

            m = frag_regex.match(current_splitted)

            if m:
                if m.group(1) == '':
                    chain = None
                else:
                    chain = m.group(1)

                if m.group(2)=='':
                    start_residue = None
                else:
                    start_residue = int(m.group(2))

                if m.group(3)=='':
                    end_residue = None
                else:
                    end_residue = int(m.group(3))

                fragments.append(PDBFragment( chain = chain, start_residue = start_residue, end_residue = end_residue ))
            else:
                print "Fragment %s does not match with regular expresion" %current_splitted

        return fragments

    fragment_parser = staticmethod(fragment_parser)


    
    

