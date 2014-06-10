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

import md5
from biana.utilities import autoincrement


class Sequence(object):
    """
    Class to represent a biologic sequencial entity, as nucleotide sequences, protein sequences, etc

    This class is suposed to be an abstract class. Only instances for their subclasses should be created
    """


    def __init__(self, sequence, sequenceMD5=None, sequenceID=None, sequence_type=None):
        """
        "sequence": the sequence itself. It is processed to eliminate non sequence elements

        "sequenceMD5" is the digested md5 code of the sequence. It usually has not to be defined as a parameter, as it implements a method to calculate it

        "sequenceID" is the unique identifier for the sequence. A sequence must have a sequenceID only when it has been inserted into database

        "sequence_type" is the type of sequence (protein, dna or rna)
        """


        self.sequence = sequence.replace(" ", "").replace("*", "").replace("\n", "").replace("\t", "").replace("\r", "").replace("_", "").upper()
        self.sequenceMD5_value = sequenceMD5

        self.sequenceID = sequenceID
        self.sequence_type = sequence_type
        
        self.fragmented_sequence = None

    def __str__(self):
        return self.sequence

    def get_length(self):
        return len(self.sequence)

    def get_sequence(self):
        return self.sequence

    def get_type(self):
        return self.sequence_type

    def _get_fragmented_sequence(self, using_size, using_translation_method, using_list):
        """
        it returns a list with the ordered indices of the fragments
        """

        divisions = len(self.sequence) / using_size
        resta = len(self.sequence) % using_size

        temp_seq = self.sequence
        
        if resta!=0:
            divisions += 1
            for less in xrange(using_size-resta):
                temp_seq += using_list[0]
            del less

        return [ using_translation_method(temp_seq[(x*using_size):(x*using_size)+using_size]) for x in xrange(divisions) ]

    def get_fragmented_sequence(self):
        return self.fragmented_sequence

    def get_sequence_MD5(self):
        """
        Now:
        SequenceMD5 is used only for internal working, and it is represented in a 16-byte string

        Previous:
        Return MD5 code for sequence "sequence"
        (MD5 hexdigestion of sequence + its leading 4 chars
        + its last 4 chars)
        """

        if self.sequenceMD5_value is None:
            sequence = self.sequence.strip()
            toconvert = md5.new(sequence)
            #digested = toconvert.digest().replace('\\','\\\\').replace('"','\\"').replace("'","\\'")
            digested = toconvert.digest()
            self.sequenceMD5_value = digested
        
        #if self.sequenceMD5_value is None:
        #    sequence = self.sequence.strip()
        #    head = sequence[:4]
        #    tail = sequence[-4:]
        #    toconvert = md5.new(sequence)
        #    digested = toconvert.hexdigest()
        #    self.sequenceMD5_value = digested + head + tail

        return self.sequenceMD5_value


    def get_sequenceID(self):
        return self.sequenceID

    def set_sequenceID(self, sequenceID, force=False):
        """
        Sets the sequenceID for a sequence.

        If it contains a sequence ID it should not be modified. In special cases where it should be able to be modified, use the parameter force=True
        """

        if( self.sequenceID is None ) or force is True:
            self.sequenceID = sequenceID
        else:
            raise ValueError("Trying to set an ID to a sequence that previously had one...)")


    def read_sequences_file(inputPath, format="fasta", sequences_type="protein"):
        """
        Returns a list of sequence objects contained in a file in the specified format

        "format" can be: fasta

        "sequences_type" can be: protein, dna or rna
        """

        if inputPath.endswith(".gz"):
            import gzip
            in_fd = gzip.open(inputPath,'r')
        else:
            in_fd = open(inputPath)

        sequences_list = []

        temp_seq = []

        format = format.lower()
        sequences_type = sequences_type.lower()

        sequenceID = None

        if format=="fasta":
            for line in in_fd:
                if line[0]==">":
                    if len(temp_seq)>0:
                        if sequences_type=="protein":
                            sequences_list.append( ProteinSequence( sequence = "".join(temp_seq), sequenceID = sequenceID ) )
                        elif sequences_type=="dna":
                            sequences_list.append( DNASequence( sequence = "".join(temp_seq), sequenceID = sequenceID ) )
                        elif sequences_type=="rna":
                            sequences_list.append( RNASequence( sequence = "".join(temp_seq), sequenceID = sequenceID ) )
                        else:
                            raise ValueError('Sequence type not recognized: %s' %sequences_type)
                    temp_seq = []
                    sequenceID=line[1:].strip()
                else:
                    temp_seq.append(line.strip())

        else:
            raise ValueError("Format not recognized: %s" %(format))

        if len(temp_seq)>0:
            if sequences_type=="protein":
                sequences_list.append( ProteinSequence( sequence = "".join(temp_seq), sequenceID = sequenceID ) )
            elif sequences_type=="dna":
                sequences_list.append( DNASequence( sequence = "".join(temp_seq), sequenceID = sequenceID ) )
            elif sequences_type=="rna":
                sequences_list.append( RNASequence( sequence = "".join(temp_seq), sequenceID = sequenceID ) )
            else:
                raise ValueError('Sequence type not recognized: %s' %sequences_type)
            
        return sequences_list

    read_sequences_file = staticmethod( read_sequences_file )


class RNASequence(Sequence):

    translation_dict = None
    dictionary = "ACUG"
    window_size = 8

    def __init__(self, sequence, sequenceMD5=None, sequenceID=None):
        """
        "sequence": the sequence itself. It is processed...
        """

        Sequence.__init__( self, sequence = sequence, sequenceMD5 = sequenceMD5, sequenceID = sequenceID, sequence_type="rna" )


    def _get_digested_code(self, window):
        """
        """

        if RNASequence.translation_dict is None:
            autoinc = autoincrement(1)
            RNASequence.translation_dict = dict( [ ("%s%s%s%s%s%s%s%s" %(a,b,c,d,e,f,g,h),autoinc.next()) for a in RNASequence.dictionary
                                           for b in RNASequence.dictionary for c in RNASequence.dictionary
                                           for d in RNASequence.dictionary for e in RNASequence.dictionary
                                           for f in RNASequence.dictionary for g in RNASequence.dictionary for h in RNASequence.dictionary ] )
        
        try:
            return RNASequence.translation_dict[window]
        except KeyError:
            return 0


    def _get_fragmented_sequence(self):

        return Sequence._get_fragmented_sequence(self, 
                                                 using_translation_method = self._get_digested_code,
                                                 using_size = RNASequence.window_size,
                                                 using_list = RNASequence.dictionary)
    

class DNASequence(Sequence):

    translation_dict = None
    dictionary = "ACTG"
    window_size = 8

    def __init__(self, sequence, sequenceMD5=None, sequenceID=None):
        """
        "sequence": the sequence itself. It is processed...
        """

        # test to summarize the sequence as a list of numbers

        Sequence.__init__( self, sequence = sequence, sequenceMD5 = sequenceMD5, sequenceID = sequenceID, sequence_type="dna" )



    def _get_digested_code(self, window):
        """
        """

        if DNASequence.translation_dict is None:
            autoinc = autoincrement(1)
            DNASequence.translation_dict = dict( [ ("%s%s%s%s%s%s%s%s" %(a,b,c,d,e,f,g,h),autoinc.next()) for a in DNASequence.dictionary
                                           for b in DNASequence.dictionary for c in DNASequence.dictionary
                                           for d in DNASequence.dictionary for e in DNASequence.dictionary
                                           for f in DNASequence.dictionary for g in DNASequence.dictionary for h in DNASequence.dictionary ] )
        
        try:
            return DNASequence.translation_dict[window]
        except KeyError:
            return 0        


    def _get_fragmented_sequence(self):

        return Sequence._get_fragmented_sequence(self, 
                                                 using_translation_method = self._get_digested_code,
                                                 using_size = DNASequence.window_size,
                                                 using_list = DNASequence.dictionary )


class ProteinSequence(Sequence):

    translation_dict = None
    dictionary = "GALMFWKQESPVICYHRNDTXZBUJ"
    window_size = 3

    aminoacid_codes_3to1 = { "GLY": 'G',
                             "ALA": 'A',
                             "LEU": 'L',
                             "MET": 'M',
                             "PHE": 'F',
                             "TRP": 'W',
                             "LYS": 'K',
                             "GLN": 'Q',
                             "GLU": 'E',
                             "SER": 'S',
                             "PRO": 'P',
                             "VAL": 'V',
                             "ILE": 'I',
                             "CYS": 'C',
                             "TYR": 'Y',
                             "HIS": 'H',
                             "ARG": 'R',
                             "ASN": 'N',
                             "ASP": 'D',
                             "THR": 'T',
                             "MSE": 'M', #Selenomethionin
                             "CSS": 'C', 
                             '2AS':'D', '3AH':'H', '5HP':'E', 'ACL':'R', 'AIB':'A',
                             'ALM':'A', 'ALO':'T', 'ALY':'K', 'ARM':'R', 'ASA':'D',
                             'ASB':'D', 'ASK':'D', 'ASL':'D', 'ASQ':'D', 'AYA':'A',
                             'BCS':'C', 'BHD':'D', 'BMT':'T', 'BNN':'A', 'BUC':'C',
                             'BUG':'L', 'C5C':'C', 'C6C':'C', 'CCS':'C', 'CEA':'C',
                             'CHG':'A', 'CLE':'L', 'CME':'C', 'CSD':'A', 'CSO':'C',
                             'CSP':'C', 'CSS':'C', 'CSW':'C', 'CXM':'M', 'CY1':'C',
                             'CY3':'C', 'CYG':'C', 'CYM':'C', 'CYQ':'C', 'DAH':'F',
                             'DAL':'A', 'DAR':'R', 'DAS':'D', 'DCY':'C', 'DGL':'E',
                             'DGN':'Q', 'DHA':'A', 'DHI':'H', 'DIL':'I', 'DIV':'V',
                             'DLE':'L', 'DLY':'K', 'DNP':'A', 'DPN':'F', 'DPR':'P',
                             'DSN':'S', 'DSP':'D', 'DTH':'T', 'DTR':'W', 'DTY':'Y',
                             'DVA':'V', 'EFC':'C', 'FLA':'A', 'FME':'M', 'GGL':'E',
                             'GLZ':'G', 'GMA':'E', 'GSC':'G', 'HAC':'A', 'HAR':'R',
                             'HIC':'H', 'HIP':'H', 'HMR':'R', 'HPQ':'F', 'HTR':'W',
                             'HYP':'P', 'IIL':'I', 'IYR':'Y', 'KCX':'K', 'LLP':'K',
                             'LLY':'K', 'LTR':'W', 'LYM':'K', 'LYZ':'K', 'MAA':'A',
                             'MEN':'N', 'MHS':'H', 'MIS':'S', 'MLE':'L', 'MPQ':'G',
                             'MSA':'G', 'MSE':'M', 'MVA':'V', 'NEM':'H', 'NEP':'H',
                             'NLE':'L', 'NLN':'L', 'NLP':'L', 'NMC':'G', 'OAS':'S',
                             'OCS':'C', 'OMT':'M', 'PAQ':'Y', 'PCA':'E', 'PEC':'C',
                             'PHI':'F', 'PHL':'F', 'PR3':'C', 'PRR':'A', 'PTR':'Y',
                             'SAC':'S', 'SAR':'G', 'SCH':'C', 'SCS':'C', 'SCY':'C',
                             'SEL':'S', 'SEP':'S', 'SET':'S', 'SHC':'C', 'SHR':'K',
                             'SOC':'C', 'STY':'Y', 'SVA':'S', 'TIH':'A', 'TPL':'W',
                             'TPO':'T', 'TPQ':'A', 'TRG':'K', 'TRO':'W', 'TYB':'Y',
                             'TYQ':'Y', 'TYS':'Y', 'TYY':'Y', 'AGM':'R', 'GL3':'G',
                             'SMC':'C', 'ASX':'B', 'CGU':'E', 'CSX':'C', 'GLX':'Z',
                             'UNK':'X'}

    def __init__(self, sequence, sequenceMD5=None, sequenceID=None, proteinMW=None, proteinIP=None):
        """
        "sequence": the sequence itself. It is processed...
        """

        self.proteinMW = proteinMW
        self.proteinIP = proteinIP

        self.biopython_protein_analyzer = None

        Sequence.__init__( self, sequence = sequence, sequenceMD5 = sequenceMD5, sequenceID = sequenceID, sequence_type="peptide" )

        


    def _get_protein_analyzer(self):

        if self.biopython_protein_analyzer is None:
            try:
                self.biopython_protein_analyzer = Bio.SeqUtils.ProtParam.ProteinAnalysis(self.get_sequence())
            except:
                pass
                #sys.stderr.write("ERROR! Bio.SeqUtils.ProtParam.ProteinAnalysis cannot be loaded.")

        return self.biopython_protein_analyzer
        

    def get_proteinMW(self):

        if self.proteinMW is None:
            try:
                # calculate the molecular weight
                self.proteinMW = self._get_protein_analyzer().molecular_weight()
            except:
                # any error in the function sets weight to 0
                # sys.stderr.write("ERROR! Assigning protein MW to 0.")
                self.proteinMW = 0
            
        return self.proteinMW


    def get_proteinIP(self):

        if self.proteinIP is None:
            try:
                self.proteinIP= analyzed_protein.isoelectric_point(correction_step = 0.001)
            except:
                #sys.stderr.write("ERROR! Assigning protein IP to 0.")
                self.proteinIP = 0

        return self.proteinIP



    def _get_digested_code(self, window):
        """
        """

        if ProteinSequence.translation_dict is None:
            autoinc = autoincrement(1)
            ProteinSequence.translation_dict = dict( [ ("%s%s%s" %(a,b,c),autoinc.next()) for a in ProteinSequence.dictionary
                                                       for b in ProteinSequence.dictionary
                                                       for c in ProteinSequence.dictionary ] )
                
        
        try:
            return ProteinSequence.translation_dict[window]
        except KeyError:
            return 0        

    def _get_fragmented_sequence(self):

        return Sequence._get_fragmented_sequence(self, 
                                                 using_translation_method = self._get_digested_code,
                                                 using_size = ProteinSequence.window_size,
                                                 using_list = ProteinSequence.dictionary )


    def get_aminoacid_code_3to1(code):
        
        try:
            return ProteinSequence.aminoacid_codes_3to1[code.upper()]
        except KeyError:
            print code
        
        return "X"

    get_aminoacid_code_3to1 = staticmethod(get_aminoacid_code_3to1)
