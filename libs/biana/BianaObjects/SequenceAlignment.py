import sets
import sys

class SequenceAlignment(object):
    """
    Object to represent a sequence alignment, multiple or not
    """
    def __init__(self):

        # The alignment is represented as a list of sequence objects
        # Each sequence in the alignment consists on a list of aligmentBlocks, that define the aligment.
        # It defines the start position in the global alignment, the start position of the fragment and the length of the fragment
        # (global_start_position, start_position, length)
        self.sequence_ids_list = []
        self.cross_ids_list = []
        self.taxonomy_ids_list = [] # List with the associated taxonomys to aligned sequences
        self.sequence_fragments = []
        self.alignment = []
        self.alignment_length = 0
        self.nseqs = 0
        self.seq_ids_position_dict = {}
        self.complete_sequences = {}    # Dictionary to store the complete sequences of the alignment (needed to print the alignment)

        # Other parameters extracted from HSSP
        self.sim = []
        self.wsim = []
        self.lali = []
        self.ngap = []
        self.lgap = []

        self.C_align = None

    def _get_C_align(self):

        if( self.C_align is None ):

            import pprc

            self.C_align = pprc.C_align(total_proteins = self.nseqs,
                                        length = self.alignment_length)

            [ self.C_align.append(self.get_aligned_sequence(seq_pos=x)) for x in xrange(self.nseqs) ]
            
        return self.C_align
        

    def set_alignment_length(self, length):
        self.alignment_length=length

    def get_alignment_length(self):
        return self.alignment_length

    def get_number_of_sequences(self):
        return self.nseqs

    def set_number_of_sequences(self, number_of_sequences):
        self.nseqs = number_of_sequences
        self.sequence_ids_list = [ None for x in xrange(number_of_sequences) ]
        self.cross_ids_list = [ None for x in xrange(number_of_sequences) ]
        self.taxonomy_ids_list = [ None for x in xrange(number_of_sequences) ]
        self.alignment = [ None for x in xrange(number_of_sequences) ]
        self.sequence_fragments = [ None for x in xrange(number_of_sequences) ]

    def get_sequence_fragments_by_id(self, seq_id):
        return self.sequence_fragments[self.seq_ids_position_dict[str(seq_id)]]

    def _get_sequence_fragments_from_ascii(self, ascii_fragments):

        return [ ((ord(ascii_fragments[x])*255+ord(ascii_fragments[x+1])),
                 (ord(ascii_fragments[x+2])*255+ord(ascii_fragments[x+3])),
                 (ord(ascii_fragments[x+4])*255+ord(ascii_fragments[x+5]))) for x in xrange(0,len(ascii_fragments),6) ]

    def set_taxIDlist(self, sequenceID, taxIDslist):

        if len(self.taxonomy_ids_list)==0:
                self.taxonomy_ids_list = [ [] for x in xrange(number_of_sequences) ]

        try:
            self.taxonomy_ids_list[self.seq_ids_position_dict[sequenceID]] = taxIDslist
        except:
            pass


    def get_sequence_fragments_in_ascii(self,seq_position):
        
        temp_list = []
        for actual_fragment in self.sequence_fragments[seq_position]:
            temp_list.append(chr(actual_fragment[0]/255))
            temp_list.append(chr(actual_fragment[0]%255))
            temp_list.append(chr(actual_fragment[1]/255))
            temp_list.append(chr(actual_fragment[1]%255))
            temp_list.append(chr(actual_fragment[2]/255))
            temp_list.append(chr(actual_fragment[2]%255))
        return "".join(temp_list).replace('\\','\\\\').replace('"','\\"')

    def append_sequence(self,sequence_id,aligned_sequence,taxID_list=[],sim=None,wsim=None,lali=None,crossID=None):
        """
        Appends the aligned sequence to the alignment
        """

        if( not self.alignment_length ):
            self.alignment_length = len(aligned_sequence)
        elif len(aligned_sequence) != self.alignment_length:
            raise ValueError("Trying to add a protein in the alignment with distinct length (%s!=%s" %(self.alignment_length,len(aligned_sequence)) )

        self.seq_ids_position_dict[sequence_id] = self.nseqs
        self.nseqs += 1
        self.sequence_ids_list.append(sequence_id)
        self.cross_ids_list.append(crossID)
        self.taxonomy_ids_list.append(taxID_list)
        self.alignment.append(aligned_sequence)
        self.sequence_fragments.append(self._get_sequence_fragments(len(self.sequence_ids_list)-1))
        self.complete_sequences[sequence_id] = aligned_sequence.replace('-','').replace('.','')

        self.sim.append(sim)
        self.wsim.append(wsim)
        self.lali.append(lali)

    def append_sequence_fragments(self,sequence_id,sequence_fragments, complete_sequence,crossID=None):

        self.seq_ids_position_dict[sequence_id] = self.nseqs
        self.nseqs += 1
        self.sequence_ids_list.append(sequence_id)
        self.cross_ids_list.append(crossID)
        self.alignment.append(None)
        self.sequence_fragments.append(sequence_fragments)
        self.complete_sequences[sequence_id] = complete_sequence

    def _get_aligned_sequence(self, sequence_id):

        seq_pos = self.seq_ids_position_dict[sequence_id]

        if self.alignment[seq_pos] is None:
            seq = []
            position_start = 0
            for actual_fragment in self.sequence_fragments[seq_pos]:
                seq.extend(["-" for x in xrange(actual_fragment[0]-position_start)])
                position_start += actual_fragment[2]
                seq.append( str(self.complete_sequences[sequence_id])[actual_fragment[1]:actual_fragment[1]+actual_fragment[2]] )
            self.alignment[seq_pos] = "".join(seq)

        return self.alignment[seq_pos]


    def add_sequence_to_position(self,sequence_id,sequence_position,taxID_list=[],sequence_fragments_ascii=None,aligned_sequence=None,complete_sequence=None, crossID=None):

        if sequence_fragments_ascii is None and aligned_sequence is None:
            raise ValueError("At least one of the parameters must be specified")
        
        self.sequence_ids_list[sequence_position] = sequence_id
        self.cross_ids_list[sequence_position] = crossID
        self.taxonomy_ids_list[sequence_position] = taxID_list
        self.seq_ids_position_dict[sequence_id] = sequence_position

        if sequence_fragments_ascii:
            self.sequence_fragments[sequence_position] = self._get_sequence_fragments_from_ascii(ascii_fragments = sequence_fragments_ascii)
            self.complete_sequences[sequence_id] = complete_sequence
            #print self.sequence_fragments[sequence_position]
        if aligned_sequence:
            #print "adding as aligned sequence to position %s" %(sequence_position)
            self.alignment[sequence_position] = aligned_sequence
            self.complete_sequences[sequence_id] = aligned_sequence.replace('-','').replace('.','')
            #print aligned_sequence

    # This could be done in a C function (using PPRC...)
    def _get_sequence_fragments(self,aln_position):

        fragments = []
        in_gap = None
        counter = 0
        aln_start = 0
        seq_start = 0
        actual_pos = 0

        if self.alignment[aln_position][0] == "-":
            in_gap = 1

        for x in xrange(len(self.alignment[aln_position])):
            if in_gap:
                if self.alignment[aln_position][x]!="-":
                    in_gap = None
                    counter = 1
                    aln_start = x
                    seq_start = actual_pos
                    actual_pos += 1
                else:
                    continue
            else:
                if self.alignment[aln_position][x]!="-":
                    counter += 1
                    actual_pos += 1
                else:
                    fragments.append((aln_start,seq_start,counter))
                    #print "Start Position: %s\tLength: %s" %(aln_start,counter)
                    in_gap = 1

        if( self.alignment[aln_position][x]!="-" ):
            fragments.append((aln_start,seq_start,counter))
            #print "Start Position: %s\tLength: %s" %(aln_start,counter)

        #print self.alignment[aln_position]
        #print fragments
        return fragments

    def get_sequence_fragments(self):
        return self.sequence_fragments

    def get_sequence_fragments_by_pos(self,seq_pos):
        return self.sequence_fragments[seq_pos]

    #def append_sequence_fragment(self,sequence_id,global_start_pos,start_pos,end_pos):
    #    self.sequence_fragments

    def get_aligned_sequence(self,seq_pos,start_pos=None,end_pos=None):
        if self.alignment[seq_pos] is None:
            #print "Trying to get position %s, and maximum is %s" %(seq_pos,self.nseqs)
            self.alignment[seq_pos] = self._get_aligned_sequence(sequence_id = self.sequence_ids_list[seq_pos])
        if start_pos is None:
            return self.alignment[seq_pos]
        else:
            #print "Getting alignment from %s to %s" %(start_pos,end_pos+1)
            return self.alignment[seq_pos][start_pos:end_pos+1]

    def print_alignment(self, format="default",fd_out = sys.stdout):
        """
        Prints the alignment with the default identifier
        """
        for x in xrange(len(self.alignment)):
            #print "%s\t%s" %(self.sequence_ids_list[x],self.get_aligned_sequence(x))
            sequence = self.get_aligned_sequence(x)
            if self.cross_ids_list[x] is not None:
                id = self.cross_ids_list[x]#[1]
                #print "%s\t%s" %(self.cross_ids_list[x][1],self.get_aligned_sequence(x))
            else:
                id = "UNKNOWN IDENTIFIER"
                #print "%s\t%s" %("UNKNOWN IDENTIFIER",self.get_aligned_sequence(x))
            print id
            if format=="fasta":
                fd_out.write(">%s\n%s\n" %(id, sequence))
            else:
                fd_out.write("%s\t%s\n" %(id, sequence))
            #control:
            #print "%s\t%s [%s]" %(self.sequence_ids_list[x],self.get_aligned_sequence(x),self.sequence_fragments[x])


    def read_alignment(fd_in, format="default"):
        
        alignment = SequenceAlignment()
        
        if format=="default":
            for line in fd_in:
                line = line.strip()
                splitted = line.split("\t")
                alignment.append_sequence(sequence_id=splitted[0],
                                          aligned_sequence=splitted[1],
                                          crossID=splitted[0])
        else:
            raise ValueError("Not recognized alignment input format")

        return alignment

    read_alignment = staticmethod(read_alignment)
        
                    
    def get_alignment(self):
        """
        Returns the alignment in the following format:
        A list where each position contains the string corresponding to the sequence in string format
        """
        return self.alignment

    def concatenate_by_crossID(self, alignment2):

        print "concatenating %s + %s" %(self.alignment_length, alignment2.alignment_length)

        new_align = SequenceAlignment()

        crossID1_to_position = {}
        crossID2_to_position = {}

        for current_position in xrange(self.nseqs):
            crossID1_to_position.setdefault(self.cross_ids_list[current_position],current_position)
        for current_position in xrange(alignment2.nseqs):
            crossID2_to_position.setdefault(alignment2.cross_ids_list[current_position],current_position)
        
        for current_crossID in crossID1_to_position:
            if current_crossID in crossID2_to_position:
                new_align.append_sequence(sequence_id = self.sequence_ids_list[crossID1_to_position[current_crossID]],
                                          aligned_sequence = self.get_aligned_sequence(seq_pos = crossID1_to_position[current_crossID])+alignment2.get_aligned_sequence(seq_pos = crossID2_to_position[current_crossID]),
                                          crossID = current_crossID,
                                          taxID_list = self.taxonomy_ids_list[crossID1_to_position[current_crossID]] )
            
        return new_align

    def concatenate_alignments(self, alignment2):
     	"""
	Concatenate two alginments. They must have the  same number of sequences
	"""

	new_align = SequenceAlignment()

	
	for current_position in xrange(self.nseqs):
            id = "%s_%s" %(self.sequence_ids_list[current_position],alignment2.sequence_ids_list[current_position])
            new_align.append_sequence(sequence_id = id,
                                      aligned_sequence = self.get_aligned_sequence(seq_pos = current_position)+alignment2.get_aligned_sequence(seq_pos=current_position),
                                      crossID = id,
                                      taxID_list = self.taxonomy_ids_list[current_position] )
            
	return new_align


    def get_species_ordered_alignment(self, alignment2, method="only_first"):
        """
        Returns an alignment object where sequence positions have been ordered according to its specie
        
        "sequenceMD5_taxID_correspondence" is a dictionary with the correspondences between sequenceIDs and taxonomy

        method. Possibilities: "all_combinations", "only_first"
        """

        al1_taxIDs = {}

        al2_taxIDs = {}

        for x in xrange(len(self.taxonomy_ids_list)):
            for actual_taxID in self.taxonomy_ids_list[x]:
                if al1_taxIDs.has_key(actual_taxID):
                    al1_taxIDs[actual_taxID].add(x)
                else:
                    al1_taxIDs[actual_taxID] = sets.Set([x])

        for x in xrange(len(alignment2.taxonomy_ids_list)):
            for actual_taxID in alignment2.taxonomy_ids_list[x]:
                if al2_taxIDs.has_key(actual_taxID):
                    al2_taxIDs[actual_taxID].add(x)
                else:
                    al2_taxIDs[actual_taxID] = sets.Set([x])

        # Create the new alignment
        new_align1 = SequenceAlignment()
        new_align2 = SequenceAlignment()

        if method == "all_combinations":
            for actual_taxID in al1_taxIDs:
                if al2_taxIDs.has_key(actual_taxID):
                    #print actual_taxID," found in both alignments!"
                    for actual_position1 in al1_taxIDs[actual_taxID]:
                        for actual_position2 in al2_taxIDs[actual_taxID]:
                            new_align1.append_sequence(sequence_id = self.sequence_ids_list[actual_position1],
                                                       aligned_sequence= self.get_aligned_sequence(seq_pos=actual_position1),
                                                       crossID = self.cross_ids_list[actual_position1],
                                                       taxID_list = self.taxonomy_ids_list[actual_position1])
                            new_align2.append_sequence(sequence_id = alignment2.sequence_ids_list[actual_position2],
                                                       aligned_sequence = alignment2.get_aligned_sequence(seq_pos=actual_position2),
                                                       crossID = alignment2.cross_ids_list[actual_position2],
                                                       taxID_list = alignment2.taxonomy_ids_list[actual_position2])

        elif method == "only_first":
            for actual_taxID in al1_taxIDs:
                if al2_taxIDs.has_key(actual_taxID):
                    for actual_position1 in al1_taxIDs[actual_taxID]:
                        for actual_position2 in al2_taxIDs[actual_taxID]:
                            new_align1.append_sequence(sequence_id = self.sequence_ids_list[actual_position1],
                                                       aligned_sequence= self.get_aligned_sequence(seq_pos=actual_position1),
                                                       crossID = self.cross_ids_list[actual_position1],
                                                       taxID_list = self.taxonomy_ids_list[actual_position1])
                            new_align2.append_sequence(sequence_id = alignment2.sequence_ids_list[actual_position2],
                                                       aligned_sequence = alignment2.get_aligned_sequence(seq_pos=actual_position2),
                                                       crossID = alignment2.cross_ids_list[actual_position2],
                                                       taxID_list = alignment2.taxonomy_ids_list[actual_position2])
                            break
                        break
                    

        return (new_align1, new_align2)
                
                
                
    def get_subalignment(self, fragments):
        """
        "fragments" is a list of tuples with the format (start_position, end_position), with the fragments to select
        """

        new_align = SequenceAlignment()

        if len(fragments)==0:
            raise ValueError("Trying to get a subalignment without specifying fragments")

        for x in xrange(self.nseqs):
            aligned_seq = "".join([self.get_aligned_sequence(seq_pos=x,start_pos=actual_fragment[0],end_pos=actual_fragment[1])
                                   for actual_fragment in fragments])
            
            #if( aligned_seq.count('-') < len(aligned_seq) ):
            new_align.append_sequence(sequence_id = self.sequence_ids_list[x],
                                      aligned_sequence = aligned_seq,
                                      crossID = self.cross_ids_list[x],
                                      taxID_list = self.taxonomy_ids_list[x])
                
        return new_align


    def get_subalignment_for_sequence_without_gaps(self, sequenceID):
        """
        Returns a new alignment, in the same order, but eliminating all the positions in which sequenceID has a gap
        """

        fragments = []

        aligned_seq = self._get_aligned_sequence(sequence_id=sequenceID)

        current_start = None
        for x in xrange(len(aligned_seq)):
            if aligned_seq[x]!='-' and aligned_seq[x]!='.' and current_start is None:
                current_start = x
                continue
            if (aligned_seq[x]=='-' or aligned_seq[x]=='.') and current_start is not None:
                fragments.append((current_start, x-1))
                current_start = None
                continue

        if current_start is not None:
            fragments.append((current_start,len(aligned_seq)-1))

        return self.get_subalignment(fragments=fragments)

    
    def clean_empty_sequences( align1, align2 ):
        """
        
        """
        
        new_align1 = SequenceAlignment()
        new_align2 = SequenceAlignment()

        if( align1.nseqs != align2.nseqs ):
            raise ValueError("Both sequences must have the same number of proteins")

        for x in xrange(align1.nseqs):
            aligned_seq1 = align1.get_aligned_sequence(seq_pos=x)
            aligned_seq2 = align2.get_aligned_sequence(seq_pos=x)
            if( aligned_seq1.count('-') < len(aligned_seq1) and aligned_seq2.count('-') < len(aligned_seq2) ):
                new_align1.append_sequence( sequence_id = align1.sequence_ids_list[x],
                                            aligned_sequence = aligned_seq1 )
                new_align2.append_sequence( sequence_id = align2.sequence_ids_list[x],
                                            aligned_sequence = aligned_seq2 )

        return (new_align1, new_align2)
