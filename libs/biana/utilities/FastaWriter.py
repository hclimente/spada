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

from math import ceil

class FastaWriter(object):
    """
        FastaWriter class to handle outputting sequences in fasta format with customized header
    """

    # Sequence field in fasta file should not contain lines more than 80 characters (79 amino acids/nucleotides + '\n')
    symbol_per_line = 79

    def __init__(self, out_method, one_line_per_sequence=False):
        """
            FastaWriter constructer, creates an instance of FastaWriter
        """
        self.out_method = out_method
	self.one_line_per_seq = one_line_per_sequence
        return

    def _convert_sequence_to_fasta(self, sequence_header, sequence):
        """
            Converts given sequence into fasta format
        """
	sequence_body = ""
	if self.one_line_per_seq:
	    sequence_body = sequence
	else:
	    number_of_lines = int(ceil(len(sequence)/float(self.symbol_per_line)))
	    #print number_of_lines, sequence
	    sequence_body = "\n".join([ sequence[i*self.symbol_per_line:(i+1)*self.symbol_per_line] for i in xrange(number_of_lines) ])
	return ">%s\n%s\n" %(sequence_header, sequence_body)


    def output_sequence(self, sequence_header, sequence):
        """
            Outputs given sequence using out_method given in the initialization
        """
        self.out_method(self._convert_sequence_to_fasta(sequence_header, sequence))
        return


