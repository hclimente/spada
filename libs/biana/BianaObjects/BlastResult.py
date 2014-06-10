import re
import sys

class BlastResult(object):
    """
    Object to store blast results

    All values are relative to absolute sequence position, starting at index 1!!!
    """

    def __init__(self, method, mode):

        self.sequenceID_A = None
        self.sequenceID_B = None
        self.e_value = None
        self.align_length = None
        self.score = None
        self.score_bits = None
        self.identities = None
        self.positives = None
        self.gaps = 0
        self.query_start = None
        self.query_end = None
        self.query_length = None
        self.sbjct_start = None
        self.sbjct_end = None
        self.sbjct_length = None
        self.method = method
        self.mode = mode
        self.query_coverage = None
        self.sbjct_coverage = None
        self.query_exact_match_list = []
        self.sbjct_exact_match_list = []
        self.query_similar_match_list = []
        self.sbjct_similar_match_list = []

    def __repr__(self):
    	return self.__str__()

    def __str__(self):

    	return "%s\n" %"\t".join(["%s" %self.sequenceID_A,
                                  "%s" %self.sequenceID_B,
                                  "%s" %self.e_value,
                                  "%s" %self.score,
                                  "%s" %self.score_bits,
                                  "%s" %self.query_start,
                                  "%s" %self.query_end,
                                  "%s" %self.sbjct_start,
                                  "%s" %self.sbjct_end,
                                  "%s" %self.identities,
                                  "%s" %self.positives,
                                  "%s" %self.gaps,
                                  "%s" %self.method,
                                  "%s" %self.mode,
                                  "%s" %self.get_query_coverage(),
			          "%s" %self.get_sbjct_coverage()])


    def get_query_coverage(self):
    	if self.query_coverage is None:
            if self.query_length is None:
                sys.stderr.write("Impossible to calculate query coverage: length not known. Query %s Subjct %s\n" %(self.sequenceID_A,self.sequenceID_B))
                sys.exit(1)
                return 0
            self.query_coverage = (int(self.query_end)-int(self.query_start))*100/float(self.query_length)
        return self.query_coverage

    def get_sbjct_coverage(self):
	if self.sbjct_coverage is None:
            if self.sbjct_length is None:
                sys.stderr.write("Impossible to calculate query coverage: length not known. Sbjct: %s\n"%self.sequenceID_B)
                return 0
            self.sbjct_coverage = (int(self.sbjct_end)-int(self.sbjct_start))*100/float(self.sbjct_length)
	return self.sbjct_coverage

    def set_evalue(self, e_value):
        if re.search("^e",e_value):
            self.e_value = float("1"+e_value)
        else:
            self.e_value = float(e_value)

    def write(self, fd_out):
     	fd_out.write(self.__str__())
    	
    def reset(self):

        self.e_value = None
        self.align_length = None
        self.score = None
        self.score_bits = None
        self.identities = None
        self.positives = None
        self.gaps = 0
        self.query_start = None
        self.query_end = None
        self.sbjct_start = None
        self.sbjct_end = None
        self.query_coverage = None
        self.sbjct_coverage = None

        
        self.query_exact_match_list = []
        self.sbjct_exact_match_list = []
        self.query_similar_match_list = []
        self.sbjct_similar_match_list = []
