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

import sys
import re
import os
import time
from sets import *
import gzip
import tempfile
import copy

from SequenceAlignment import *
import biana.biana_globals as biana_globals
from BlastResult import BlastResult


# TO CHECK
def get_clustalw_alignment(sequencesList):
    """
    Executes clustalw to get the multiple sequence alignment between the proteins in the list
    """

    temp_file = tempfile.NamedTemporaryFile(bufsize=0)
    [ temp_file.write(">%s\n%s\n" %(actual_sequence.get_sequenceID(), actual_sequence.get_sequence())) for actual_sequence in sequencesList ]

    # Execute clustalw
    command = "%s %s > /dev/null" %(biana_globals.CLUSTALW_EXEC, temp_file.name)
    #command = "%s %s" %(biana_globals.CLUSTALW_EXEC, temp_file.name)
    os.system(command)
    #print command
    
    temp_file.close()

    # Read results
    alignment = {}   # Ids are stored as keys
    sequence_ids_list = []
    
    result_re = re.compile("^([^S]+)\s+([\w\-]+)$")

    clustalw_output = open(temp_file.name+".aln",'r')

    seq_done = Set()

    for line in clustalw_output:
        m = result_re.match(line)
        if m:
            if m.group(1) not in seq_done:
                sequence_ids_list.append(m.group(1))
                seq_done.add(m.group(1))
            try:
                alignment[m.group(1)].append(m.group(2).strip())
            except:
                alignment[m.group(1)] = [m.group(2).strip()]

    clustalw_output.close()

    #delete files
    os.system("rm %s.*" %temp_file.name)
    temp_file.close()
    
    alignmentObject = SequenceAlignment()

    for actual_id in sequence_ids_list:
        alignmentObject.append_sequence( sequence_id = actual_id,
                                         aligned_sequence = "".join(alignment[actual_id]),
                                         crossID = actual_id )

    return alignmentObject


# TO CHECK
def get_cd_hit_clusters(fasta_file, output_path="./", sequence_identity_threshold = 0.95):
    """
    Executes CD-HIT with the sequences fasta file and saves the results in the ouput_path

    # TODO: Automate for executing in cluster (now it is done manually)
    """

    if( sequence_identity_threshold >= 0.7 ):
        n = 5
    elif( sequence_identity_threshold >= 0.6 ):
        n = 4
    elif( sequence_identity_threshold >= 0.5 ):
        n = 3
    elif( sequence_identity_threshold >= 0.4 ):
        n = 2

    if( sequence_identity_threshold > 1.0 or sequence_identity_threshold<0.4 ):
        raise ValueError("sequence_identity_threshold must be between 0.4 or 1.0")
    
    command = "%scd-hit -i %s -o %s -c %s -n %s -B 1 -p 1 -M 700" %(biana_globals.CD_HIT_PATH, fasta_file, output_path, sequence_identity_threshold, n)
    sys.stderr.write(command)

    cd_hit_result = os.system(command)
    
    if cd_hit_result:
        raise ValueError("CD-HIT has cannot be executed correctly")
    
        
    

# TO CHECK
def blast_cd_hit_clusters(cd_hit_clusters_file, output_fd, dbaccess, length_blast_db = None, effective_length_space_search = None):
    """
    Performs a blast between all the proteins belonging to the same cd-hit cluster
    """

    ## EXAMPLE OF Clusters output file
    ## >Cluster 0
    ## 0       36805aa, >9662380a987baf844... *
    ## >Cluster 1
    ## 0       35213aa, >0623d1531507bb004... *
    ## 1       90aa, >07b1a5a6aab138163... at 1:90:11834:11923/100%
    ## 2       247aa, >5e00433d0ae984091... at 1:247:12464:12710/100%
    ## 3       153aa, >d845d402ebfa7c203... at 1:153:11430:11582/100%
    ## 4       100aa, >ea8147fd16954563f... at 1:100:11483:11582/100%
    ## 5       183aa, >fab09a87d7f461e75... at 1:183:10973:11155/100%
    ## >Cluster 2
    ## 0       391aa, >1cacc168fca71118c... at 1:391:10887:11277/99%
    ## 1       34942aa, >d861c67fa2f88e4cd... *

    initial_time = time.time()

    input_file_fd = file(cd_hit_clusters_file, 'r')
    processed_clusters = 0
    new_cluster=0

    num_clusters=0

    fd_output_file = output_fd

    for line in input_file_fd:


        if re.search(">Cluster \d+",line):
            if new_cluster:
                processed_clusters += 1
                _calculate_similarities(dbaccess=dbaccess, sequenceID_list=cluster_sequences, fd_output_file=fd_output_file, length_blast_db = length_blast_db, effective_length_space_search = effective_length_space_search )

            new_cluster=1
            cluster_sequences = []
            num_clusters += 1
            if num_clusters%100==0:
                sys.stderr.write("%s clusters done in %seconds\n" %(num_clusters, time.time()-initial_time))

            # Control of the progress of the process
            #if processed_clusters >= number_clusters_percentage:
            #    # New percentage achieved
            #    actual_percentage += 1
            #    processed_clusters = 0
            #    if time_control:
            #        sys.stderr.write("%s per 10 mil clusters done in %s seconds\n" %(actual_percentage,time.time()-initial_time) )

        else:
            # Search for the sequenceID in the sequence line and append it to the cluster sequences list
            sequence = re.search("\d+\s+\d+aa,\s+\>(\w+)\.+",line)
            if sequence:
                cluster_sequences.append(int(sequence.group(1)))
            else:
                sys.stderr.write("%s" %line)


    if new_cluster:
        # Process the last line
        _calculate_similarities(sequenceID_list=cluster_sequences,dbaccess=dbaccess,fd_output_file=fd_output_file, length_blast_db = length_blast_db, effective_length_space_search = effective_length_space_search )





def blast_sequence(blastDatabase, sequenceObject, temporalOutputPath=None):
    """
    Does a blast of sequence sequenceObject against database "blastDatabase"

    ATTENTION: "temporalOutputPath" is the temporal gzipped file where the blast results are stored. If it exists, blast is not calculated!
    If it is None, results are not saved.
    """

    # First check if blast results are stored in file temporalOutputPath
    if os.path.exists(temporalOutputPath):
        if temporalOutputPath.endswith(".gz"):
            file_fd = gzip.open(temporalOutputPath,'r')
        else:
            file_fd = open(temporalOutputPath, 'r')
        blast_results = parse_blastall_output(file_fd)
        return blast_results
    
    command = biana_globals.BLASTALL_EXEC

    # TO CHECK: Do iterations???

    # Arguments for the command. set -F T if filter must be used, -F F otherwise
    # Argument -v 0 is to not print summary results in the output
    # -z is used to specify database size (in order to calculate e-value)
    args = ["-p","blastp","-d",blastDatabase,"-F","F","-v","0","-b","10000"]

    # Prepare pipes
    p_readfd, c_writefd = os.pipe()
    c_readfd, p_writefd = os.pipe()
            
    pid1 = os.fork()
    if pid1:
        # Parent
        for fd in (c_readfd, c_writefd, p_writefd):
            os.close(fd)
        # Convert the pipe fd to a file object, so we can use its
        # read() method to read all data.
        fp = os.fdopen(p_readfd, 'r')
        # result = fp.read()
        if temporalOutputPath.endswith(".gz"):
            blast_results_file = gzip.open(temporalOutputPath,'w')
        else:
            blast_results_file = open(temporalOutputPath,'w')
        blast_results = parse_blastall_output(fp, blast_results_file)
        
        fp.close()                      # Will close p_readfd.
        #sys.stderr.write(result)
        #os.wait()
    else:
        # Child
        try:
            pid = os.fork()
            if pid:
                # Still the same child
                # for sequence in sequenceObjectDict.values():
                # os.write(p_writefd, proteinSequence+"\n\n")
                # os.write(p_writefd, "%s\n\n" %sequence.get_
                os.write(p_writefd, ">%s\n%s\n\n" %(sequenceObject.get_sequenceID(), sequenceObject.get_sequence()))
                os.close(p_writefd)
            else:
                # Grandchild
                try:
                    # Redirect the pipe to stdin.
                    os.close(0)
                    os.dup(c_readfd)
                    # Redirect stdout to the pipe.
                    os.close(1)
                    os.dup(c_writefd)
                    # Now close unneeded descriptors.
                    for fd in (c_readfd, c_writefd, p_readfd, p_writefd):
                        os.close(fd)
                    # Finally, execute the external command.
                    os.execv(command, [command] + args)
                except:
                    safe_traceback()
                    os._exit(127)
        except:
            safe_traceback()
            os._exit(127)
        else:
            os._exit(0)

    return blast_results


def _calculate_similarities(dbaccess, length_blast_db = None, effective_length_space_search = None, sequenceID_list=[], fd_output_file=sys.stdout, representant=None):
    """
    Calculates the similarity between all proteins in the list "cluster_sequences".

    "cluster_sequences" must be a list of the sequenceIDs

    It uses bl2seq or blastall to calculate them

    Results are printed to "fd_output_file"

    """

    if len(sequenceID_list)<=1:
        # It is not necessary to calculate anything, as the cluster has an unique sequence
        return
    elif len(sequenceID_list)>10:
        # If the number of sequences in the cluster is larger than this cutoff, use blastall

        def safe_traceback():
            # Child processes catch exceptions so that they can exit using
            # os._exit() without fanfare.  They use this function to print
            # the traceback to stderr before dying.
            import traceback
            sys.stderr.write("Error in child process, pid %d.\n" %os.getpid())
            sys.stderr.flush()
            traceback.print_exc()
            sys.stderr.flush()
        
        # Make a database with this cluster, and then run blastall with the sequence against this database
        file_prefix = "./cluster_%s" %sequenceID_list[0]
        file_name = file_prefix+".fa"
        out_fasta_file = open(file_name, 'w')
        
        # Obtain the sequences for all the cluster
        # proteinSequences = [piana_access.get_protein_sequence(proteinPiana=x) for x in cluster_sequences]
        #proteinSequences = [piana_access.get_sequence_from_sequenceID(sequenceID = x) for x in cluster_sequences]
        load_time = time.time()
        sequenceObjectDict = dbaccess._load_sequences( sequenceIdList = sequenceID_list, type = "proteinsequence" )
        #print "%s sequences loaded in %s seconds" %(len(sequenceID_list), time.time()-load_time)
        #print_sequences_in_fasta_format( outmethod = out_fasta_file, sequenceObjList = sequenceObjectDict.values() )
        [ out_fasta_file.write(">%s\n%s\n" %(actual_sequence.get_sequenceID(), actual_sequence.get_sequence())) for actual_sequence in sequenceObjectDict.itervalues() ]
        out_fasta_file.close()

        # Generate the fasta file
        #output_file_fd = file(file_name, "w")

        #for i_sequenceID in xrange(len(cluster_sequences)):
        #    output_file_fd.write(">%s\n%s\n" %(cluster_sequences[i_sequenceID],proteinSequences[i_sequenceID]) )

        #output_file_fd.close()

        # Format the database

        command = "%s -t %s -i %s -l %s.log -p T -a F -o F" %(biana_globals.FORMATDB_EXEC, file_name, file_name, file_prefix)

        format_db = os.system(command)

        # Execute blastall with all sequences against the database
        
        #command = "blastall"
        command = biana_globals.BLASTALL_EXEC

        # Arguments for the command. set -F T if filter must be used, -F F otherwise
        # Argument -v 0 is to not print summary results in the output
        # -z is used to specify database size (in order to calculate e-value)
        args = ["-p","blastp","-d",file_name,"-z","1720800858","-F","F","-v","0","-b","1000000"]
        if length_blast_db is not None:
            args.append("-z")
            args.append(length_blast_db)
        if effective_length_space_search is not None:
            args.append("-Y")
            args.append(effective_length_space_search)

        args = map(str,args)

        # Prepare pipes
        p_readfd, c_writefd = os.pipe()
        c_readfd, p_writefd = os.pipe()
            
        pid1 = os.fork()
        if pid1:
            # Parent
            for fd in (c_readfd, c_writefd, p_writefd):
                os.close(fd)
            # Convert the pipe fd to a file object, so we can use its
            # read() method to read all data.
            fp = os.fdopen(p_readfd, 'r')
            #result = fp.read()
            blastall_iterator = BlastallParserIterator( fd_blastall_output = fp, parse_detailed_alignments = False )
            #blast_results = parse_blastall_output(fp)
            #[ fd_output_file.write(x.__str__()) for x in blast_results ]

            try:
                while True:
                    fd_output_file.write(blastall_iterator.next().__str__())
            except:
                pass
            
            fp.close()                      # Will close p_readfd.
            #sys.stderr.write(result)
            os.wait()   #UNCOMMENTED TO TEST... JAVI


            os.remove(file_prefix+".fa")
            os.remove(file_prefix+".log")
            os.remove(file_prefix+".fa.phr")
            os.remove(file_prefix+".fa.pin")
            os.remove(file_prefix+".fa.psq")
        else:
            # Child
            try:
                pid = os.fork()
                if pid:
                    # Still the same child
                    #for sequence in sequenceObjectDict.values():
                        #os.write(p_writefd, proteinSequence+"\n\n")
                    #    os.write(p_writefd, "%s\n\n" %sequence.get_sequence())
                    if representant is None:
                        [ os.write(p_writefd, ">%s\n%s\n\n" %(x.sequenceID,x.get_sequence())) for x in sequenceObjectDict.values() ]
                    else:
                        raise NotImplementedError("To implement representat similarity searches")
                        #os.write(p_writefd, "%s\n\n" %representant)
                    #os.wait()  # UNCOMMENTED TO TEST
                    os._exit(0) # ADDED TO TEST
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
                        os.execv(command, [command] + args)
                        
                    except:
                        safe_traceback()
                        os._exit(127)
            except:
                safe_traceback()
                os._exit(127)
            else:
                os._exit(0)

    else:
        # Run bl2seq for all these sequences combinations

        # Create the temporary files
        temp_file_name = "./temp_seqfile_"

        temporal_files = []

        sequenceObjectDict = dbaccess._load_sequences( sequenceIdList = sequenceID_list, type = "proteinsequence" )

        # Obtain and save the sequences for the cluster
        #for sequenceID in cluster_sequences:
        #for sequenceID in sequenceObjectDict:
        for sequenceID in sequenceID_list:
            #sequence = piana_access.get_sequence_from_sequenceID(sequenceID = sequenceID)
            # Save it into a temporary file
            file_name = temp_file_name+str(sequenceID)
            temporal_files.append(file_name)
            output_file_fd = file(file_name, "w")
            output_file_fd.write(">%s\n%s\n" %(sequenceID,sequenceObjectDict[sequenceID].get_sequence()) )
            output_file_fd.close()
        
        for i in xrange(len(sequenceID_list)):
            for j in xrange(i+1,len(sequenceID_list)):
                # execute bl2seq for this two sequences
                #command = "%s -i %s -j %s -p blastp -d 1720800858 -F F" %(bl2seq, file_name+sequenceID_list[i],file_name+sequenceID_list[j])
                #command = "%s -i %s -j %s -p blastp -d 1720800858 -F F" %(biana_globals.BL2SEQ_EXEC, temporal_files[i],temporal_files[j])
                
                command_list = [biana_globals.BL2SEQ_EXEC,"-i",temporal_files[i],"-j",temporal_files[j],"-p","blastp","-F","F"]

                if length_blast_db is not None:
                    command_list.append("-d")
                    command_list.append(length_blast_db)
                if effective_length_space_search is not None:
                    command_list.append("-Y")
                    command_list.append(effective_length_space_search)

                command = " ".join(map(str,command_list))


                # TEST OF LOOP TO AVOID "Resource temporarily unavailable" EXCEPTION
                # uncommented because this was really ugly..
                #while(True):
                    #try:
                bl2seq_f = os.popen(command,"r")
                bl2seq_data = bl2seq_f.read()
                bl2seq_f.close()
                #break
                    #except:
                    #    time.sleep(1)

                #sys.stderr.write(bl2seq_data)
                
                parse_bl2seq_output(sequenceID_A = sequenceID_list[i], sequenceID_B = sequenceID_list[j], bl2seq_output=bl2seq_data,fd_output_file=fd_output_file)

        # Deleting temporary files
        map(os.remove, temporal_files)


def self_blast_fasta_file(fasta_file, fd_output_file):
    """
    Does a blast within all sequences in a file

    If file is very big, it won't work. It is intended for small sets
    """

    def safe_traceback():
        # Child processes catch exceptions so that they can exit using
        # os._exit() without fanfare.  They use this function to print
        # the traceback to stderr before dying.
        import traceback
        sys.stderr.write("Error in child process, pid %d.\n" %os.getpid())
        sys.stderr.flush()
        traceback.print_exc()
        sys.stderr.flush()

    # Read all the sequences in the file
    sequences_list = []
    ids_list = []

    input_file_fd = file(fasta_file, 'r')

    for line in input_file_fd:
        if line[0] == ">":
            ids_list.append(line[0:].strip())
        else:
            sequences_list.append(line.strip())

    input_file_fd.close()

    # Format the database
    #import tempfile

    #fasta_file = tempfile.NamedTemporaryFile(bufsize=0,prefix="biana_fasta")
    #log_file = tempfile.NamedTemporaryFile(bufsize=0,prefix="biana_log_fasta")

    command = "%s -t clusterFASTA -i %s -l clusterFASTA.log -p T -a F -o F" %(biana_globals.FORMATDB_EXEC, fasta_file)
    
    format_db = os.system(command)   
    
    command = biana_globals.BLASTALL_EXEC
    args = ["-p","blastp","-d",fasta_file,"-z","1720800858","-F","F","-v","0"]


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
        #result = fp.read()
        blast_results = parse_blastall_output(fp)
        [ fd_output_file.write(x.__str__()) for x in blast_results ]
        fp.close()                      # Will close p_readfd.
        #sys.stderr.write(result)
    else:
        # Child
        try:
            if os.fork():
                # Still the same child
                for current_sequence in sequences_list:
                    os.write(p_writefd, "%s\n\n" %current_sequence)
            else:
                # Grandchild
                try:
                    # Redirect the pipe to stdin.
                    os.close(0)
                    os.dup(c_readfd)
                    # Redirect stdout to the pipe.
                    os.close(1)
                    os.dup(c_writefd)
                    # Now close unneeded descriptors.
                    for fd in (c_readfd, c_writefd, p_readfd, p_writefd):
                        os.close(fd)
                    # Finally, execute the external command.
                    os.execv(command, [command] + args)
                except:
                    safe_traceback()
                    os._exit(127)
        except:
            safe_traceback()
            os._exit(127)
        else:
            os._exit(0)


class BlastallParserIterator(object):
    """
    """

    query_re = re.compile("Query=\s*(.+)\s*")
    letters_re = re.compile("\(\s*([\,\d]+)\s*letters\s*")
    sbjct_re = re.compile(">([\w\d\_\.\|]+)")
    length_re = re.compile("Length \= ([\,\d]+)")
    score_re = re.compile("Score\s+=\s+([\.\d]+)\s+bits\s+\((\d+)\),\s+Expect\s+=\s+([\d\.e\-]+)")
    score_option2_re = re.compile("Score\s+=\s+([\.\d]+)\s+\(([\.\d]+)\s+bits\)\,\s+Expect\s+=\s+([\d\.e\-]+)")
    identities_re = re.compile("Identities\s+=\s+\d+\/(\d+)\s+\((\d+)%\)")
    positives_re = re.compile("Positives\s+=\s+\d+\/\d+\s+\((\d+)%\)")
    gaps_re = re.compile("Gaps\s+=\s+\d+\/\d+\s+\((\d+)%\)")
    intervals_query_re = re.compile("Query:\s+(\d+)\s+(\S+)\s+(\d+)$")
    sbjct_intervals_re = re.compile("Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)$")

    def __init__(self, fd_blastall_output, parse_detailed_alignments=False):

        self.fd = fd_blastall_output
        
        #self.current_sequenceID_A = None
        #self.current_query_length = None
        #self.current_sequenceID_B = None
        #self.current_sbjct_length = None
        self.current_blastResult_obj = None
        self.parse_lines = parse_detailed_alignments

    def __iter__(self):
        return self


    def parse_alignment_line(self, alignment_line_list, aligned_query, aligned_sbjct):
        """
        """
      
 
        # BUG AMB ELS GAPS!!!!!!!!!!
        print "Parsing line"
        print alignment_line_list
        alignment_line = "".join(alignment_line_list)
        aligned_query = "".join(aligned_query)
        aligned_sbjct = "".join(aligned_sbjct)

        if len(alignment_line) != len(aligned_query) or len(aligned_query) != len(aligned_sbjct):
            print aligned_query
            print alignment_line
            print aligned_sbjct
            raise ValueError("Alignments must be of the same size")

        query_gaps = 0
        sbjct_gaps = 0

        for x in xrange(len(alignment_line)):
            value = alignment_line[x]
            if aligned_query[x]=="-":
                query_gaps += 1
            if aligned_sbjct[x]=="-":
                sbjct_gaps += 1
            if value == " ":
                continue
            else:
                if value != "+":
                    self.current_blastResult_obj.query_exact_match_list.append(x+self.current_blastResult_obj.query_start-query_gaps)
                    self.current_blastResult_obj.sbjct_exact_match_list.append(x+self.current_blastResult_obj.sbjct_start-sbjct_gaps)
                    self.current_blastResult_obj.query_similar_match_list.append(x+self.current_blastResult_obj.query_start-query_gaps)
                    self.current_blastResult_obj.sbjct_similar_match_list.append(x+self.current_blastResult_obj.sbjct_start-sbjct_gaps)
                else:
                    self.current_blastResult_obj.query_similar_match_list.append(x+self.current_blastResult_obj.query_start-query_gaps)
                    self.current_blastResult_obj.sbjct_similar_match_list.append(x+self.current_blastResult_obj.sbjct_start-sbjct_gaps)

    def next(self):

        # Temporal variables to store information to read exact alignment
        alignment_start_index = None
        capture_matching_line = False
        sbjct_matching = False
        alignment_summary = []
        aligned_query = []
        aligned_sbjct = []
        
        #blastResult_obj = BlastResult(method="blastall",mode="F")
        #blastResult_obj.sequenceID_A = self.current_sequenceID_A
        #blastResult_obj.query_length = self.current_query_length
        #blastResult_obj.sequenceID_B = self.current_sequenceID_B
        #blastResult_obj.sbjct_length = self.current_sbjct_length

        for line in self.fd:
            #print line

            if capture_matching_line:
                alignment_summary.append(line[alignment_start_index:alignment_start_index+subalignment_length])
                capture_matching_line = False
                continue
                                                 
            m = BlastallParserIterator.query_re.search(line)

            if m:
                self.current_blastResult_obj = BlastResult(method="blastall",mode="F")
                self.current_blastResult_obj.sequenceID_A = m.group(1)

            m = BlastallParserIterator.letters_re.search(line)
            if m:
                self.current_blastResult_obj.query_length = int(m.group(1).replace(",",''))

            sequenceID_B_search = BlastallParserIterator.sbjct_re.search(line)

            if sequenceID_B_search:
                #print line
                # New sequenceID_B:
                if self.current_blastResult_obj.e_value is not None:
                    if self.current_blastResult_obj.sequenceID_A != self.current_blastResult_obj.sequenceID_B:
                        if self.parse_lines:
                            self.parse_alignment_line(alignment_summary, aligned_query, aligned_sbjct)
                        t = self.current_blastResult_obj
                        self.current_blastResult_obj = BlastResult(method="blastall",mode="F")
                        self.current_blastResult_obj.sequenceID_A = t.sequenceID_A
                        self.current_blastResult_obj.sequenceID_B = sequenceID_B_search.group(1)
                        self.current_blastResult_obj.query_length = t.query_length
                        return t
                else:
                    self.current_blastResult_obj.sequenceID_B = sequenceID_B_search.group(1)
                                
            m = BlastallParserIterator.length_re.search(line)
            if m:
                self.current_blastResult_obj.sbjct_length = int(m.group(1).replace(",",''))

            if re.search("^Matrix",line):
                # Query finished
                if self.current_blastResult_obj.e_value is not None:
                    if self.current_blastResult_obj.sequenceID_A != self.current_blastResult_obj.sequenceID_B:
                        if self.parse_lines:
                            self.parse_alignment_line(alignment_summary, aligned_query, aligned_sbjct)
                        t = self.current_blastResult_obj
                        self.current_blastResult_obj = BlastResult(method="blastall",mode="F")
                        return t

            else:
                get_evalue = BlastallParserIterator.score_re.search(line)
                if not get_evalue:
                	get_evalue = BlastallParserIterator.score_option2_re.search(line)
                if get_evalue:
                    # New hit found
                    if self.current_blastResult_obj.e_value is not None:    ## Check if there were other hits before
                        if self.current_blastResult_obj.sequenceID_A != self.current_blastResult_obj.sequenceID_B:
                            if self.parse_lines:
                                self.parse_alignment_line(alignment_summary, aligned_query, aligned_sbjct)
                            t = self.current_blastResult_obj
                            self.current_blastResult_obj = BlastResult(method="blastall",mode="F")
                            self.current_blastResult_obj.sequenceID_A = t.sequenceID_A
                            self.current_blastResult_obj.query_length = t.query_length
                            self.current_blastResult_obj.set_evalue(get_evalue.group(3))
                            self.current_blastResult_obj.score_bits = str(get_evalue.group(1))
                            self.current_blastResult_obj.score = str(get_evalue.group(2))
                            self.current_blastResult_obj.sequenceID_B = t.sequenceID_B
                            self.current_blastResult_obj.sbjct_length = t.sbjct_length
                            return t

                    self.current_blastResult_obj.set_evalue(get_evalue.group(3))
                    self.current_blastResult_obj.score_bits = str(get_evalue.group(1))
                    self.current_blastResult_obj.score = str(get_evalue.group(2))

                get_identities = BlastallParserIterator.identities_re.search(line)

                if get_identities:
                    self.current_blastResult_obj.align_length = get_identities.group(1)
                    self.current_blastResult_obj.identities= str(get_identities.group(2))

                get_positives = BlastallParserIterator.positives_re.search(line)

                if get_positives:
                    self.current_blastResult_obj.positives = str(get_positives.group(1))

                get_gaps = BlastallParserIterator.gaps_re.search(line)
                if get_gaps:
                    self.current_blastResult_obj.gaps = str(get_gaps.group(1))

                get_intervals_query = BlastallParserIterator.intervals_query_re.search(line)

                if get_intervals_query:
                    if self.current_blastResult_obj.query_start is None:
                        self.current_blastResult_obj.query_start = int(get_intervals_query.group(1))
                        alignment_start_index = line.index(get_intervals_query.group(2))
                    subalignment_length = len(get_intervals_query.group(2))
                    capture_matching_line = True
                    sbjct_matching = True
                    aligned_query.append(get_intervals_query.group(2))
                    self.current_blastResult_obj.query_end = int(get_intervals_query.group(3))

                get_intervals_Sbjct = BlastallParserIterator.sbjct_intervals_re.search(line)

                if get_intervals_Sbjct and sbjct_matching:
                    if self.current_blastResult_obj.sbjct_start is None:
                        self.current_blastResult_obj.sbjct_start = int(get_intervals_Sbjct.group(1))
                    self.current_blastResult_obj.sbjct_end = int(get_intervals_Sbjct.group(3))
                    aligned_sbjct.append(get_intervals_Sbjct.group(2))
                    sbjct_matching = False

        raise StopIteration


# CHECKED
def parse_blastall_output(fd_blastall_output, temporalOutputFile_fd=None, return_only_ids=False, limit_to_sequenceIDs=sets.Set()):
    """
    "fd_blastall_output" is the output fd of the blast process (input for this method)
    
    "temporalOutputFile" is a file where all the input of fd_blastall_output is saved

    "return_only_ids" is used to store only ids, not complete blast results

    "limit_to_sequenceIDs" is used to filter blast parsing to only those sequenceids
    """

    blast_results = []

    query_re = re.compile("Query=\s*(.+)\s*")    # IT WAS INCORRECT.... DID IT AFFECT ANY RESULT???
    letters_re = re.compile("\(\s*([\,\d]+)\s*letters\s*\)")
    sbjct_re = re.compile("^>([\w\d\_\.\|]+)")
    length_re = re.compile("Length \= ([\,\d]+)")
    score_re = re.compile("Score\s+=\s+([\.\d]+)\s+bits\s+\((\d+)\),\s+Expect\s+=\s+([\d\.e\-]+)")
    identities_re = re.compile("Identities\s+=\s+\d+\/(\d+)\s+\((\d+)%\)")
    positives_re = re.compile("Positives\s+=\s+\d+\/\d+\s+\((\d+)%\)")
    gaps_re = re.compile("Gaps\s+=\s+\d+\/\d+\s+\((\d+)%\)")
    intervals_query_re = re.compile("Query:\s+(\d+)\s+(\S+)\s+(\d+)$")
    sbjct_intervals_re = re.compile("Sbjct:\s+(\d+)\s+(\S+)\s+(\d+)$")

    # Temporal variables to store information to read exact alignment
    alignment_start_index = None
    capture_matching_line = False
    sbjct_matching = False
    alignment_summary = []
    
    blastResult_obj = BlastResult(method="blastall",mode="F")

    def parse_alignment_line(alignment_line_list, blastResult_obj, aligned_query, aligned_sbjct):
        """
        
        """
        
        # BUG AMB ELS GAPS!!!!!!!!!!

        alignment_line = "".join(alignment_line_list)
        aligned_query = "".join(aligned_query)
        aligned_sbjct = "".join(aligned_sbjct)

        if len(alignment_line) != len(aligned_query) or len(aligned_query) != len(aligned_sbjct):
            print aligned_query
            print alignment_line
            print aligned_sbjct
            raise ValueError("Alignments must be of the same size")

        query_gaps = 0
        sbjct_gaps = 0

        for x in xrange(len(alignment_line)):
            value = alignment_line[x]
            if aligned_query[x]=="-":
                query_gaps += 1
            if aligned_sbjct[x]=="-":
                sbjct_gaps += 1
            if value == " ":
                continue
            else:
                if value != "+":
                    blastResult_obj.query_exact_match_list.append(x+blastResult_obj.query_start-query_gaps)
                    blastResult_obj.sbjct_exact_match_list.append(x+blastResult_obj.sbjct_start-sbjct_gaps)
                    blastResult_obj.query_similar_match_list.append(x+blastResult_obj.query_start-query_gaps)
                    blastResult_obj.sbjct_similar_match_list.append(x+blastResult_obj.sbjct_start-sbjct_gaps)
                else:
                    blastResult_obj.query_similar_match_list.append(x+blastResult_obj.query_start-query_gaps)
                    blastResult_obj.sbjct_similar_match_list.append(x+blastResult_obj.sbjct_start-sbjct_gaps)

        #print blastResult_obj.query_similar_match_list
        #print blastResult_obj.query_exact_match_list
        #print blastResult_obj.sbjct_exact_match_list

    for line in fd_blastall_output:

        if temporalOutputFile_fd:
            temporalOutputFile_fd.write(line)

        if capture_matching_line:
            alignment_summary.append(line[alignment_start_index:alignment_start_index+subalignment_length])
            capture_matching_line = False
            continue
                                                 
        m = query_re.search(line)
        if m:
            sequenceID_A = m.group(1)
            blastResult_obj.sequenceID_A = sequenceID_A

        m = letters_re.search(line)
        if m:
            blastResult_obj.query_length = int(m.group(1).replace(",",''))

        sequenceID_B_search = sbjct_re.search(line)
        if sequenceID_B_search:

            # New sequenceID_B:
            if blastResult_obj.e_value is not None:
                if blastResult_obj.sequenceID_A != blastResult_obj.sequenceID_B:
                    parse_alignment_line(alignment_summary, blastResult_obj, aligned_query, aligned_sbjct)
                    if return_only_ids:
                        blast_results.append(blastResult_obj.sequenceID_B)
                    else:
                        if len(limit_to_sequenceIDs)==0 or blastResult_obj.sequenceID_B in limit_to_sequenceIDs:
                            blast_results.append(copy.copy(blastResult_obj))

            alignment_start_index = None
            alignment_summary = []
            aligned_query = []
            aligned_sbjct = []

            blastResult_obj.reset()
            blastResult_obj.sequenceID_B = sequenceID_B_search.group(1)

        m = length_re.search(line)
        if m:
            blastResult_obj.sbjct_length = int(m.group(1).replace(",",''))

        if re.search("^Matrix",line):
            # Query finished
            if blastResult_obj.e_value is not None:
                if blastResult_obj.sequenceID_A != blastResult_obj.sequenceID_B:
                    parse_alignment_line(alignment_summary, blastResult_obj, aligned_query, aligned_sbjct)
                    if return_only_ids:
                        blast_results.append(blastResult_obj.sequenceID_B)
                    else:
                        if len(limit_to_sequenceIDs)==0 or blastResult_obj.sequenceID_B in limit_to_sequenceIDs:
                            blast_results.append(blastResult_obj)

            alignment_start_index = None
            alignment_summary = []
            aligned_query = []
            aligned_sbjct = []
            
            blastResult_obj = BlastResult(method="blastall",mode="F")

        else:
            get_evalue = score_re.search(line)
            if get_evalue:
                # New hit found
                if blastResult_obj.e_value is not None:    ## Check if there were other hits before
                    if blastResult_obj.sequenceID_A != blastResult_obj.sequenceID_B:
                        parse_alignment_line(alignment_summary, blastResult_obj, aligned_query, aligned_sbjct)
                        if return_only_ids:
                            blast_results.append(blastResult_obj.sequenceID_B)
                        else:
                            if len(limit_to_sequenceIDs)==0 or blastResult_obj.sequenceID_B in limit_to_sequenceIDs:
                                blast_results.append(copy.copy(blastResult_obj))

                blastResult_obj.reset()

                alignment_start_index = None
                alignment_summary = []
                aligned_query = []
                aligned_sbjct = []

                blastResult_obj.set_evalue(get_evalue.group(3))
                blastResult_obj.score_bits = str(get_evalue.group(1))
                blastResult_obj.score = str(get_evalue.group(2))

            get_identities = identities_re.search(line)

            if get_identities:
                blastResult_obj.align_length = get_identities.group(1)
                blastResult_obj.identities= str(get_identities.group(2))

            get_positives = positives_re.search(line)

            if get_positives:
                blastResult_obj.positives = str(get_positives.group(1))

            get_gaps = gaps_re.search(line)
            if get_gaps:
                blastResult_obj.gaps = str(get_gaps.group(1))

            get_intervals_query = intervals_query_re.search(line)

            if get_intervals_query:
                if blastResult_obj.query_start is None:
                    blastResult_obj.query_start = int(get_intervals_query.group(1))
                    alignment_start_index = line.index(get_intervals_query.group(2))
                subalignment_length = len(get_intervals_query.group(2))
                capture_matching_line = True
                sbjct_matching = True
                aligned_query.append(get_intervals_query.group(2))
                blastResult_obj.query_end = int(get_intervals_query.group(3))

            get_intervals_Sbjct = sbjct_intervals_re.search(line)

            if get_intervals_Sbjct and sbjct_matching:
                if blastResult_obj.sbjct_start is None:
                    blastResult_obj.sbjct_start = int(get_intervals_Sbjct.group(1))
                blastResult_obj.sbjct_end = int(get_intervals_Sbjct.group(3))
                aligned_sbjct.append(get_intervals_Sbjct.group(2))
                sbjct_matching = False

    return blast_results


def parse_bl2seq_output(sequenceID_A, sequenceID_B, bl2seq_output=None, fd_output_file=None):

    score_re = re.compile("Score\s+=\s+([\.\d]+)\s+bits\s+\((\d+)\),\s+Expect\s+=\s+([\d\.e\-]+)")
    identities_re = re.compile("Identities\s+=\s+\d+\/(\d+)\s+\((\d+)%\)")
    positives_re = re.compile("Positives\s+=\s+\d+\/\d+\s+\((\d+)%\)")
    gaps_re = re.compile("Gaps\s+=\s+\d+\/\d+\s+\((\d+)%\)")
    intervals_query_re = re.compile("Query:\s+(\d+)\s+.+\s(\d+)$")
    sbjct_intervals_re = re.compile("Sbjct:\s+(\d+)\s+.+\s(\d+)$")
    #intervals_query_re = re.compile("Query:\s+(\d+)\s+\S+\s+(\d+)$")
    #sbjct_intervals_re = re.compile("Sbjct:\s+(\d+)\s+\S+\s+(\d+)$")
    letters_re = re.compile("\(\s*([\d\,]+)\s*letters\s*\)")
    length_re = re.compile("Length\s+\=\s+([\d\,]+)")

    if fd_output_file is None:
        fd_output_file = sys.stdout

    if bl2seq_output is None:
        return
    else:

        # Split the output in lines
        bl2seq_lines = bl2seq_output.split("\n")

        blastResult_obj = BlastResult(method="bl2seq",mode="F")
        blastResult_obj.sequenceID_A = sequenceID_A
        blastResult_obj.sequenceID_B = sequenceID_B

        for line in bl2seq_lines:
            if re.search("Lambda",line):
                # Useful information is finished
                # Appending the last result
                if blastResult_obj.e_value is not None:
                    if blastResult_obj.e_value < 0.1:
                        fd_output_file.write(str(blastResult_obj))
                        # blastResult_obj.write(fd_output_file)
                        blastResult_obj.reset()
                    
            else:
                get_evalue = score_re.search(line)
                if get_evalue:
                    # New hit found
                    if blastResult_obj.e_value is not None:
                        if blastResult_obj.e_value < 0.1:
                            #blastResult_obj.write(fd_output_file)
                            fd_output_file.write(str(blastResult_obj))
                            blastResult_obj.reset()

                    blastResult_obj.set_evalue(get_evalue.group(3))
                    blastResult_obj.score= get_evalue.group(2)
                    blastResult_obj.score_bits = get_evalue.group(1)

                m = letters_re.search(line)
                if m:
                    blastResult_obj.query_length = int(m.group(1).replace(',',''))

                m = length_re.search(line)
                if m:
                    blastResult_obj.sbjct_length = int(m.group(1).replace(',',''))

                get_identities = identities_re.search(line)
                if get_identities:
                    blastResult_obj.align_length = get_identities.group(1)
                    blastResult_obj.identities= get_identities.group(2)

                get_positives = positives_re.search(line)
                if get_positives:
                    blastResult_obj.positives = get_positives.group(1)

                get_gaps = gaps_re.search(line)
                if get_gaps:
                    blastResult_obj.gaps = get_gaps.group(1)

                get_intervals_query = intervals_query_re.search(line)
                if get_intervals_query:
                    if blastResult_obj.query_start is None:
                        blastResult_obj.query_start = get_intervals_query.group(1)
                    blastResult_obj.query_end = get_intervals_query.group(2)

                get_intervals_Sbjct = sbjct_intervals_re.search(line)
                if get_intervals_Sbjct:
                    if blastResult_obj.sbjct_start is None:
                        blastResult_obj.sbjct_start = get_intervals_Sbjct.group(1)
                    blastResult_obj.sbjct_end = get_intervals_Sbjct.group(2)

                    

 






