#!/soft/devel/python-2.7/bin/python

import loxun
import interface.iLoops_parser as parser
from os import listdir
from fnmatch import filter
import tarfile
import sys

class iLoopsOutput_pruner:
	def __init__(self, transcript, workingDirectory):
		self.transcriptName = transcript
		self.wd 			= workingDirectory + "/"
		
		self._init_file()

	def getWd(self):		return self.wd
	def getTxName(self):	return self.transcriptName
	def getWriter(self):	return self.writer
	def getOutFile(self):	return self.outFile

	def _init_file(self):
		self.outFile 		= open(self.getWd() + self.getTxName() + ".ips", "w")
		self.writer			= loxun.XmlWriter(self.outFile)
		self.getWriter().startTag("xml")
	def _end_file(self): 
		self.getWriter().endTag('xml')
		self.getWriter().close()
	def _change_pretty_state(self, state):
		self.getWriter()._pretty = state

	def _open_close_tag_with_info(self, tag_name, value):
		self.getWriter().startTag(tag_name)
		self._change_pretty_state(False)

		self.getOutFile().write(self.getWriter()._indent * (len(self.getWriter()._elementStack)-1))
		self.getOutFile().flush()
		self.getWriter().text(str(value))
		self.getWriter().endTag(tag_name)
		self.getWriter().newline()
		self._change_pretty_state(True)

	def makeLiteVersion(self):
		xmlFile = self.getWd() + self.getTxName() + "_raw.ips"
		myParser = parser.iLoopsParser()
		myParser.makeProteinsLite(
									xmlOutput					  = self,
									xmlOriginal					  = xmlFile,
									output_proteins               = True, 
									output_alignments             = False,
									output_domain_mappings        = False,
									output_protein_features       = True,
									output_domain_assignations    = True,
									output_interactions           = False,
									output_interaction_signatures = False,
									output_RF_results             = False,
									output_RF_precisions          = False
								 )

		myParser.makeInteractionsLite(
										xmlOutput					  = self,
										xmlOriginal					  = xmlFile,
										output_proteins               = False, 
										output_alignments             = False,
										output_domain_mappings        = False,
										output_protein_features       = False,
										output_domain_assignations    = False,
										output_interactions           = True,
										output_interaction_signatures = False,
										output_RF_results             = True,
										output_RF_precisions          = True
									  )
		self._end_file()
		self.checkConsistency()
		self.makeTarfile()

		return True

	def checkConsistency(self):
		outFile = self.getWd() + self.getTxName() + ".ips"

		myParser = parser.iLoopsParser()
		myParser.results_parser(
										xml_file					  = outFile, 
										report_level				  = 0, 
										output_proteins               = True, 
										output_alignments             = True,
										output_domain_mappings        = True,
										output_protein_features       = True,
										output_domain_assignations    = True,
										output_interactions           = True,
										output_interaction_signatures = True,
										output_RF_results             = True,
										output_RF_precisions          = True
									  )

	def joinFiles(self):
		with open(self.getWd() + self.getTxName() + "_raw.ips", "w") as RAW_XML_JOIN:
			RAW_XML_JOIN.write("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
			RAW_XML_JOIN.write("<xml>\n")

			maxBatch = len( filter( listdir(self.getWd()), self.getTxName() + "_[0-9]") )
				
			for batch in [ str(n) for n in range(1, maxBatch + 1) ]:
				candidateOut = self.getWd() + self.getTxName() + "_" + batch + "/sge_output"
				for number in [ str(x) + str(y) for x in [""] + range(10) for y in range(10) if str(x) + str(y) != "00" ]:
					for nodeMap in filter( listdir( candidateOut ), "*.assignation." + number + ".xml" ):
						proteinsFile = candidateOut + "/" + nodeMap
						with open(proteinsFile, "r") as PROTEINS_FILE:
							for line in PROTEINS_FILE:
								if line.strip() != "<?xml version=\"1.0\" encoding=\"utf-8\"?>" and line.strip() != "<xml>" and line.strip() != "</xml>":
									RAW_XML_JOIN.write(line)

			for batch in [ str(n) for n in range(1, maxBatch + 1) ]:
				candidateOut = self.getWd() + self.getTxName() + "_" + batch + "/sge_output"
				for number in [ str(x) + str(y) for x in [""] + range(10) for y in range(10) if str(x) + str(y) != "00" ]:
					for nodeInt in filter(listdir( candidateOut ), "*.scoring." + number + ".xml" ):
						interactionsFile = candidateOut + "/" + nodeInt
						with open(interactionsFile, "r") as INTERACTIONS_FILE:
							for line in INTERACTIONS_FILE:
								if line.strip() != "<?xml version=\"1.0\" encoding=\"utf-8\"?>" and line.strip() != "<xml>" and line.strip() != "</xml>":
									RAW_XML_JOIN.write(line)

			RAW_XML_JOIN.write("</xml>\n")

	def makeTarfile(self):
		tar = tarfile.open(self.getWd() + self.getTxName() + ".tar.gz", "w:gz")
		tar.add(self.getWd() + "/" + self.getTxName() + ".ips", arcname = self.getTxName() + ".ips")
		tar.close()

if __name__ == '__main__':
	
	transcript = sys.argv[1]
	workingDirectory = sys.argv[2]

	r = iLoopsOutput_pruner(transcript, workingDirectory)
	r.joinFiles()
	r.makeLiteVersion()