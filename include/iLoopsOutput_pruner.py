#!/soft/devel/python-2.7/bin/python

import loxun
import custom_iLoops_xml_parser as parser
from os import listdir
from fnmatch import filter
import tarfile

class xmlLite:
	def __init__(self, transcript, workingDirectory):
		self.transcriptName = transcript
		self.wd 			= workingDirectory + "/"
		self.outFile 		= open(self.wd + self.transcriptName + ".xml", "w")
		self.writer			= loxun.XmlWriter(self.outFile)

		self._init_file()

	def getWd(self):		return self.wd
	def getTxName(self):	return self.transcriptName
	def getWriter(self):	return self.writer
	def _init_file(self): self.getWriter().startTag("xml")
	def _end_file(self): 
		self.getWriter().endTag('xml')
		self.getWriter().close()

	def _open_close_tag_with_info(self, tag_name, value):
		self.getWriter().startTag(tag_name)
		self.getWriter().text(str(value))
		self.getWriter().endTag(tag_name)

	def makeLiteVersion(self):
		xmlFile = self.getWd() + self.transcriptName + "_raw.ips"
		self.getWriter().startTag("xml")
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
		self.getWriter().endTag("xml")
		self.getWriter().close()

	def joinFiles(self):
		path = self.getWd() + "/" + self.getTxName() + "/"
		with open(self.getWd() + self.transcriptName + "_raw.ips", "w") as RAW_XML_JOIN:
			RAW_XML_JOIN.write("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
			RAW_XML_JOIN.write("<xml>\n")

			maxBatch = len( filter( listdir(path), transcript + "_[0-9]") )
				
			for batch in [ str(n) for n in range(1, maxBatch + 1) ]:
				candidateOut = path + transcript + "_" + batch + "/sge_output"
				for number in [ str(x) + str(y) for x in [""] + range(10) for y in range(10) if str(x) + str(y) != "00" ]:
					for nodeMap in filter( listdir( candidateOut ), "*.assignation." + number + ".xml" ):
						proteinsFile = candidateOut + "/" + nodeMap
						with open(proteinsFile, "r") as PROTEINS_FILE:
							for line in PROTEINS_FILE:
								if line.strip() != "<?xml version=\"1.0\" encoding=\"utf-8\"?>" and line.strip() != "<xml>" and line.strip() != "</xml>":
									RAW_XML_JOIN.write(line)

			for batch in [ str(n) for n in range(1, maxBatch + 1) ]:
				candidateOut = path + transcript + "_" + batch + "/sge_output"
				for number in [ str(x) + str(y) for x in [""] + range(10) for y in range(10) if str(x) + str(y) != "00" ]:
					for nodeInt in filter(listdir( candidateOut ), "*.scoring." + number + ".xml" ):
						interactionsFile = candidateOut + "/" + nodeInt
						with open(interactionsFile, "r") as INTERACTIONS_FILE:
							for line in INTERACTIONS_FILE:
								if line.strip() != "<?xml version=\"1.0\" encoding=\"utf-8\"?>" and line.strip() != "<xml>" and line.strip() != "</xml>":
									RAW_XML_JOIN.write(line)

			RAW_XML_JOIN.write("</xml>\n")

	def makeTarfile(self):
		with tarfile.open(self.getWd() + self.transcriptName + ".tar.gz", "w:gz") as tar:
			tar.add(self.getWd() + "/" + self.transcriptName + ".xml", arcname = self.transcriptName + ".xml")

if __name__ == '__main__':
	
	transcript = "testTranscript"
	workingDirectory = "/home/hector/Desktop"

	r = xmlLite(transcript, workingDirectory)
	#r.joinFiles()
	r.makeLiteVersion()
	r.makeTarfile()