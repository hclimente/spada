from libs import iLoops_xml_parser
from libs import options

import tarfile
import os

class iLoopsParser(iLoops_xml_parser.ILXMLParser):
	def custom_protein_output(self, protein_object, **kwds): 
		return protein_object

	def custom_interaction_output(self, interaction_object, **kwds): 
		return interaction_object
	
	def parseLoops(self, xmlOutput, **kwds):
		parsedLoops = {}

		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
			if isinstance(resultItem, iLoops_xml_parser.ILXMLProtein):
				if resultItem.get_loops():
					loopList = [ aLoop.get_code() for aLoop in resultItem.get_loops() ]
					loopList.sort() # sorted list by loop name
					parsedLoops[resultItem.get_name()] = ";".join(loopList)

		return parsedLoops

	def parseNoLoops(self, iLoopsPath, xmlOutput, **kwds):
		noLoops = []
		
		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
			if isinstance(resultItem, iLoops_xml_parser.ILXMLProtein):
				noLoops.append(resultItem.get_name())

		return noLoops

	def parseInteractions(self,tx,expressedTranscripts):
		interactions = {}
		partner = ""
		tmpRC = 0

		tarFile = "{0}iLoops/{1}/{2}/{3}.tar.gz".format(options.Options().wd, 
													 options.Options().inputType, 
													 options.Options().iLoopsVersion,tx)
		xmlFile = "{0}{1}.ips".format(options.Options().qout,tx)

		tar = tarfile.open(tarFile)
		tar.extract(tx+".ips", path=options.Options().qout)

		tar.close()
		with open(xmlFile) as FILE:
			for line in FILE:
				strippedLine = line.strip()
				if strippedLine == "<interaction>":
					partner = ""
					tmpRC = 0
				elif "RF_cost" in strippedLine:
					if int(strippedLine[9:-10]) > tmpRC:
						tmpRC = int(strippedLine[9:-10])
				elif "RF_prediction" in strippedLine:
					interactions.setdefault(partner, 0)
					if "True" in strippedLine and tmpRC > interactions[partner]:
						interactions[partner] = tmpRC
				elif ("P1ID" in strippedLine or "P2ID" in strippedLine):
					if not partner or partner == tx:
						partner = strippedLine[6:-7]

		os.remove(xmlFile)

		expressedInteractions = {k: interactions[k] for k in expressedTranscripts if k in interactions }

		return expressedInteractions

	def makeProteinsLite(self, xmlOutput, xmlOriginal, **kwds):
		for resultItem in self.results_parser(xml_file=xmlOriginal, report_level=0, **kwds):
			if not isinstance(resultItem, iLoops_xml_parser.ILXMLProtein):continue

			xmlOutput.getWriter().startTag("protein")

			xmlOutput._open_close_tag_with_info(tag_name='name', value=resultItem.get_name())
			xmlOutput._open_close_tag_with_info(tag_name='MD5', value=resultItem.get_MD5())
			xmlOutput._open_close_tag_with_info(tag_name='sequence', value=resultItem.get_sequence())

			for feature in resultItem.get_loops():
				xmlOutput.getWriter().startTag("feature")
				xmlOutput._open_close_tag_with_info(tag_name="type", value="loop")
				xmlOutput._open_close_tag_with_info(tag_name="loop_subClass", value=feature.get_subClass())
				xmlOutput._open_close_tag_with_info(tag_name="loop_match", value=feature.get_loopMatch())
				xmlOutput._open_close_tag_with_info(tag_name="align_ini", value=feature.get_ali_ini())
				xmlOutput._open_close_tag_with_info(tag_name="align_end", value=feature.get_ali_end())
				xmlOutput._open_close_tag_with_info(tag_name="query_feature_ali_seq", value=feature.get_QAliSeq())
				xmlOutput._open_close_tag_with_info(tag_name="target_feature_ali_seq", value=feature.get_HAliSeq())

				if len(feature.get_mappings()) < 1: xmlOutput.getWriter().tag('SCOP_assignations')
				else: 
					xmlOutput.getWriter().startTag('SCOP_assignations')

					scop_dom = feature.get_mappings()[0]
			 		if not scop_dom.get_ID() == "out_of_domain" and scop_dom.get_ID() is not 'None':

			 			xmlOutput.getWriter().startTag('SCOP_assignation')
			 			xmlOutput._open_close_tag_with_info(tag_name='SCOP_id', value=scop_dom.get_ID())
			 			xmlOutput._open_close_tag_with_info(tag_name='SCOP_ini', value=scop_dom.get_ini())
			 			xmlOutput._open_close_tag_with_info(tag_name='SCOP_end', value=scop_dom.get_end())
			 			xmlOutput.getWriter().endTag('SCOP_assignation')
			 				
			 		xmlOutput.getWriter().endTag('SCOP_assignations')

				xmlOutput.getWriter().endTag('feature')

			xmlOutput.getWriter().endTag("protein")

	def makeInteractionsLite(self, xmlOutput, xmlOriginal, **kwds):
		for resultItem in self.results_parser(xml_file=xmlOriginal, report_level=0, **kwds):
			if not isinstance(resultItem, iLoops_xml_parser.ILXMLInteraction):
				continue

			xmlOutput.getWriter().startTag("interaction")

			xmlOutput._open_close_tag_with_info(tag_name='P1ID', value=resultItem.get_i1name())
			xmlOutput._open_close_tag_with_info(tag_name='P2ID', value=resultItem.get_i2name())
			xmlOutput._open_close_tag_with_info(tag_name='Splus', value=resultItem.get_Splus())
			xmlOutput._open_close_tag_with_info(tag_name='Sminus', value=resultItem.get_Sminus())
			xmlOutput._open_close_tag_with_info(tag_name='LSR', value=resultItem.get_LSR())
			xmlOutput._open_close_tag_with_info(tag_name='LpVR', value=resultItem.get_LpVR())
			
			xmlOutput.getWriter().startTag('RF_results')
			for rfResultItem in resultItem.get_RFResults():
				xmlOutput.getWriter().startTag('RF_result')
				xmlOutput._open_close_tag_with_info(tag_name='RF_cost', value=rfResultItem.get_cost() )
				xmlOutput._open_close_tag_with_info(tag_name='RF_prediction', value=rfResultItem.get_prediction() )
				xmlOutput._open_close_tag_with_info(tag_name='RF_score', value=rfResultItem.get_RFscore() )
				xmlOutput.getWriter().endTag('RF_result')
			xmlOutput.getWriter().endTag('RF_results')

			xmlOutput.getWriter().endTag("interaction")