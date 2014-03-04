import iLoops_xml_parser

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
					loopList = []
					for aLoop in resultItem.get_loops():
						loopList.append(aLoop.get_code())
					loopList.sort()
					parsedLoops[resultItem.get_name()] = ";".join(loopList)

		return parsedLoops

	def parseNoLoops(self, iLoopsPath, xmlOutput, **kwds):
		noLoops = []
		
		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
			if isinstance(resultItem, iLoops_xml_parser.ILXMLProtein):
				noLoops.append(resultItem.get_name())

	def parseInteractions(self, thisCandidate, xmlOutput, **kwds):
		interactions = []

		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
			if isinstance(resultItem, parser.ILXMLInteraction):
				if thisCandidate == resultItem.get_interactor1_ID():
					print resultItem.get_interactor2_ID()
					interactions.append(resultItem.get_interactor2_ID())
				elif thisCandidate == resultItem.get_interactor2_ID():
					print resultItem.get_interactor1_ID()
					interactions.append(resultItem.get_interactor1_ID())

		return interactions