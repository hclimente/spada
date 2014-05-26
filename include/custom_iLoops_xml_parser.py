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
						loopList.append( (aLoop.get_code(), [mapping.get_ID() for mapping in aLoop.get_mappings()]) )
					loopList.sort(key=lambda x: x[0]) # sorted list by loop name
					parsedLoops[resultItem.get_name()] = ";".join([ x[0]+"_"+"?".join([dom for dom in x[1]]) for x in loopList])

		return parsedLoops

	def parseNoLoops(self, iLoopsPath, xmlOutput, **kwds):
		noLoops = []
		
		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
			if isinstance(resultItem, iLoops_xml_parser.ILXMLProtein):
				noLoops.append(resultItem.get_name())

		return noLoops

	def parseInteractions(self, thisCandidate, expressedIsoforms, xmlOutput, **kwds):
		maxCost = {}

		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds): 
			if isinstance(resultItem, iLoops_xml_parser.ILXMLInteraction):
				intPartner = resultItem.get_i2name()
				if intPartner == thisCandidate:
					intPartner = resultItem.get_i1name()

				if intPartner not in expressedIsoforms:
					continue
				
				maxCost.setdefault(intPartner, 0)
				for RFResult in resultItem.get_RFResults(): 
					if RFResult.get_prediction() and RFResult.get_cost() > maxCost[intPartner]:
						maxCost[intPartner] = int(RFResult.get_cost())

		return maxCost