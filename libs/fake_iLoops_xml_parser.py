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

		return noLoops

	def parseInteractions(self, thisCandidate, xmlOutput, **kwds):

		KK = open(thisCandidate + ".tsv", "w")
		KK.write("Partner\tSplus\tSminus\tLSR\tLpVR\tmaxCost\n")

		for resultItem in self.results_parser(xml_file=xmlOutput, report_level=0, **kwds):
			mc = 0
			if isinstance(resultItem, iLoops_xml_parser.ILXMLInteraction):
				intPartner = resultItem.get_i2name()
				KK.write( str(intPartner) + "\t" + str( resultItem.get_Splus() ) + "\t" + str( resultItem.get_Sminus() ) + "\t" + str( resultItem.get_LpVR() ) + "\t" + str( resultItem.get_LSR() ) + "\t")

				for RFResult in resultItem.get_RFResults(): 
					if RFResult.get_prediction() and RFResult.get_cost() > mc:
						mc = int(RFResult.get_cost())

				KK.write( str(mc) + "\n")