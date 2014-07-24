from collections import Counter
import logging

class AminoAcid:
	def __init__(self, resNum, resName):
		self.logger 			= logging.getLogger(__name__)
		self._num 				= resNum
		self._res 				= resName
		self._tag 				= []
		self._isoformSpecific 	= False
		self._pdbMapping 		= {}
		self._genomicPosition 	= None

	@property
	def res(self): return self._res
	@property
	def genomicPosition(self): return self._genomicPosition
	def setGenomicPosition(self,genomicPosition): self._genomicPosition = genomicPosition
	@property
	def tag(self): 
		"""Returns the mode from the list of tags."""
		frequency = Counter(self._tag).most_common(1)
		if frequency:
			return frequency[0][0]
		else:
			return None

	def setTag(self,tag): self._tag.append(tag)

	@property
	def isoformSpecific(self): 	return self._isoformSpecific
	
	def setIsoformSpecific(self,isoformSpecific): self._isoformSpecific=isoformSpecific
