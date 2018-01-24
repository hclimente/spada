from collections import Counter

class AminoAcid:
	def __init__(self, resNum, resName):
		"""
		pdbMapping: dictionary with the topological information (B, IS or NIS) of that residue (value)
		in a particular PDB (key).
		"""

		self._num 				= resNum
		self._res 				= resName
		self._tag 				= []
		self._isoformSpecific 	= False
		self._pdbMapping 		= {}
		self._genomicPosition 	= None
		self._disordered 		= False
		self._kinases 	 		= set()
		self._features			= set()

	@property
	def res(self): return self._res
	@property
	def num(self): return self._num
	@property
	def genomicPosition(self): return self._genomicPosition
	def setGenomicPosition(self,genomicPosition): self._genomicPosition = genomicPosition
	@property
	def disordered(self): return self._disordered
	def setDisordered(self, disordered): self._disordered = disordered

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

	def inFeature(self, f):
		return [ x[1] for x in self._features if f == x[0] ]
