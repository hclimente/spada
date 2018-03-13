from collections import Counter

class Aminoacid:
	def __init__(self, resNum, resName):

		self._num 				= resNum
		self._res 				= resName
		self._isoformSpecific 	= False
		self._genomicPosition 	= None
		self._disordered 		= False
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
	def isoformSpecific(self): 	return self._isoformSpecific
	def setIsoformSpecific(self,isoformSpecific): self._isoformSpecific=isoformSpecific

	def inFeature(self, f):
		return [ x[1] for x in self._features if f == x[0] ]
