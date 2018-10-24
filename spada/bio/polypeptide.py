class Polypeptide:
	def __init__(self, segment):

		self._structure	= segment

	def __repr__(self):
		return 'Segment of {} residues: {}'.format(len(self), ''.join([ x.res for x in self._structure ]) )

	def __len__(self):
		return len(self._structure)

	def __and__(self, other):
		return Polypeptide(list(set(self._structure) & set(other._structure)))

	def __or__(self, other):
		return Polypeptide(list(set(self._structure) | set(other._structure)))

	def __iter__(self):
		for residue in self._structure:
			yield residue