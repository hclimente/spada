#!/soft/devel/python-2.7/bin/python

from libs.utils import *

class Method:
	def __init__(self, path="", samples=None, candidates={}):
		self._path			= path
		self._samples		= samples
		self._candidates	= candidates

	def path(self):			return self._path
	def samples(self):		return self._samples
	def candidates(self):	return self._candidates