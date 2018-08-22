import logging
import warnings

class SpadaError(Exception):
	def __init__(self, value):
		self.logger = logging.getLogger(__name__)
		self.value = value
		self.logger.exception(value)
	def __str__(self):
		return repr(self.value)

class SpadaWarning(RuntimeWarning):
	def __init__(self, value):
		self.logger = logging.getLogger(__name__)
		self.value = value
		self.logger.warning(value) 
	def __str__(self):
		return repr(self.value)
