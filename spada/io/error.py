import logging

class SpadaError(Exception):
	def __init__(self, value):
		self.logger = logging.getLogger(__name__)
		self.value = value
		self.logger.exception(value) 
	def __str__(self):
		self.logger.exception(value)
		return repr(self.value)
