import logging

class SpadaError(Exception):
	def __init__(self, value):
		self.logger = logging.getLogger('spada error')
		self.value = value
		self.logger.exception(value)
	def __str__(self):
		return repr(self.value)
