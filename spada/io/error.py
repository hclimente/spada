import logging

class SpadaError(Exception):
	def __init__(self, value):
		self.logger = logging.getLogger('spada error')
		self.value = value
	def __str__(self):
		self.logger.exception(value)
		return repr(self.value)
