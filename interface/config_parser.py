import argparse
from libs import utils as ut

class ConfigParser():
	def parseOptions(self):
		pass

	def parse_arguments():
		options = a.parse_args()
		return options
	
	def convert_arg_line_to_args(self, arg_line):
		for arg in arg_line.split("="):
			if not arg.strip():
				continue
			yield arg

if __name__ == '__main__':
	pass