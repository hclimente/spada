import argparse
from libs import utils as ut

class ConfigParser():
	def parseOptions(self):
		self._parser = argparse.ArgumentParser(prog="SmartAS.py", description = "Find significant alternative splicing switches. Analyze their functional impact.", epilog= "Hector Climente, 2014")

		self._parser.add_argument('-f','--config-file', dest='config_file', action = 'store', default='Parameters.cfg',
						  help = 'Input file with the parameters')
		self._parser.add_argument('-s', '--initial-step', dest='initial_step', action='store', default='0',
						  type=int, help='Where should SmartAS start.')
		self._parser.add_argument('-wd', '--working-directory', dest='wd', action='store', default='/home/hector/SmartAS/',
						  help='Root file of SmartAS folder in the current machine.')
		self._parser.add_argument('-gwd', '--gaudi-wd', dest='gwd', action='store', default='/sbi/users/hectorc/SmartAS',
						  help='Root file of SmartAS folder in Gaudi.')
		self._parser.add_argument('-m', '--minimum-expression', dest='minExpression', action='store', default='0',
						  type=float, help='Minimum expression to consider a transcript not residual.')
		self._parser.add_argument('-i', '--input-type', dest='inputType', action='store', default='TCGA',
						  help='Origin of the data.')
		self._parser.add_argument('-r', '--replicates', dest='replicates', action='store', default='0',
						  type=int, help='Number of patients or biological replicates.')
		self._parser.add_argument('-v', '--iloops-version', dest='iLoopsVersion', action='store', default='iLoops13',
						  help='Version of iLoops to run.')
		self._parser.add_argument('-u', '--unpaired-replicates', dest='unpairedReplicates', action='store', default='0',
						  type=int, help='Number of unpaired samples.')
		self._parser.add_argument('-t', '--tag', dest='tag', action='store', default='20',
						  help='Tag of the files.')

		options = self._parser.parse_args()
		repr(options)
		#"Conditions" : ["N", "T"]

	def parse_arguments():
		options = a.parse_args()
		return options
	
	def convert_arg_line_to_args(self, arg_line):
		for arg in arg_line.split("="):
			if not arg.strip():
				continue
			yield arg

if __name__ == '__main__':
	x = ConfigParser()