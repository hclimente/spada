from spada.methods import create_network
import os

scriptPath = os.path.realpath(__file__)
dataPath = os.path.dirname(scriptPath) + "/data"

c = create_network.CreateNetwork("test", "gencode")
c.run("{}/gtf".format(dataPath),
	  "{}/expression".format(dataPath),
	  "{}/expression".format(dataPath),
	  -3.3,
	  "{}/fasta".format(dataPath),
	  "{}/mitab".format(dataPath),
	  "{}/drivers".format(dataPath),
	  "{}/features".format(dataPath))
