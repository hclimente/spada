# Preparing to run SPADA

The first step is making sure you meet the requirements i.e. you installed all SPADA's dependencies, and downloaded all the required databases. Read more [here](requirements.md).

# Import your data
In order to import and analyze your externally calculated switches, you need the following things:

* A list of your switches
* A parent SmartAS experiment to copy the networks from, with the same notation
* A list of expressed transcripts

The list of switches is a tab-delimited table with a header and a line per switch containing the gene name, the tumor isoform and the normal isoform. Its path is specified using the argument external-switches e.g. external-switches=my_switches.tsv. The table has the following format:

```
#!text

# One line header e.g. Gene	Tumor_isoform	NormalIsoform
Gene1	tumor_isoform	normal_isoform
Gene2	tumor_isoform	normal_isoform

```
The parent SmartAS experiment, identified by a tag, needs to be located under the same annotation directory in the working directory. It must me specified using the argument parent-tag e.g. parent-tag=brca. The new experiment will copy both the gene (geneNetwork.pkl) and the transcript (txNetwork.pkl) networks from its parent.
