#!/bin/bash

for folder in  `ls Results/TCGA/*_mE0.0_mCE-1.0/RWorkspaces/2_GetCandidates.RData`
do
	tag=`echo $folder | cut -d'/' -f3 | cut -d'_' -f1`
	Pipeline/scripts/prepareListForDAVID.r $tag $folder
	mv "$tag"_expressedGenes.lst ~/Dropbox/SmartAS/stuff
	mv "$tag"_candidateGenes.lst ~/Dropbox/SmartAS/stuff
done
