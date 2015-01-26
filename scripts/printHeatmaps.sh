#!/bin/bash

for folder in  `ls SmartAS/testResults/TCGA/*_mE0.0_mCE-1.0/RWorkspaces/2_GetCandidates.RData`
do
	tag=`echo $folder | cut -d'/' -f4 | cut -d'_' -f1`
	echo $tag
	SmartAS/Pipeline/scripts/printHeatmaps.r $tag $folder
	mv "$tag"_heatmap.png ~/Dropbox/SmartAS/stuff
done