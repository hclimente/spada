load("~/SmartAS/Results/TCGA/luad_mE0.0_mCE-1.0/RWorkspaces/2_GetCandidates.RData")

library(ggplot2)

tumCands <- intraReplicate[[8]]$Transcript %in% candidates[[8]]$mindPSI
norCands <- intraReplicate[[8]]$Transcript %in% candidates[[8]]$maxdPSI

png(paste0("~/luad_rep8_allBlack.png"), width=1920, height=1080)
ggplot(NULL) + geom_point(data=intraReplicate[[8]], aes(x=la_tTPM, y=deltaPSI)) + theme_minimal(base_size=30) + guides(fill=FALSE) + ylab("ΔPSI") + xlab("Average (log tTPM)")
graphics.off()

png(paste0("~/luad_rep8_kansurFirebrick2_normalSteelBlue3.png"), width=1920, height=1080)
ggplot(NULL) + geom_point(data=intraReplicate[[8]], colour="gray80", aes(x=la_tTPM, y=deltaPSI)) + geom_point(data=intraReplicate[[8]][tumCands,], colour="firebrick2", aes(x=la_tTPM, y=deltaPSI)) + geom_point(data=intraReplicate[[8]][norCands,], colour="steelblue3", aes(x=la_tTPM, y=deltaPSI))  + theme_minimal(base_size=30) + guides(fill=FALSE) + ylab("ΔPSI") + xlab("Average (log tTPM)")
graphics.off()