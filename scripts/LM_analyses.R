library(ggplot2)
require("grid")

setwd("C:/Users/Héctor/Documents/Trabajo/SmartAS/data/mutations")

source("../../pipeline/scripts/variablesAndFunctions.r")

mutual_exclusion <- read.delim("tables/gene_all_mutations_all_switches_allCancers.txt",row.names=NULL)
colnames(mutual_exclusion) <- c("GeneId","Symbol","Normal_transcript","Tumor_transcript","MS","M","S","N","H","p_me","p.adj_me","p_o","p.adj_o","OR")

switches <- read.delim("../switches/tables/candidateList_allCancers_models_notNoise.txt",row.names=NULL)

switches <- merge(switches,mutual_exclusion)

ggplot(switches) + geom_boxplot(aes(x=as.character(Driver),y=log2(OR)))
ggplot(switches) + geom_histogram(aes(x=log2(OR+0.00001),fill=as.character(Driver)), alpha=0.5)
ggplot(switches) + stat_density(aes(x=-log10(p_o),fill=as.character(Driver)), alpha=0.5)

ggplot(switches) + geom_boxplot(aes(x=as.character(Driver),y=p_o))
ggplot(switches, aes(x=as.character(Driver),y=p_me)) + geom_boxplot(outlier.colour = NA)  + geom_point(position = position_jitter(width = 0.2))
ggplot(switches, aes(x=as.character(Driver),y=p_me)) + geom_violin()
ggplot(switches, aes(x=as.character(Driver),y=p_o)) + geom_violin()
ggplot(switches, aes(x=as.character(Driver),y=log2(OR))) + geom_violin()
ggplot(switches) + geom_histogram(aes(x=-log10(p_o),fill=as.character(Driver)), alpha=0.5)

ggplot(switches) + geom_point(aes(x=-log10(p_o),y=PatientNumber,color=as.character(Driver)))
ggplot(switches) + geom_point(aes(x=-log10(p_me),y=PatientNumber,color=as.character(Driver)))

switches$J <- switches$MS/(switches$MS+switches$M+switches$S)

ggplot(switches) + geom_boxplot(aes(x=as.character(Driver),y=J))
ggplot(switches) + geom_histogram(aes(x=J,fill=as.character(Driver)), alpha=0.5)
ggplot(switches) + stat_density(aes(x=J,fill=as.character(Driver)), alpha=0.5) + scale_y_log10()

ggplot(switches,aes(x=log10(M+MS),y=log10(S+MS),size=log10(MS),shape=as.character(Driver),color=(OR > 1))) + geom_point()

# final plot
df <- switches

df$Predominant <- factor(NA,levels=c("Mutual exclusion","Co-ocurrence"))
df$Predominant[df$p_me<df$p_o] <- "Mutual exclusion"
df$Predominant[df$p_me>df$p_o] <- "Co-ocurrence"

df$p <- NA
df$p[df$p_me<df$p_o] <- df$p_me[df$p_me<df$p_o]
df$p[df$p_me>df$p_o] <- df$p_o[df$p_me>df$p_o]

# use or
p <- ggplot(subset(df, Driver==1),aes(x=log10(M+MS),y=log10(S+MS),size=abs(log2(OR)),color=(OR > 1))) +
  geom_point() +
  labs(x=bquote(log[10] ~ "(# mutations)"), y=bquote( log[10] ~ "(# switches)")) +
  scale_color_discrete("Predominant mechanism", labels=c("Mutual exclusion","Co-ocurrence")) +
  scale_size_continuous(bquote(bold("|"~log[2](OR)~"|"))) +
  smartas_theme() +
  theme(legend.position="bottom")
ggsave("figures/ms_vs_o_OR.png",p)

p <- ggplot(subset(df, Driver==1 & IsRelevant==1),aes(x=log10(M+MS),y=log10(S+MS),size=abs(log2(OR)),color=(OR > 1))) +
  geom_point() +
  labs(x=bquote(log[10] ~ "(# mutations)"), y=bquote( log[10] ~ "(# switches)")) +
  scale_color_discrete("Predominant mechanism", labels=c("Mutual exclusion","Co-ocurrence")) +
  scale_size_continuous(bquote(bold("|"~log[2](OR)~"|"))) +
  smartas_theme() +
  theme(legend.position="bottom")
ggsave("figures/ms_vs_o_OR_functional.png",p)

# use p
p <- ggplot() +
  geom_point(data=subset(df, Driver==1),aes(x=log10(M+MS),y=log10(S+MS),size=-log10(p),color=Predominant)) +
  geom_text(data=subset(df, Driver==1 & p < 0.05),aes(x=log10(M+MS)+0.1,y=log10(S+MS)+0.1,label=Symbol,color=Predominant)) +
  labs(x=bquote(log[10] ~ "(# mutations)"), y=bquote( log[10] ~ "(# switches)")) +
  scale_color_discrete("Predominant mechanism", labels=c("Mutual exclusion","Co-ocurrence")) +
  scale_size_continuous(bquote(bold(-log[10]~"(p)  "))) +
  smartas_theme() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(override.aes = list(shape = 15, size = 9)))
ggsave("figures/ms_vs_o_p.png",p)

p <- ggplot() +
  geom_point(data=subset(df, Driver==1 & IsRelevant==1),aes(x=log10(M+MS),y=log10(S+MS),size=-log10(p),color=Predominant)) +
  geom_text(data=subset(df, Driver==1 & p < 0.05 & IsRelevant==1),aes(x=log10(M+MS)+0.1,y=log10(S+MS)+0.1,label=Symbol,color=Predominant)) +
  labs(x=bquote(log[10] ~ "(# mutations)"), y=bquote( log[10] ~ "(# switches)")) +
  scale_color_discrete("Predominant mechanism", labels=c("Mutual exclusion","Co-ocurrence")) +
  scale_size_continuous(bquote(bold(-log[10]~"(p)  "))) +
  smartas_theme() +
  theme(legend.position="bottom") +
  guides(col = guide_legend(override.aes = list(shape = 15, size = 9)))
ggsave("figures/ms_vs_o_p_functional.png",p)