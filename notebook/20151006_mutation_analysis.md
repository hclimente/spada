-   Specific switches: overlap in features affected by mutations and switches
    -   Analysis per tumor type
    -   Analysis aggregating all tumors
    -   Comparison with meScores
-   General overlap in features affected by mutations and switches
    -   View of the dat
        -   Analysis per tumor type
    -   Analysis aggregating all tumors
    -   Cluster in domain function

    ## Loading required package: grid
    ## Loading required package: quadprog

    ## Warning: package 'gridExtra' was built under R version 3.2.2

We studied how switches and mutations affect the same features of the genes. In this study, we will consider protein-affecting mutations only. The goal is to find features that might be activated or inactivated in tumor both through mutations and splicing changes.

Specific switches: overlap in features affected by mutations and switches
=========================================================================

Make a binomial test to find switched features that are frequently mutated. The expected frequency is the relative size of the feature in the isoform (between 0 and 1); the obtained frequency is the proportion of mutations in the gene that affect that feature.

Analysis per tumor type
-----------------------

``` r
mut_feat_overlap <- read.delim(paste0(wd,"mutation_switch_feature_overlap.txt"), header=TRUE)
colnames(mut_feat_overlap) <- c("Gene","Symbol","Cancer","Normal_transcript","Tumor_transcript","What","FeatureType","Feature","Ratio","Driver","TotalMutations","MutationsInFeature","FeatureSize")
pvals <- apply(mut_feat_overlap[,c("MutationsInFeature","TotalMutations","FeatureSize")],1, function(x){ 
  if (x[2]==0){
    p <- 1
  } else {
    p <- binom.test(x[1],x[2],x[3],"greater")
    p <- p$p.value
  }
  p
} )
adjpvalues <- p.adjust(pvals)

mut_feat_overlap$p <- pvals
mut_feat_overlap$adjp <- adjpvalues
mut_feat_overlap <- mut_feat_overlap[order(mut_feat_overlap$p),]

write.table(mut_feat_overlap,paste0(wd,'tables/mutation_switch_feature_overlap_withPVals.txt'),quote=F,col.names=T,sep="\t",row.names=FALSE)
```

| Cancer | Symbol    | What              | Feature                                                                   | Driver |      p|   adjp|
|:-------|:----------|:------------------|:--------------------------------------------------------------------------|:-------|------:|------:|
| luad   | GLG1      | Lost\_in\_tumor   | PS51289|GLG1\_C\_RICH\_Cysteine-rich\_GLG1\_repeat\_profile.              | False  |  0.000|  0.125|
| hnsc   | OGT       | Lost\_in\_tumor   | PS50005|TPR\_TPR\_repeat\_profile.                                        | False  |  0.000|  0.854|
| hnsc   | OGT       | Lost\_in\_tumor   | PF13414|TPR\_repeat                                                       | False  |  0.000|  1.000|
| luad   | RAB11FIP4 | Gained\_in\_tumor | PS50222|EF\_HAND\_2\_EF-hand\_calcium-binding\_domain\_profile.           | False  |  0.000|  1.000|
| brca   | CTCF      | Lost\_in\_tumor   | PF13465|Zinc-finger\_double\_domain                                       | True   |  0.000|  1.000|
| brca   | APAF1     | Lost\_in\_tumor   | PF00400|WD\_domain,\_G-beta\_repeat                                       | True   |  0.000|  1.000|
| coad   | DUS3L     | Lost\_in\_tumor   | PS50103|ZF\_C3H1\_Zinc\_finger\_C3H1-type\_profile.                       | False  |  0.000|  1.000|
| kirc   | ARHGAP5   | Lost\_in\_tumor   | PS51676|FF\_FF\_domain\_profile.                                          | False  |  0.000|  1.000|
| brca   | CHD6      | Lost\_in\_tumor   | PF00385|Chromo\_(CHRromatin\_Organisation\_MOdifier)\_domain              | False  |  0.001|  1.000|
| kirc   | ZNF644    | Lost\_in\_tumor   | PF13894|C2H2-type\_zinc\_finger                                           | False  |  0.001|  1.000|
| kirc   | ZNF644    | Lost\_in\_tumor   | PS50157|ZINC\_FINGER\_C2H2\_2\_Zinc\_finger\_C2H2\_type\_domain\_profile. | False  |  0.001|  1.000|
| coad   | ADAMTSL1  | Lost\_in\_tumor   | PF00090|Thrombospondin\_type\_1\_domain                                   | False  |  0.001|  1.000|
| hnsc   | LGR6      | Lost\_in\_tumor   | PF13855|Leucine\_rich\_repeat                                             | True   |  0.002|  1.000|
| kirc   | MEGF6     | Gained\_in\_tumor | PS50026|EGF\_3\_EGF-like\_domain\_profile.                                | False  |  0.003|  1.000|
| kirc   | SETDB1    | Lost\_in\_tumor   | PF05033|Pre-SET\_motif                                                    | True   |  0.003|  1.000|
| brca   | LIMS2     | Lost\_in\_tumor   | PF00412|LIM\_domain                                                       | False  |  0.004|  1.000|
| brca   | WDR59     | Lost\_in\_tumor   | PF00400|WD\_domain,\_G-beta\_repeat                                       | False  |  0.004|  1.000|
| luad   | PPP1R9A   | Gained\_in\_tumor | PF00595|PDZ\_domain\_(Also\_known\_as\_DHR\_or\_GLGF)                     | False  |  0.004|  1.000|
| lihc   | GNB5      | Lost\_in\_tumor   | PF00400|WD\_domain,\_G-beta\_repeat                                       | False  |  0.004|  1.000|
| kirc   | RAP1GDS1  | Lost\_in\_tumor   | PF00514|Armadillo/beta-catenin-like\_repeat                               | True   |  0.004|  1.000|
| kich   | ARHGAP5   | Lost\_in\_tumor   | PS51676|FF\_FF\_domain\_profile.                                          | False  |  0.004|  1.000|
| luad   | OAS2      | Lost\_in\_tumor   | PF10421|2'-5'-oligoadenylate\_synthetase\_1,\_domain\_2,\_C-terminus      | False  |  0.005|  1.000|
| brca   | WDR59     | Lost\_in\_tumor   | PS50082|WD\_REPEATS\_2\_Trp-Asp\_(WD)\_repeats\_profile.                  | False  |  0.006|  1.000|
| luad   | ABCA5     | Lost\_in\_tumor   | PF12698|ABC-2\_family\_transporter\_protein                               | False  |  0.006|  1.000|
| kirc   | EGFR      | Lost\_in\_tumor   | PF07714|Protein\_tyrosine\_kinase                                         | True   |  0.006|  1.000|
| kirc   | SDK1      | Lost\_in\_tumor   | PS50835|IG\_LIKE\_Ig-like\_domain\_profile.                               | False  |  0.006|  1.000|
| luad   | C2        | Lost\_in\_tumor   | PF00084|Sushi\_domain\_(SCR\_repeat)                                      | False  |  0.006|  1.000|
| brca   | PREX1     | Lost\_in\_tumor   | PF00610|Domain\_found\_in\_Dishevelled,\_Egl-10,*and\_Pleckstrin*(DEP)    | False  |  0.006|  1.000|
| brca   | ABR       | Lost\_in\_tumor   | PF00168|C2\_domain                                                        | False  |  0.007|  1.000|
| coad   | KIAA0892  | Lost\_in\_tumor   | PS50310|ALA\_RICH\_Alanine-rich\_region\_profile.                         | False  |  0.007|  1.000|

Some interesting cases:

-   EGFR, due to the loss of the kinase domain. Also MEGF6.

-   Thrombospondins are secreted proteins with antiangiogenic abilities.

-   Some domains related to transcription

> The **FF domain** is present in a variety of nuclear transcription and splicing factors, as well as the p190 family of RhoGAPs.

> **Zinc finger (Znf) domains** are relatively small protein motifs which contain multiple finger-like protrusions that make tandem contacts with their target molecule. [...] [T]hey are now recognised to bind DNA, RNA, protein and/or lipid substrates

> **WD40-repeat** proteins are a large family found in all eukaryotes and are implicated in a variety of functions ranging from signal transduction and transcription regulation to cell cycle control, autophagy and apoptosis.

Analysis aggregating all tumors
-------------------------------

``` r
mut_feat_overlap_agg <- ddply(mut_feat_overlap,
                              .(Gene,Symbol,Normal_transcript,Tumor_transcript,What,
                                FeatureType,Feature,Driver,FeatureSize), 
                              summarise, inMut=sum(MutationsInFeature), 
                              totalMut=sum(TotalMutations) )
mut_feat_overlap_agg$Ratio = 100 * mut_feat_overlap_agg$inMut/mut_feat_overlap_agg$totalMut

pvals <- apply(mut_feat_overlap_agg[,c("inMut","totalMut","FeatureSize")],1, function(x){ 
  if (x[2]==0){
    p <- 1
  } else {
    p <- binom.test(x[1],x[2],x[3],"greater")
    p <- p$p.value
  }
  p
} )
adjpvalues <- p.adjust(pvals)

mut_feat_overlap_agg$p_mutation_feature_overlap <- pvals
mut_feat_overlap_agg$adjp_mutation_feature_overlap <- adjpvalues
mut_feat_overlap_agg <- mut_feat_overlap_agg[order(mut_feat_overlap_agg$p_mutation_feature_overlap),]

write.table(mut_feat_overlap_agg,paste0(wd,'tables/mutation_switch_feature_overlap_allCancers_withPVals.txt'),quote=F,col.names=T,sep="\t",row.names=FALSE)
```

    ## Warning: Removed 652 rows containing missing values (geom_point).

![*Mutations that affect a particular structural feature vs size of that particular feature. In purple, significant cases before multitesting correction; in orange, significant cases after multitesting correction.*](20151006_mutation_analysis_files/figure-markdown_github/unnamed-chunk-3-1.png)

Comparison with meScores
------------------------

There is a big overlap with the previous list regarding the top genes. It will be interesting to check if these mutations tend to appear in the same patients as the switch of in different ones. In order to do that, we will show the meScore for those genes:

| Symbol   | What              | Feature                                                                                  | Driver |   Score|
|:---------|:------------------|:-----------------------------------------------------------------------------------------|:-------|-------:|
| EPS15    | Lost\_in\_tumor   | PS50031|EH\_EH\_domain\_profile.                                                         | True   |   1.257|
| TDRD7    | Lost\_in\_tumor   | PS51644|HTH\_OST\_OST-type\_HTH\_domain\_profile.                                        | False  |   0.923|
| MCC      | Gained\_in\_tumor | PS50324|SER\_RICH\_Serine-rich\_region\_profile.                                         | False  |   0.816|
| CNTN4    | Lost\_in\_tumor   | PF07679|Immunoglobulin\_I-set\_domain                                                    | False  |   0.725|
| OAS2     | Lost\_in\_tumor   | PF10421|2'-5'-oligoadenylate\_synthetase\_1,\_domain\_2,\_C-terminus                     | False  |   0.631|
| CD300A   | Lost\_in\_tumor   | PS50835|IG\_LIKE\_Ig-like\_domain\_profile.                                              | False  |   0.627|
| CD300A   | Lost\_in\_tumor   | PF07686|Immunoglobulin\_V-set\_domain                                                    | False  |   0.627|
| TSC22D1  | Gained\_in\_tumor | PS50322|GLN\_RICH\_Glutamine-rich\_region\_profile.                                      | False  |   0.616|
| GNB5     | Lost\_in\_tumor   | PS50082|WD\_REPEATS\_2\_Trp-Asp\_(WD)\_repeats\_profile.                                 | False  |   0.510|
| GNB5     | Lost\_in\_tumor   | PF00400|WD\_domain,\_G-beta\_repeat                                                      | False  |   0.510|
| C2       | Lost\_in\_tumor   | PS50923|SUSHI\_Sushi/CCP/SCR\_domain\_profile.                                           | False  |  -0.509|
| C2       | Lost\_in\_tumor   | PF00084|Sushi\_domain\_(SCR\_repeat)                                                     | False  |  -0.509|
| DST      | Lost\_in\_tumor   | PF00435|Spectrin\_repeat                                                                 | False  |  -0.510|
| DST      | Gained\_in\_tumor | PF00681|Plectin\_repeat                                                                  | False  |  -0.510|
| ARHGAP5  | Lost\_in\_tumor   | PS51676|FF\_FF\_domain\_profile.                                                         | False  |  -0.516|
| ADAMTSL4 | Lost\_in\_tumor   | PF00090|Thrombospondin\_type\_1\_domain                                                  | False  |  -0.526|
| HIRIP3   | Lost\_in\_tumor   | PS50313|GLU\_RICH\_Glutamic\_acid-rich\_region\_profile.                                 | False  |  -0.541|
| ERCC8    | Lost\_in\_tumor   | PS50082|WD\_REPEATS\_2\_Trp-Asp\_(WD)\_repeats\_profile.                                 | False  |  -0.593|
| ERCC8    | Lost\_in\_tumor   | PF00400|WD\_domain,\_G-beta\_repeat                                                      | False  |  -0.593|
| ALCAM    | Lost\_in\_tumor   | PS50835|IG\_LIKE\_Ig-like\_domain\_profile.                                              | False  |  -0.612|
| BRF1     | Lost\_in\_tumor   | PF00382|Transcription\_factor\_TFIIB\_repeat                                             | False  |  -0.619|
| NR2F2    | Lost\_in\_tumor   | PS51030|NUCLEAR\_REC\_DBD\_2\_Nuclear\_hormone\_receptors\_DNA-binding\_domain\_profile. | False  |  -0.663|
| NR2F2    | Lost\_in\_tumor   | PF00105|Zinc\_finger,*C4\_type*(two\_domains)                                            | False  |  -0.663|
| FLT1     | Lost\_in\_tumor   | PF13895|Immunoglobulin\_domain                                                           | True   |  -0.679|
| RABGAP1L | Lost\_in\_tumor   | PF00566|Rab-GTPase-TBC\_domain                                                           | False  |  -0.689|
| SP140L   | Lost\_in\_tumor   | PS50016|ZF\_PHD\_2\_Zinc\_finger\_PHD-type\_profile.                                     | False  |  -0.701|
| SP140L   | Lost\_in\_tumor   | PF00628|PHD-finger                                                                       | False  |  -0.701|
| TAF5L    | Lost\_in\_tumor   | PS50082|WD\_REPEATS\_2\_Trp-Asp\_(WD)\_repeats\_profile.                                 | False  |  -0.710|
| TAF5L    | Lost\_in\_tumor   | PF00400|WD\_domain,\_G-beta\_repeat                                                      | False  |  -0.710|
| MYOM3    | Lost\_in\_tumor   | PF07679|Immunoglobulin\_I-set\_domain                                                    | False  |  -0.764|
| LAMA3    | Lost\_in\_tumor   | PF00053|Laminin\_EGF-like\_(Domains\_III\_and\_V)                                        | False  |  -0.819|
| ANK3     | Lost\_in\_tumor   | PF12796|Ankyrin\_repeats\_(3\_copies)                                                    | True   |  -0.828|
| COL12A1  | Gained\_in\_tumor | PF00041|Fibronectin\_type\_III\_domain                                                   | False  |  -0.970|
| ANKRD29  | Lost\_in\_tumor   | PS50088|ANK\_REPEAT\_Ankyrin\_repeat\_profile.                                           | False  |  -1.012|
| ANKRD29  | Lost\_in\_tumor   | PF00023|Ankyrin\_repeat                                                                  | False  |  -1.012|
| KLF3     | Lost\_in\_tumor   | PS50157|ZINC\_FINGER\_C2H2\_2\_Zinc\_finger\_C2H2\_type\_domain\_profile.                | False  |  -2.168|

![*meScore vs fold-change, measured as the ratio between the proportion of mutations affecting the feature and the expected proportion (relative size of the feature). In purple, significant cases before multitesting correction; in orange, significant cases after multitesting correction.*](20151006_mutation_analysis_files/figure-markdown_github/unnamed-chunk-4-1.png)

It seems that the genes that accumulate most mutations in a splicing gained/lost feature have meScores close to 0. In the case that there are mutations equivalent to splicing alterations (meScore \> 0), they seem to be less concentrated. This suggest that they are more specific, tumorogenic mutations, rather than random alterations.

Looking more in detail to the genes with meScre \> 0.7, we find:

-   EPS15: annotated as driver. Losing an EH domain, a [substrate for the tyrosine kinase activity of the epidermal growth factor receptor](http://www.ncbi.nlm.nih.gov/pubmed/11911876).

> To date, several EH-containing and EH-binding proteins have been identified, which establish in the cell a network of protein:protein interactions, defined as the EH network. This network coordinates cellular functions connected with endocytosis, actin remodeling and intracellular transduction of signals.

No I3D interaction is describing any change in this gene.

-   TDRD7
-   MCC: or Mutated In Colorectal Cancers...should be annotated as driver? ([candidate colorectal tumor suppressor gene that is thought to negatively regulate cell cycle progression](http://www.genecards.org/cgi-bin/carddisp.pl?gene=MCC)). The region might contain several phosphorylation sites.
-   CNTN4

General overlap in features affected by mutations and switches
==============================================================

We look for structural features that are frequently mutated and we compare them to the structural features that are usually switched. The interest is to look for splicing changes that are equivalent to mutation changes. We use the following parameters: \* Mutation frequency: proportion of all the functional mutations that affect one particular feature. Measured taking into account only the features in the most abundant isoform of each gene. \* Switch frequency: sum of all the patients that, per switch, have a feature affected. \* Domain frequency: frequency of a feature in the proteome, counting the features in the most abundant isoform for each gene. Used to normalize the previous quantities.

We annotated the domains relatively highly mutated and highly switched (highly meaning \> quantile 95).

View of the dat
---------------

### Analysis per tumor type

``` r
feat_enrich <- read.delim(paste0(wd,'feature_enrichment.txt'), header=TRUE)
feat_enrich <- feat_enrich[feat_enrich$DomainFrequency != 0,]
feat_enrich$MutTotal <- feat_enrich$MutIn + feat_enrich$MutOut
feat_enrich$SwitchesTotal <- feat_enrich$SwitchesIn + feat_enrich$SwitchesOut
feat_enrich$MutFreq <- feat_enrich$MutIn/feat_enrich$MutTotal
feat_enrich$SwitchFreq <- feat_enrich$SwitchesIn/feat_enrich$SwitchesTotal

for (cancer in cancerTypes){
  
  cat(paste0(cancer,"\n"))
  cancer.feat_enrich <- feat_enrich[feat_enrich$Cancer==cancer,]
  
  # plot normalized mutation frequency and switched frequency
  minMut <- quantile(cancer.feat_enrich$MutFreq/cancer.feat_enrich$DomainFrequency,0.95)
  minSwitch <- quantile(cancer.feat_enrich$SwitchFreq/cancer.feat_enrich$DomainFrequency,0.95)
  
  # plot raw data
  r <- cor(cancer.feat_enrich$MutFreq/cancer.feat_enrich$DomainFrequency,
           cancer.feat_enrich$SwitchFreq/cancer.feat_enrich$DomainFrequency)
  
  p <- ggplot(cancer.feat_enrich,aes(MutFreq/DomainFrequency,SwitchFreq/DomainFrequency)) + 
    geom_point() + 
    geom_point(data=subset(cancer.feat_enrich, MutFreq/DomainFrequency > minMut & SwitchFreq/DomainFrequency > minSwitch),aes(MutFreq/DomainFrequency,SwitchFreq/DomainFrequency,color=Domain)) + 
    smartas_theme() + 
    geom_smooth(method=lm) + 
    geom_text(x=20,y=100,label=paste0("R = ",round(r,2)))
  p <- direct.label(p)
  
  print(p)  
}
```

    ## brca

    ## Loading required package: proto

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-1.png)

    ## coad

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-2.png)

    ## hnsc

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-3.png)

    ## kich

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-4.png)

    ## kirc

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-5.png)

    ## kirp

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-6.png)

    ## lihc

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-7.png)

    ## luad

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-8.png)

    ## lusc

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-9.png)

    ## prad

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-10.png)

    ## thca

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_mutVsSwitch-11.png)

The non-logarithmic scale leads to think that usually changes achieved by splicing are not achieved by mutations and viceversa. Still, there are some outlier domains that are a bit outside the L-shape, suggesting that some features can either be affected by any mechanism. To find out about those, we remove the deatures that have 0 on either SwitchFreq or MutFreq, and perform a log2 to better observe them.

``` r
  # plot log representations and remove those that have 0 on either side
for (cancer in cancerTypes){
  
  cat(paste0(cancer,"\n"))
  cancer.feat_enrich <- feat_enrich[feat_enrich$Cancer==cancer,]
  
  # plot normalized mutation frequency and switched frequency
  minMut <- quantile(cancer.feat_enrich$MutFreq/cancer.feat_enrich$DomainFrequency,0.95)
  minSwitch <- quantile(cancer.feat_enrich$SwitchFreq/cancer.feat_enrich$DomainFrequency,0.95)

  # filter data and get correlation
  r.df <- subset(cancer.feat_enrich, MutFreq != 0 & SwitchFreq != 0)
  r <- cor(log2(r.df$MutFreq/r.df$DomainFrequency),
           log2(r.df$SwitchFreq/r.df$DomainFrequency))
  
  p <- ggplot(r.df,aes(log2(MutFreq/DomainFrequency),log2(SwitchFreq/DomainFrequency))) + 
    geom_point() + 
    geom_point(data=subset(r.df, MutFreq/DomainFrequency > minMut & SwitchFreq/DomainFrequency > minSwitch),aes(log2(MutFreq/DomainFrequency),log2(SwitchFreq/DomainFrequency),color=Domain)) + 
    smartas_theme() + 
    geom_smooth(data=subset(r.df,MutFreq>0 & SwitchFreq>0),method=lm) + 
    geom_text(x=-2,y=5,label=paste0("R = ",round(r,2)))
  p <- direct.label(p)
  
  print(p)
}
```

    ## brca

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-1.png)

    ## coad

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-2.png)

    ## hnsc

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-3.png)

    ## kich

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-4.png)

    ## kirc

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-5.png)

    ## kirp

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-6.png)

    ## lihc

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-7.png)

    ## luad

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-8.png)

    ## lusc

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-9.png)

    ## prad

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-10.png)

    ## thca

![](20151006_mutation_analysis_files/figure-markdown_github/tumor_log_mutVsSwitch-11.png)

When applying the log2 scale, not only some correlation appears, but some domains that appear to be interesting.

Analysis aggregating all tumors
-------------------------------

We apply the same process to aggregated data from all tumors.

``` r
feat_enrich_agg <- ddply(feat_enrich,.(Domain), summarise, 
                        MutIn=sum(MutIn), MutTotal=sum(MutTotal),
                        SwitchesIn=sum(SwitchesIn), SwitchesTotal=sum(SwitchesTotal),
                        MeanDomainFrequency=mean(DomainFrequency))
feat_enrich_agg$MutFreq <- feat_enrich_agg$MutIn/feat_enrich_agg$MutTotal
feat_enrich_agg$SwitchFreq <- feat_enrich_agg$SwitchesIn/feat_enrich_agg$SwitchesTotal

minMut <- quantile(feat_enrich_agg$MutFreq/feat_enrich_agg$MeanDomainFrequency,0.95)
minSwitch <- quantile(feat_enrich_agg$SwitchFreq/feat_enrich_agg$MeanDomainFrequency,0.95)

# plot original data
r <- cor(feat_enrich_agg$MutFreq/feat_enrich_agg$MeanDomainFrequency,
         feat_enrich_agg$SwitchFreq/feat_enrich_agg$MeanDomainFrequency)

p <- ggplot(data=feat_enrich_agg,
            aes(MutFreq/MeanDomainFrequency,SwitchFreq/MeanDomainFrequency)) + 
  geom_point() + 
  geom_point(data=subset(feat_enrich_agg,
                         MutFreq/MeanDomainFrequency > minMut & 
                         SwitchFreq/MeanDomainFrequency > minSwitch),
             aes(MutFreq/MeanDomainFrequency,
                 SwitchFreq/MeanDomainFrequency,
                 color=Domain)) + 
  smartas_theme() + 
  geom_smooth(method=lm) + 
  geom_text(x=20,y=100,label=paste0("R = ",round(r,2)))
p <- direct.label(p)

p
```

![](20151006_mutation_analysis_files/figure-markdown_github/aggregated_mutVsSwitch-1.png)

``` r
# plot log2 representations and remove those that have 0 on either side
r.df <- subset(feat_enrich_agg, MutFreq != 0 & SwitchFreq != 0)
r <- cor(log2(r.df$MutFreq/r.df$MeanDomainFrequency),
         log2(r.df$SwitchFreq/r.df$MeanDomainFrequency))

p <- ggplot(data=r.df,
            aes(log2(MutFreq/MeanDomainFrequency),log2(SwitchFreq/MeanDomainFrequency))) + 
  geom_point() + 
  geom_point(data=subset(r.df,
                         MutFreq/MeanDomainFrequency > minMut & 
                         SwitchFreq/MeanDomainFrequency > minSwitch),
             aes(log2(MutFreq/MeanDomainFrequency),
                 log2(SwitchFreq/MeanDomainFrequency),
                 color=Domain)) + 
  smartas_theme() + 
  geom_smooth(data=subset(r.df,MutFreq>0 & SwitchFreq>0),method=lm) + 
  geom_text(x=-3,y=2,label=paste0("R = ",round(r,2)))
p <- direct.label(p)

p
```

![](20151006_mutation_analysis_files/figure-markdown_github/aggregated_log_mutVsSwitch-1.png)

Similar results appear, but the correlation is quite lost, suggesting that the domains that need to be affected in each cancer type are tumor specific.

Cluster in domain function
--------------------------

We wanted to know if there are some common functions between
