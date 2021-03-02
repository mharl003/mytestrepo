#!/usr/bin/env Rscript

library ('DESeq2')

#code to run
directory<-'/bigdata/jinkersonlab/shared/space/sandbox/RNASeq_template_Heinz/results/mutant_vs_wildtype_leaf'
sampleFiles<-grep('count', list.files(directory), value=TRUE)
sampleFiles
sampleCondition<-c('M', 'M', 'M', 'M', 'W', 'W', 'W', 'W', 'W')
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)

sampleTable

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)

ddsHTSeq

colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('M','W'))

dds<-DESeq(ddsHTSeq)

res<-results(dds)
res<-res[order(res$padj),]
head(res)

plotMA(dds,ylim=c(-3,3),main='DESeq2')
dev.copy(png,'volcano_diff_expression_W_vs_M.png')
dev.off()

mcols(res,use.names=TRUE)

write.csv(as.data.frame(res),file='SPACE_tomato_Mutant_vs_WT.csv')

summary(res)
