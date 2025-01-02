# DESeq2 워크플로우

## Load package
```R
suppressMessages({
    library("DESeq2")
    library(pheatmap)
    library(withr)
    #library(tidyverse)
    library(RColorBrewer)
    library(gplots)
    library(dplyr)
    })
```

## Set path
```R
setwd("/data-blog/bi1")
getwd()
```
```plain text
'/data-blog/bi1'
```

## Run DESeq2
```r
dir <- 'res'
S1 <- '33'
S2 <- '150'

countdata <- read.csv("results.csv", header=TRUE, sep=',')
colnames(countdata) <- c('GeneID','150-1','150-2','150-3','33-1','33-2','33-3','con-1','con-2','con-3')
countdata <- countdata[, c('GeneID','150-1','150-2','150-3','33-1','33-2','33-3','con-1','con-2','con-3')]

selected_columns <- paste(c('GeneID',paste0(S2,"-1"), paste0(S2,"-2"), paste0(S2,"-3"),paste0(S1,"-1"), paste0(S1,"-2"), paste0(S1,"-3")), sep="")
countdata <- countdata[, selected_columns]
countdata <- countdata[rowSums(countdata[, -1]) != 0, ]

sample.names <- paste(c(paste0(S2,"-1"), paste0(S2,"-2"), paste0(S2,"-3"),paste0(S1,"-1"), paste0(S1,"-2"), paste0(S1,"-3")), sep="")
conditions <- factor(c(S2,S2,S2,S1,S1,S1))
metadata <- data.frame(Sample = sample.names, group = conditions)

metadata

N <- dim(countdata)[[2]]
cData = countdata[,2:N]
GeneID = countdata[,1]
rownames(cData) = GeneID

dds <- DESeqDataSetFromMatrix(countData = cData,
                                colData = metadata,
                                design = ~group)
dds$group <- relevel(dds$group, ref = S1)
colData(dds)

dds <- DESeq(dds)  
vsd <- vst(dds, blind=FALSE)
rld <- rlogTransformation(dds, blind=FALSE)

res <- results(dds, contrast = c("group", S2, S1))
res_tbl <- as.data.frame(res)
res_tbl$GeneID <- rownames(res_tbl)
res_tbl <- res_tbl[order(res_tbl$padj), ]

NM_no_NA <- na.omit(res)
res_cut <- NM_no_NA[NM_no_NA$padj<0.05,]
res_cut

# padj
val_str <- 'padj'
cutoff <- 0.05
cutoff_str <- as.character(cutoff)

sig_res <- dplyr::filter(res_tbl, padj < cutoff)
sig_res <- dplyr::arrange(sig_res, padj)
sig_res_file <- paste0(dir, '/res_', S2, '_', S1, '_', val_str, cutoff_str, '.csv')
#write.csv(sig_res, file = sig_res_file)

print(paste0(S2, ' vs ', S1, ' | ', val_str, '<', cutoff_str))
sig_idx <- res$padj < cutoff & !is.na(res$padj)
sig_dat <- assay(rld)[sig_idx, ]
annC <- data.frame(condition = conditions)
rownames(annC) <- colnames(sig_dat)
heat_colors <- brewer.pal(6, "RdYlGn")
heat_colors_reversed <- rev(heat_colors)
ann_colors <- list(condition = setNames(c("#F7819F", "#58D3F7"), c(S2, S1)))
```
```plain text
A data.frame: 6 x 2
Sample	group
<chr>	<fct>
150-1	150
150-2	150
150-3	150
33-1	33
33-2	33
33-3	33
DataFrame with 6 rows and 2 columns
           Sample    group
      <character> <factor>
150-1       150-1      150
150-2       150-2      150
150-3       150-3      150
33-1         33-1      33 
33-2         33-2      33 
33-3         33-3      33

estimating size factors

estimating dispersions

gene-wise dispersion estimates

mean-dispersion relationship

final dispersion estimates

fitting model and testing

log2 fold change (MLE): group 150 vs 33 
Wald test p-value: group 150 vs 33 
DataFrame with 205 rows and 6 columns
         baseMean log2FoldChange     lfcSE      stat      pvalue        padj
        <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
ABHD2      50.721       1.352060  0.431587   3.13276 1.73168e-03 4.00143e-02
ADAM12    706.120      -0.571960  0.168494  -3.39454 6.87431e-04 2.03489e-02
ADD2     1819.643       0.868228  0.148791   5.83521 5.37230e-09 9.73246e-07
AIF1L     144.513       1.168923  0.283764   4.11935 3.79938e-05 2.07318e-03
AKAP5    1042.005      -0.637445  0.202189  -3.15271 1.61761e-03 3.81572e-02
...           ...            ...       ...       ...         ...         ...
ZNF655   774.2459      -0.910286  0.198632  -4.58277 4.58855e-06 3.52229e-04
ZNF682    59.7573      -1.382049  0.438336  -3.15295 1.61632e-03 3.81572e-02
ZNF77     76.0271      -1.231188  0.388382  -3.17004 1.52417e-03 3.71126e-02
ZRANB3   536.2301      -0.878732  0.179932  -4.88367 1.04128e-06 9.22422e-05
ZSCAN25 1257.3596      -0.460023  0.149161  -3.08408 2.04184e-03 4.57797e-02
[1] "150 vs 33 | padj<0.05"
```
