###
## @author [Wankun Deng]
## @email [dengwankun@hotmail.com]
## @create date 2022-03-15 17:30:31
## @modify date 2022-03-15 17:30:31
## @desc [description]
###
library(readr)
library(GenomicRanges)
library(BisRNA)
library(idr)
library(scales)

idr_one_rbp<-function(peak1, peak2, tissue)
{
    print(tissue)
    file1<-paste0(peak1)
    file2<-paste0(peak2)

    df1 <- read_delim(file1, delim="\t", col_names=FALSE)
    df2 <- read_delim(file2, delim="\t", col_names=FALSE)
    peak1 <- GRanges(df1$X1, IRanges(df1$X2, df1$X3), score=df1$X8,strand=df1$X6)
    peak2 <- GRanges(df2$X1, IRanges(df2$X2, df2$X3), score=df2$X8,strand=df2$X6)
    peak1 <- keepStandardChromosomes(peak1, pruning.mode="coarse")
    peak2 <- keepStandardChromosomes(peak2, pruning.mode="coarse")
    fo <- findOverlaps(peak1, peak2,minoverlap=2)
    fo <- as.data.frame(fo)
    fo <- fo[!duplicated(fo$queryHits) & !duplicated(fo$subjectHits),]

    # use -log10(pvals) as score, for peaks with pvals==0, score=max(-log10(pval)) + rescaled signal value 
    y1 <- peak1$score[fo[,1]]
    y1_fc<-rescale((df1$X7[fo[,1]])[y1==0])
    y1<- -log10(y1)
    y1[!is.finite(y1)]<-max(y1[is.finite(y1)])+y1_fc

    y2 <- peak2$score[fo[,2]]
    y2_fc<-rescale((df2$X7[fo[,2]])[y2==0])
    y2<- -log10(y2)
    y2[!is.finite(y2)]<-max(y2[is.finite(y2)])+y2_fc

    dat <- cbind(y1, y2)

    res <- est.IDR(dat, mu=3, sigma=1, rho=.9, p=.5)
    df1_1<-df1[fo[res$idr<0.01,1],]
    df2_2<-df2[fo[res$idr<0.01,2],]

    df3<-df1_1
    df3$X8<-rowMeans(data.frame(P1=df1_1$X7,P2=df2_2$X7))
    df3$X8<-apply(data.frame(P1=df1_1$X8,P2=df2_2$X8),1,fisher.method)
    df3$X9<-res$idr[res$idr<0.01]
    write.table(df3,file=paste0('idr_peak/',tissue,'_IDR.bed'),sep='\t',row.names = F,col.names=F,quote=F)
}
args = commandArgs(trailingOnly=TRUE)
idr_one_rbp(args[1],args[2],args[3])