###
## @author [Wankun Deng]
## @email [dengwankun@gmail.com]
## @create date 2023-03-06 01:13:15
## @modify date 2023-03-06 01:13:15
## @desc [description]
###
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(colorspace))
suppressMessages(library(circlize))
suppressMessages(library(Hmisc))
suppressMessages(library(scales))
suppressMessages(library(Rtsne))
suppressMessages(library(corrplot))
suppressMessages(library(preprocessCore))
suppressMessages(library(locfdr))
suppressMessages(library(dendextend))
suppressMessages(library(ggpubr))
theme_set(theme_grey(base_size=20))
suppressMessages(library(dplyr))
suppressMessages(library(MASS))
suppressMessages(library(reshape2))
suppressMessages(library(ggseqlogo))
suppressMessages(library(stringr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(networkD3))
suppressMessages(library(htmlwidgets))
suppressMessages(library(webshot))
suppressMessages(library(ggseqlogo))
suppressMessages(library(rlist))
suppressMessages(library(WGCNA))
suppressMessages(library(RUVSeq))
suppressMessages(library(edgeR))
suppressMessages(library(EDASeq))
suppressMessages(library(tximport))
suppressMessages(library(RColorBrewer))
suppressMessages(library(DESeq2))

# suppressMessages(library(psycModel))
set.seed(42)
ht_opt$message = FALSE
## code for word cloud
#source("https://gist.githubusercontent.com/jokergoo/bfb115200df256eeacb7af302d4e508e/raw/14f315c7418f3458d932ad749850fd515dec413b/word_cloud_grob.R")
## code for word cloud saved to local
source("scripts/word_cloud.R")
`%notin%` <- Negate(`%in%`)
font_size<-16
cv<-function(vec){
  vec<-vec[!is.na(vec)]
  sd(vec)/mean(vec)
}

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

theme_Publication <- function(base_size=20, base_family="") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.5, "cm"),
               legend.margin = unit(0, "cm"),
#                legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506",
                                                                "#a6cee3","#fb9a99","#984ea3","#ffff33",'#6060f4','#ad27ad')), ...)
}

scale_fill_Publication_continuous <- function(...){
      library(scales)
      continuous_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506",
                                                                "#a6cee3","#fb9a99","#984ea3","#ffff33",'#6060f4','#ad27ad')), ...)
}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506",
                                                                  "#a6cee3","#fb9a99","#984ea3","#ffff33",'#6060f4','#ad27ad')), ...)
}
tissue_name<-setNames(
    c('CirculatorySystem','CirculatorySystem','CirculatorySystem','CirculatorySystem','CirculatorySystem',
      'CirculatorySystem','CirculatorySystem','CirculatorySystem','DigestiveSystem','DigestiveSystem','DigestiveSystem',
      'DigestiveSystem','DigestiveSystem','DigestiveSystem','DigestiveSystem','DigestiveSystem','DigestiveSystem',
      'DigestiveSystem','DigestiveSystem','DigestiveSystem','DigestiveSystem','GlandularTissues','GlandularTissues',
      'GlandularTissues','GlandularTissues','GlandularTissues','GlandularTissues','GlandularTissues','GlandularTissues',
      'GlandularTissues','GlandularTissues','GlandularTissues','GlandularTissues','ImmuneSystem','ImmuneSystem',
      'NervousSystem','NervousSystem','NervousSystem','NervousSystem','NervousSystem','NervousSystem','NervousSystem',
      'NervousSystem','NervousSystem','NervousSystem','NervousSystem','NervousSystem','NervousSystem','NervousSystem',
      'Other','Other','Other','Other','Other','ReproductiveSystem','ReproductiveSystem','ReproductiveSystem',
      'ReproductiveSystem','ReproductiveSystem','ReproductiveSystem','RespiratorySystem','RespiratorySystem',
      'RespiratorySystem','UrinarySystem','UrinarySystem','UrinarySystem','Control'),
    c('AAO','AAD','AAS','AHE','ALV','APE','ARV','FHE','AAP','ACA','ACL','ADU','AES','AIC','AIM','AJE','ASI','AST','ACS',
      'ASO','ATG','AAC','AAG','AMG','APA','APG','ASG','ASP','ATY','ATN','FAG','FSP','FTH','ALN','ALE','ABR','ACN','ACE',
      'ACC','ACO','AFL','AHI','AMO','APO','ASC','ATL','ATA','FBR','FSC','AAT','ALI','ASM','FWH','FLI','AEP','AOV','APL',
      'APR','ATE','AUT','ALU','ATR','FLU','ABL','AKI','FKI','Control')
    )

gm_mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

ClusterZ <- function(t, file, TF1, TF2, start, end){

# t: distance that we want to test a clustering pattern between two TF bindings
# file: the data set with two TF binding sites
# TF1: name of the first TF
# TF2: name of the second TF
# start: start location to TSS of upstream region that is used in data
# end: end location to TSS of upstream region


TF <- read.table(file, header=T, sep="")

ind <- TF$count1*TF$count2
TF <- subset(TF, ind!=0)

n1 <- grep(TF1,names(TF)); l1 <- length(n1) 
n2 <- grep(TF2,names(TF)); l2 <- length(n2)
L <- end-start


# Kernel Density Estimation #

 kernel1<-c() 
 for (i in 1:l1) {kernel1 <-c(kernel1 ,na.omit(TF[,n1[i]]))}

 kernel2<-c() 
 for (i in 1:l2) {kernel2 <-c(kernel2,na.omit(TF[,n2[i]]))}


# Find a bin size #

 y <- seq_len(L)
 bin0 <- y[ (L/100)%%y == 0 ]*100 # possible bin sizes #
 
 i <- 0; prodn <- 0
 while(prodn <= 0){
  i <- i+1
  bin <- bin0[i]
  grid <- seq(start,end,bin)
  TF1_ni <- table(cut(kernel1, grid))
  TF2_ni <- table(cut(kernel2, grid))
  prodn <- prod(TF1_ni)*prod(TF2_ni)
 }


# lambda(x) function #

lambda1 <- function(x){
	lc <- c()
	for (i in 1:length(x)){
	if (start >= x[i]) {lc[i]<-(1/(bin*length(kernel1)))*TF1_ni[cut(start+1,grid)]}
	if (x[i] > start) {lc[i]<-(1/(bin*length(kernel1)))*TF1_ni[cut(x[i],grid)]}
	}
	return(lc)
}
 
lambda2 <- function(x){
	lc <- c()
	for (i in 1:length(x)){
	if (start >= x[i]) {lc[i]<-(1/(bin*length(kernel2)))*TF2_ni[cut(start+1,grid)]}
	if (x[i] > start) {lc[i]<-(1/(bin*length(kernel2)))*TF2_ni[cut(x[i],grid)]}
	}
	return(lc)
}



# bivariate K-function # 
  
bivk <- function(t,g1,g2){
    ww <- 0
	for (i in 1:length(g1)){
		for (j in 1:length(g2)){
			d <- abs(g1[i]-g2[j])
			if (d < t){ 
				if ((g1[i]-d) >= start & (g1[i]+d) <= end) w <- 1
    			if ((g1[i]-d) < start | (g1[i]+d) > end) w <- 2 
				ww <- ww + (w/(lambda1(g1[i])*lambda2(g2[j])))
			    }
			}
		}
	return(ww/(10000*length(g1)*length(g2)))
}




# Variacne of K(t) #

seq2 <- grid[-1]
f1 <- function(x) 1/lambda1(x) 
cf1<-c(); for (i in 1:length(seq2)){cf1[i]<-f1(seq2[i])} 
f2 <- function(x) 1/lambda2(x) 
cf2<-c(); for (i in 1:length(seq2)){cf2[i]<-f2(seq2[i])} 



vark <- function(t){
	I1 <- 0; I2 <-0; I3 <-0
	x1 <- seq(start, start+t)[-1]
	y1.1 <- function(x) return(seq(start,-start+2*x)[-1])
	y1.2 <- function(x) return(seq(-start+2*x, x+t)[-1])
	for (i in 1:length(x1)){ I1 <- I1 + f1(x1[i])*(sum(table(cut(y1.1(x1[i]),grid))*cf2) + sum(table(cut(y1.2(x1[i]),grid))*cf2*4))}
	x2 <- seq(start+t, end-t)[-1]
	y2 <- function(x) return(seq(x-t, x+t)[-1])
	for (i in 1:length(x2)){ I2 <- I2 + f1(x2[i])*sum(table(cut(y2(x2[i]),grid))*cf2) }
	x3 <- seq(end-t, end)[-1]
	y3.1 <- function(x) return(seq(2*x-end,end)[-1])
	y3.2 <- function(x) return(seq(x-t, 2*x-end)[-1])
	for (i in 1:length(x3)){ I3 <- I3 + f1(x3[i])*(sum(table(cut(y3.1(x3[i]),grid))*cf2) + sum(table(cut(y3.2(x3[i]),grid))*cf2*4))}
	return(I1+I2+I3-(2*L*t)^2)
}


varkt <- vark(t)
Ixt <- (2*t)^2*(sum(cf1*bin)-10000^2)
Iyt <- (2*t)^2*(sum(cf2*bin)-10000^2)


var_k <- function(nx,ny){
	(varkt+(ny-1)*Ixt+(nx-1)*Iyt)/(100000000*nx*ny)
	}


# Calculating Zc and its p-value # 

kt <- rep(0,nrow(TF))
vpm <- rep(0,nrow(TF))

for (i in 1:nrow(TF)){
	   g1 <- na.omit(as.numeric(c(TF[i,(n1[1]:n1[length(n1)])])))
	   g2 <- na.omit(as.numeric(c(TF[i,(n2[1]:n2[length(n2)])])))
       nx <- length(g1); ny <- length(g2)
       kt[i] <- bivk(t,g1,g2)
       vpm[i] <- var_k(nx,ny)
}



Zc <-sqrt(nrow(TF))*mean((kt-2*t)/sqrt(vpm)) 
pvalue <- 2*(1-pnorm(abs(Zc)))

return(list(Zc=Zc, pvalue=pvalue))

}