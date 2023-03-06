###
## @author [Wankun Deng]
## @email [dengwankun@hotmail.com]
## @create date 2022-04-04 21:07:03
## @modify date 2022-04-04 21:07:03
## @desc [description]
###
library(ggplot2)
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

args=commandArgs(trailingOnly=T)
data_fn=args[1]
bin_len=as.numeric(args[2])
fo=args[3]
pdf(fo,height=4.5)

data<-read.table(paste0(data_fn),sep='\t',header=T,stringsAsFactors=T)
sep1=bin_len+0.5
sep2=bin_len*2+.5
p<-ggplot(data,aes(x=X,y=Y,color=Group,group=Group))+geom_line(size=1)+
    geom_vline(xintercept=c(sep1, sep2), linetype='dashed', size = 0.5) +
  annotate('text', label=c("5'UTR", "CDS", "3'UTR"), size = 7,x=c( sep1/2, (sep1+sep2)/2, sep2+sep1/2 ),
           y=max(data[,3])*0.9)  +
    ylab('Proportion of peaks')+xlab('Bin')+
    scale_colour_Publication()+ theme_Publication()+
    theme(panel.spacing = unit(0.5, 'lines'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=15),
        legend.position = c(0.17,0.5))+guides(color=guide_legend(ncol=1,title.position='top'))
print(p)

dev.off()