###
## @author [Wankun Deng]
## @email [dengwankun@hotmail.com]
## @create date 2022-04-04 21:06:54
## @modify date 2022-04-04 21:06:54
## @desc [description]
###
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)

args=commandArgs(trailingOnly=TRUE)
figure_out=args[1]

set1 <- read.csv(args[2],header=FALSE,sep="\t")[,1]
set2 <- read.csv(args[3],header=FALSE,sep="\t")[,1]
if(length(args)>3){
set3 <- read.csv(args[4],header=FALSE,sep="\t")[,1]
myCol <- brewer.pal(length(args)-1, "Pastel2")

venn.diagram(
        x = list(set1, set2, set3),
        category.names = c("Rep1" , "Rep2" , "Rep3"),
        filename = figure_out,
        output=TRUE,
        
        # Output features
        imagetype="tiff" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = myCol,
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        rotation = 1
)}else{
    venn.diagram(
        x = list(set1, set2),
        category.names = c("Rep1" , "Rep2"),
        filename = figure_out,
        output=TRUE,
        
        # Output features
        imagetype="tiff" ,
        height = 480 , 
        width = 480 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 2,
        lty = 'blank',
        # fill = myCol,
        
        # Numbers
        cex = .6,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 0.6,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        # cat.pos = c(-27, 27, 135),
        # cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans"#,rotation = 1
        )
}