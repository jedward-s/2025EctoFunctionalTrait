#------------------start--------------------------------


library(ggplot2)
library(plyr)
library(gridExtra)
library(car)
library(visreg)
library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)
library(reshape2)
library(pbkrtest)
library(MuMIn)
library(ade4)
library(nlme)
library(car)
library(MuMIn)
library(MASS)
library(ape)
library(egg)
library(svglite)
library(dplyr)
library(hillR)
library(ggExtra)
library(ggside)
library(ggrepel)
library(tidyverse)

mytheme <-   theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
                   panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
                   panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
                   panel.border = element_rect(colour = "black", fill=NA, size=1),
                   axis.line = element_line(colour = 'black', linewidth = 0.5),#sets the axis line size
                   axis.ticks=element_line(colour = 'black', linewidth = 0.5), #sets the tick lines
                   axis.title.x = element_text( face="bold", size=14, color="black"), #size of x-axis title
                   axis.title.y = element_text(face="bold", size=14, color="black"), #size of y-axis title
                   axis.text.x = element_text(size=10, color="black"), #size of x-axis text
                   axis.text.y = element_text( size=10, color="black"),
                   legend.position = "none")#size of y-axis text


#--------readdata--------
meta <- read.csv("ECM_meta.csv", row.names = 1)

soil <- read.csv("ECM_soil.csv", row.names = 1)

#calculate SOM index
soil.st <- scale(soil) - min(scale(soil))

ddats <- as.matrix(vegdist(soil.st, method = "euclid"))
soil.db <- dbrda(ddats ~1 + Condition(as.matrix(meta[,c("PointX", "PointY")])), method="euclid", add = TRUE)

meta$soil1 <- scores(soil.db, display="sites")[,1]
meta$soil2 <- scores(soil.db, display="sites")[,2]




fungi <- read.csv("Fungi.csv", header = TRUE, row.names = 1)
 

#calculate ECM fungal alpha diversity

ecm <- t(subset(fungi, Lifestyle == "Ectomycorrhizal")[,0:75])

meta$ecmq1 <- apply(ecm,1,hill_taxa,q=1)

meta$ox <- ((meta$PO + meta$PX))

#-----------Figure 2--------------

mytheme1 <-   theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
                    panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
                    panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
                    #panel.border=element_blank(), #gets rid of square going around the entire graph
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line = element_line(colour = 'black', linewidth = 0.5),#sets the axis line size
                    axis.ticks=element_line(colour = 'black', linewidth = 0.5), #sets the tick lines
                    axis.title.x = element_text( face="bold", size=14, color="black"), #size of x-axis title
                    axis.title.y = element_text(face="bold", size=14, color="black"), #size of y-axis title
                    axis.text.x = element_text(size=7, color="black"), #size of x-axis text
                    axis.text.y = element_text( size=7, color="black"),
                    plot.margin = unit(c(-.5,-.5,-.5,-.5), "cm"),
                    legend.position = "none")#size of y-axis text

mytheme2 <-   theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
                    panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
                    panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
                    #panel.border=element_blank(), #gets rid of square going around the entire graph
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line = element_line(colour = 'black', linewidth = 0.5),#sets the axis line size
                    axis.ticks=element_line(colour = 'black', linewidth = 0.5), #sets the tick lines
                    axis.title.x = element_text( face="bold", size=14, color="black"), #size of x-axis title
                    axis.title.y = element_text(face="bold", size=14, color="black"), #size of y-axis title
                    axis.text.x = element_text(size=7, color="black"), #size of x-axis text
                    axis.text.y = element_text( size=7, color="black"),
                    plot.margin = unit(c(0,-.5,-.5,-.5), "cm"),
                    legend.position = "none")#size of y-axis text


mytheme3 <-   theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
                    panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
                    panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
                    #panel.border=element_blank(), #gets rid of square going around the entire graph
                    panel.border = element_rect(colour = "black", fill=NA, size=1),
                    axis.line = element_line(colour = 'black', linewidth = 0.5),#sets the axis line size
                    axis.ticks=element_line(colour = 'black', linewidth = 0.5), #sets the tick lines
                    axis.title.x = element_text( face="bold", size=14, color="black"), #size of x-axis title
                    axis.title.y = element_text(face="bold", size=14, color="black"), #size of y-axis title
                    axis.text.x = element_text(size=6.5, color="black"), #size of x-axis text
                    axis.text.y = element_text( size=7, color="black"),
                    plot.margin = unit(c(0,0,0,0), "cm"),
                    legend.position = "none")#size of y-axis text

meta.fun <- meta

meta.fun$G <- ordered(meta.fun$G, levels = c("Q", "C", "T"))


LCN <- ggplot(meta.fun, aes(x = G, y = leaf_CN, color = G,  fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme2  

LM <- ggplot(meta.fun, aes(x = MAOM_d13C, y = leaf_CN)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_x_continuous(limits = c(-27.01, -24.5), breaks = c(-27, -26, -25), expand = c(0,0))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

LLL <- ggplot(meta.fun, aes(x = leaf_d15N, y = leaf_CN)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

LLecm <- ggplot(meta.fun, aes(x = ecmq1, y = leaf_CN)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

LLox <- ggplot(meta.fun, aes(x = ox, y = leaf_CN)) +
  geom_point(alpha = 0.8, size =2,  col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1


LLs <- ggplot(meta.fun, aes(x = soil1, y = leaf_CN)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1


#MAOM_d13C

MM <- ggplot(meta.fun, aes(x = G, y = MAOM_d13C, color = G,  fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_y_continuous(limits = c(-27.01, -24.5), breaks = c(-27, -26, -25), expand = c(0,0))+
  ylab("")+ 
  xlab("")+ 
  mytheme2

ML <- ggplot(meta.fun, aes(x = leaf_d15N, y = MAOM_d13C)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_y_continuous(limits = c(-27.01, -24.5), breaks = c(-27, -26, -25), expand = c(0,0))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

Mecm <- ggplot(meta.fun, aes(x = ecmq1, y = MAOM_d13C)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_y_continuous(limits = c(-27.01, -24.5), breaks = c(-27, -26, -25), expand = c(0,0))+
  ylab("")+ 
  xlab("")+ 
  mytheme1


Mox <- ggplot(meta.fun, aes(x = ox, y = MAOM_d13C)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_y_continuous(limits = c(-27.01, -24.5), breaks = c(-27, -26, -25), expand = c(0,0))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

Ms <- ggplot(meta.fun, aes(x = soil1, y = MAOM_d13C)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_y_continuous(limits = c(-27.01, -24.5), breaks = c(-27, -26, -25), expand = c(0,0))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

#leaf_d15N


LL <- ggplot(meta.fun, aes(x = G, y = leaf_d15N, color = G, fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme2


Lecm <- ggplot(meta.fun, aes(x = ecmq1, y = leaf_d15N)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

Lox <- ggplot(meta.fun, aes(x = ox, y = leaf_d15N)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

Ls <- ggplot(meta.fun, aes(x = soil1, y = leaf_d15N)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1


#ecmq diversity

ECM <- ggplot(meta.fun, aes(x = G, y = ecmq1, color = G,fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme2

Eox <- ggplot(meta.fun, aes(x = (ox), y = ecmq1)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1


Es <- ggplot(meta.fun, aes(x = soil1, y = ecmq1)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = "lm", col = "black", se = FALSE)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1


#MBC

Ox <- ggplot(meta.fun, aes(x = G, y = (ox), color = G, fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme2

Os <- ggplot(meta.fun, aes(x = soil1, y = ox)) +
  geom_point(alpha = 0.8, size =2, col="black", shape = 21, aes(fill = G)) +
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1


S1 <- ggplot(meta.fun, aes(x = G, y = soil1, color = G,  fill = G)) +
  geom_boxplot(linewidth = 1, width = .9, alpha = 0.8, fill = NA)+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme1

blanktheme <-   theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
                      panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
                      panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
                      #panel.border=element_blank(), #gets rid of square going around the entire graph
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      axis.line = element_line(colour = 'black', linewidth = 0.5),#sets the axis line size
                      axis.ticks=element_line(colour = 'black', linewidth = 0.5), #sets the tick lines
                      axis.title.x = element_text( face="bold", size=14, color="black"), #size of x-axis title
                      axis.title.y = element_text(face="bold", size=14, color="black"), #size of y-axis title
                      axis.text.x = element_text(size=7, color="black"), #size of x-axis text
                      axis.text.y = element_text( size=7, color="black"),
                      plot.margin = unit(c(0,-.5,0,-.5), "cm"),
                      legend.position = "none")#size of y-axis text




blank <- ggplot() +
  geom_blank(
    mapping = NULL,
    data = NULL,
    stat = "identity",
    position = "identity",
    show.legend = NA,
    inherit.aes = TRUE
  )+
  blanktheme


library(gridExtra)
library(grid)


quartz()
Fig2 <- ggarrange(LCN, LM, LLL, LLecm, LLox, LLs, 
                  blank, MM, ML, Mecm, Mox, Ms, 
                  blank, blank,LL, Lecm, Lox, Ls,
                  blank, blank, blank, ECM, Eox, Es,
                  blank, blank, blank, blank,Ox, Os,
                  blank, blank, blank,blank, blank, S1, ncol = 6, nrow = 6)





#------ Figure 3-----------

isoC <- meta[,c("G","leaf_d13C", "MAOM_d13C", "bulk_d13C", "POM_d13C")]
isoN <- meta[,c("G","leaf_d15N", "MAOM_d15N", "bulk_d15N", "POM_d15N")]


isoC$G <- ordered(isoC$G, levels = c("Q", "C", "T"))
isoN$G <- ordered(isoN$G, levels = c("Q", "C", "T"))

isoC.melt <- melt(isoC)
isoN.melt <- melt(isoN)


isoC.melt$G <- ordered(isoC.melt$G, levels = c("Q", "C", "T"))

C <- ggplot(isoC.melt, aes(x = variable, y = value, color = G, fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  #geom_point(alpha = 0.3, size =3, col="black", shape = 21, position=position_jitter()) +
  #geom_text(aes(label=ID), col = "black",hjust=1, vjust=1)+
  #geom_point(shape=1, size = 3.25, colour = "black")+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  #scale_y_continuous(limits=c(0.075, 0.30))+
  #scale_y_log10()+
  ylab("")+ 
  xlab("")+ 
  mytheme

isoN.melt$G <- ordered(isoN.melt$G, levels = c("Q", "C", "T"))
N <- ggplot(isoN.melt, aes(x = variable, y = value, color = G, fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  #geom_point(alpha = 0.3, size =3, col="black", shape = 21, position=position_jitter()) +
  #geom_text(aes(label=ID), col = "black",hjust=1, vjust=1)+
  #geom_point(shape=1, size = 3.25, colour = "black")+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  #scale_y_continuous(limits=c(0.075, 0.30))+
  #scale_y_log10()+
  ylab("")+ 
  xlab("")+ 
  mytheme

isoC$MLC <- isoC$MAOM_d13C - isoC$leaf_d13C
isoN$MLN <- isoN$MAOM_d15N - isoN$leaf_d15N
quartz()
MLC <- ggplot(isoC, aes(x = G, y = MLC, color = G, fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  #geom_point(alpha = 0.3, size =3, col="black", shape = 21, position=position_jitter()) +
  #geom_text(aes(label=ID), col = "black",hjust=1, vjust=1)+
  #geom_point(shape=1, size = 3.25, colour = "black")+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  #scale_y_continuous(limits=c(0.075, 0.30))+
  #scale_y_log10()+
  ylab("")+ 
  xlab("")+ 
  mytheme
MLN <-  ggplot(isoN, aes(x = G, y = MLN, color = G, fill = G)) +
  geom_boxplot(size = 1, width = .9, alpha = 0.8, fill = NA)+
  #geom_point(alpha = 0.3, size =3, col="black", shape = 21, position=position_jitter()) +
  #geom_text(aes(label=ID), col = "black",hjust=1, vjust=1)+
  #geom_point(shape=1, size = 3.25, colour = "black")+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  #scale_y_continuous(limits=c(0.075, 0.30))+
  #scale_y_log10()+
  ylab("")+ 
  xlab("")+ 
  mytheme


ml.melt <- melt(data.frame(isoC[,c("G", "MLC")], isoN$MLN))

MLCN <-  ggplot(ml.melt, aes(x = variable, y = value, color = G, fill = G)) +
  geom_boxplot(size = 1, width = .8, alpha = 0.8, fill = NA)+
  #geom_point(alpha = 0.3, size =3, col="black", shape = 21, position=position_jitter()) +
  #geom_text(aes(label=ID), col = "black",hjust=1, vjust=1)+
  #geom_point(shape=1, size = 3.25, colour = "black")+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  #scale_y_continuous(limits=c(0.075, 0.30))+
  #scale_y_log10()+
  ylab("")+ 
  xlab("")+ 
  mytheme

MLCN.corr <-  ggplot(isoC, aes(x = MLC, y = isoN$MLN)) +
  geom_point(alpha = 0.8, size =3, col="black", shape = 21, aes(fill = G)) +
  stat_smooth(method = lm , se = FALSE, col = "black")+
  scale_color_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  scale_fill_manual(values=c("forestgreen","goldenrod",  "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme


Fig3 <- grid.arrange(                      
  arrangeGrob(C,  N,  ncol = 1),
  arrangeGrob(MLCN,  MLCN.corr,  ncol = 1),
  widths = c(2,1.5), 
  nrow = 1)  

ggsave(Figure3, filename= "Figure3.svg")



#------------ Figure 4-------------


soil.st <- (data.frame(scale(soil) - min(scale(soil))))

soil.db <- capscale(soil.st~1 + Condition(as.matrix(meta[,c("PointX", "PointY")])), method="euclid", add = TRUE, binary = FALSE)

summary(soil.db)

soil.sp <- data.frame(scores(soil.db, display = "species"))
soil.sp$Env <- row.names(soil.sp)


SOM <- ggplot() +
  geom_point(data = meta.fun, aes(x =soil1, y = soil2, fill = G, col = G), size =3, alpha = 0.3) +
  geom_point(data = meta.fun, aes(x =soil1, y = soil2), size =3, col="black", alpha = 0.4, shape = 21) +
  stat_ellipse(data = meta.fun, aes(x =soil1, y = soil2, col = G), alpha = 0.5)+
  scale_color_manual(values=c("forestgreen", "goldenrod", "royalblue"))+
  scale_fill_manual(values=c("forestgreen", "goldenrod", "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  scale_x_continuous(expand = c(0.15,0.1))+
  geom_text_repel(data = soil.sp, aes(x = (MDS1*2), 
                                   y = (MDS2*2), 
                                   label = Env), size = 2.5, box.padding = 0)+
  mytheme

Fig4 <- ggMarginal(SOM, groupColour = TRUE, groupFill = TRUE, type = "boxplot")





#--------- Figure 5--------

#Average fungal abundance at genus level to account for overly variable ASV ID's

fun.tax.in <- fungi[,76:ncol(fungi)]

fun.tax.in$genus.c <- gsub(";s__.*$", "", fun.tax.in$taxonomy)

fungi_asv_abun <- fungi[,0:75]

fungi_asv_abun$genus <- fun.tax.in[rownames(fungi_asv_abun), "genus.c"]
fun_genus_sum <- aggregate(. ~ genus, data=fungi_asv_abun, FUN=sum)

fun.gen <- t(fun_genus_sum[,c(-1)])

fun.rel <- decostand(t(fun_genus_sum[,-1] ), method = "hellinger")


ddatf <- as.matrix(vegdist(fun.rel, method = "jaccard"))

vars <- meta.fun[,c("leaf_CN", "MAOM_d13C", "leaf_d15N", "ecmq1", "ox", "soil1")]


modelfg <- dbrda(ddatf ~ G + Condition(as.matrix(meta.fun[,c("PointX", "PointY")])), data = meta.fun, method="jaccard", add = TRUE)
modelfv <- dbrda(ddatf ~ . + Condition(as.matrix(meta.fun[,c("PointX", "PointY")])), data = vars, method="jaccard", add = TRUE)
modelf <- dbrda(ddatf ~1 + Condition(as.matrix(meta.fun[,c("PointX", "PointY")])), data = meta.fun, method="jaccard", add = TRUE)

anova(modelfg, by = "margin", type = 3)
anova(modelfv, by = "margin", type = 3)

summary(modelfg)
summary(modelfv)
summary(modelf)

ordistep(modelfv)                                                                                    

modelfv <- dbrda(formula = ddatf ~ MAOM_d13C + leaf_d15N + ecmq1 + soil1 + Condition(as.matrix(meta.fun[,c("PointX", "PointY")])), data = vars, add = TRUE, method = "jaccard")

scf <- data.frame(scores(modelf))



fo <- ggplot(meta, aes(x = scf[,1], y = scf[,2])) +
  #geom_polygon(data = hulls, alpha = 0.1, size= 3, fill=NA) +
  geom_point(size =5, alpha = 0.5, aes(color = meta.fun$G,  fill = meta.fun$G)) +
  #geom_segment(data=ecm.vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
  #             arrow = arrow(length = unit(0.5, "cm")),colour="black", size = 1.3, alpha = 0.5) + 
  stat_ellipse(aes(col = G), size = 3, alpha = 0.8)+
  #geom_text(data=ecm.vec.sp.df, aes(x=MDS1,y=MDS2,label=Env),size=4, nudge_y = -0.02, nudge_x = -0.02)+
  #geom_point(shape=1, size = 3.25, colour = "black")+
  scale_color_manual(values=c("goldenrod", "forestgreen", "royalblue"))+
  scale_fill_manual(values=c("goldenrod", "forestgreen", "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme 



scbd <- data.frame(scores(modelfv, display = "sites"))
scbp <- data.frame(scores(modelfv, display = "bp"))


d <- t(data.frame("MAOM δC**", "Leaf δN*", "EMF diversity**", "Soil index*"))

scbp$Env <- d[,1]



quartz()
fb <- ggplot(meta, aes(x = scbd[,1], y = scbd[,2])) +
  geom_point(size =5, aes(color = G,  fill = G), alpha = 0.5) +
  stat_ellipse(aes(col = G), size = 3, alpha = 0.6)+
  geom_text(data=scbp, aes(x=(dbRDA1*2),y=(dbRDA2*2),label=Env),size=5)+
  scale_color_manual(values=c("goldenrod", "forestgreen", "royalblue"))+
  scale_fill_manual(values=c("goldenrod", "forestgreen", "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme 



plot(modelfv)

#bacteria

data.bac.t <- read.csv("16S_ASV_t.csv", row.names = 1)


data.bac1 <- subset(data.bac.t, rowSums(data.bac.t) >= 10)
#return to useful form
data.bac.raw <- t(data.bac1)



data.bac <- decostand(data.bac.raw, method="hellinger")

bac.tax.raw <- read.csv("16S_tax.csv", row.names = 1)


bac.all <- merge(t(data.bac), bac.tax.raw, by.x = 'row.names', 
                   by.y = 'row.names', all.x = TRUE, all.y = FALSE)

bac.tax <- data.frame(taxonomy = bac.all$taxonomy)


ddatb <- as.matrix(vegdist(data.bac, method = "jaccard"))


modelbg <- dbrda(ddatb ~ G + Condition(as.matrix(meta.fun[,c("PointX", "PointY")])), data = meta.fun, method="jaccard", add = TRUE)
modelbv <- dbrda(ddatb ~ . + Condition(as.matrix(meta.fun[,c("PointX", "PointY")])), data = vars, method="jaccard", add = TRUE)
modelb <- dbrda(ddatb ~1 + Condition(as.matrix(meta.fun[,c("PointX", "PointY")])), data = meta.fun, method="jaccard", add = TRUE)

anova(modelbg, by = "margin", type = 3)

RsquareAdj(modelbg)
anova(modelbv, by = "margin", type = 3)

summary(modelbg)
summary(modelbv)
summary(modelb)


ordistep(modelbv)

modelbv <- dbrda(formula = ddatb ~ MAOM_d13C + ecmq1 + Condition(as.matrix(meta.fun[, c("PointX", "PointY")])), data = vars,
                 add = TRUE, method = "jaccard")

plot(modelbv)

scb <- data.frame(scores(modelb))

bo <- ggplot(meta, aes(x = scb[,1], y = scb[,2])) +
  #geom_polygon(data = hulls, alpha = 0.1, size= 3, fill=NA) +
  geom_point(size =5, aes(color = G,  fill = G), alpha = 0.5) +
  #geom_segment(data=ecm.vec.sp.df,aes(x=0,xend=MDS1,y=0,yend=MDS2),
  #             arrow = arrow(length = unit(0.5, "cm")),colour="black", size = 1.3, alpha = 0.5) + 
  stat_ellipse(aes(col = G), size = 3, alpha = 0.6)+
  #geom_text(data=ecm.vec.sp.df, aes(x=MDS1,y=MDS2,label=Env),size=4, nudge_y = -0.02, nudge_x = -0.02)+
  #geom_point(shape=1, size = 3.25, colour = "black")+
  scale_color_manual(values=c("goldenrod", "forestgreen", "royalblue"))+
  scale_fill_manual(values=c("goldenrod", "forestgreen", "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme 



scbd.b <- data.frame(scores(modelbv, display = "sites"))
scbp.b <- data.frame(scores(modelbv, display = "bp"))[c("MAOM_d13C", "ecmq1"),]



c <- t(data.frame("MAOM δC***", "EMF diversity"))




scbp.b$Env <- c[,1]



bd <- ggplot(meta, aes(x = scbd.b[,1], y = scbd.b[,2])) +
  geom_point(size =5, aes(color = G,  fill = G), alpha = 0.5) +
  stat_ellipse(aes(col = G), size = 3, alpha = 0.6)+
  geom_text(data=scbp.b, aes(x=(dbRDA1*2),y=(dbRDA2*2),label=Env),size=5)+
  scale_color_manual(values=c("goldenrod", "forestgreen", "royalblue"))+
  scale_fill_manual(values=c("goldenrod", "forestgreen", "royalblue"))+
  ylab("")+ 
  xlab("")+ 
  mytheme 

quartz()


Fig5 <- ggarrange(fo, bo, fb, bd, ncol = 2)



#-------------- Table 1 -------------------


# Calculate distance matrix and inverse distance matrix
dist <- as.matrix(dist(cbind(meta$PointX, meta.fun$PointY)))
dist.inv <- 1/dist
diag(dist.inv) <- 0


morandf <- data.frame(matrix(NA, nrow = 32, ncol = 5))
names(morandf) <- c("Variable", "Moran_I", "unadjusted p-value", "F-value", "adjusted p-value")


var_names <- c("pH", "dw", "ppmC", "ppmN", "MBC", "MBN", "NH4", "NO3", "InorgN", 
               "pct_POM", "MAOM_d13C", "pct_MAOM_C", "MAOM_d15N", "pct_MAOM_N", "MAOM_CN", 
               "POM_d13C", "pct_POM_C", "POM_d15N", "pct_POM_N", "POM_CN", "bulk_d13C", 
               "bulk_pctC", "bulk_d15N", "bulk_pctN", "bulk_CN", "leaf_d15N", "leaf_pctN", 
               "leaf_d13C", "leaf_pctC", "leaf_CN", "CEC", "Pct_base_sat")


morandf$Variable <- var_names

#------- Moran's I p-value calculation and model fitting----------
for (i in 1:length(var_names)) {
  var <- var_names[i]
  
  # Moran's I p-value
  morandf[i, 2] <- Moran.I(meta.fun[[var]], dist.inv, na.rm = T)$p.value
  
  # LME model and ANOVA p-value
  model <- lme(fixed = as.formula(paste(var, "~ G")), data = meta.fun, random = ~ 1 | Type, method = "ML")
  morandf[i, 3] <- anova(model)[2, 4]  # unadjusted p-value
  morandf[i, 4] <- anova(model)[2, 3]  # F-value
  
  # Spatial correlation model (Gaussian correlation)
  spatial_model <- update(model, correlation = corGaus(form = ~ PointX + PointY), method = "ML")
  morandf[i, 5] <- anova(spatial_model)[2, 4]  # adjusted p-value
}



