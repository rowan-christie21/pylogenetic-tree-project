#5/2/2020
#Calculating CG and amino acid content of tetM sequences

#-------------------------------------------
# load libraries

library(ggplot2)
library(dplyr)
library(wesanderson)
library(RColorBrewer)
library(MuMIn)
library(car)
library(nlme)  
library(ncf)
library(lme4)
library(Matrix)
library(coefplot)
library(MASS)
library(genbankr)


#-------------------------------------------

tetm_strains <- readxl::read_xlsx("D:/Bioinformatics-Projects/phylogentic-tree-project-02292020/bacterial tetM datasets calculated 5_3_2020.xlsx", sheet=1)
tetm_strains <- tetm_strains[1:37,]

cor(tetm_strains$`Percent CG`, tetm_strains$`Percent identity`)
#-0.4651386

cor(tetm_strains$`Percent CG`, tetm_strains$`tetracycline µg/ml`)
#-0.859289

t.test(tetm_strains$`Percent CG`, tetm_strains$`tetracycline µg/ml`)
#t = -4.3667, df = 36, p-value = 0.0001021

plot(tetm_strains$`Percent CG`, tetm_strains$`tetracycline µg/ml`)

cg_mic <- ggplot(tetm_strains, aes(`Percent CG`, `tetracycline µg/ml`)) +
  geom_point() +
  geom_smooth(method=lm, size=2) +
  scale_x_continuous(name = "Percent CG content", expand = c(0,0)) +
  scale_y_continuous(name = "Tetracycline µg/ml", expand = c(0,0)) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(size = 27, margin=margin(0,0,15,0)),
        axis.text.x=element_text(colour="black", size = 18),
        axis.text.y=element_text(colour="black", size = 18),
        axis.title.x = element_text(size = 23, margin=margin(15,0,0,0)),
        axis.title.y = element_text(size = 23, margin=margin(0,15,0,0)))

cg_mic

png(filename = paste("D:/Bioinformatics-Projects/phylogentic-tree-project-02292020/cg_mic ",Sys.Date(),".png", sep = ''), width = 600, height = 454)
cg_mic
dev.off()

cor(tetm_strains$`Percent identity`, tetm_strains$`tetracycline µg/ml`)
#0.2859238

cor(tetm_strains$year, tetm_strains$`tetracycline µg/ml`, use = "complete.obs")
#-0.4358124

species <- ggplot(tetm_strains, aes(Species, `Percent CG`, size=`tetracycline µg/ml`, color=Species)) +
  geom_point() +
  #geom_smooth(method=lm, size=2) +
  #scale_x_continuous(name = "Percent CG content", expand = c(0,0)) +
  scale_y_continuous(name = "Percent CG Content", expand = c(0,0)) +
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title=element_text(size = 27, margin=margin(0,0,15,0)),
        axis.text.x=element_text(colour="black", size = 15, angle = 45, hjust = 1),
        axis.text.y=element_text(colour="black", size = 18),
        axis.title.x = element_text(size = 23, margin=margin(15,0,0,0)),
        axis.title.y = element_text(size = 23, margin=margin(0,15,0,0)))

species

png(filename = paste("D:/Bioinformatics-Projects/phylogentic-tree-project-02292020/species ",Sys.Date(),".png", sep = ''), width = 900, height = 600)
species
dev.off()


