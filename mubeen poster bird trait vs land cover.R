setwd("~/Desktop/iiser/9th Semester/Mubeen Poster NSAB22")
commmatrix<- read.csv("combinedallbirdslessthn4.csv")
birdtraits <-read.csv("Occupancy bird traits.xlsx - AVONET BirdLife.csv")
vegdata<- read.csv("Vegetation only 1ha for mubeen's poster.csv")

library(tidyverse)

"databird<-read.csv("~/Desktop/iiser/9th Semester/occupancy bird/collated bird data/Bird_Data_combined for TN report.csv")
str(databird)
sort(unique(databird$Site))
data1<-read.csv("~/Desktop/iiser/8 sem/occupancy/plot centroids _forest plots(1).csv")"

"colnames(data1)[1]<- "Site"
final_bird_total <- merge(databird,data1, by=c("Site"))
unique(final_bird_total$Site)
write.csv(final_bird_total, "~/Desktop/iiser/9th Semester/occupancy bird/collated bird data/finalirddata_compiled_for_ebird.csv")"

vegmatrix<- read.csv("2017_only_nilgiri_anamalai_palaniwithbufferdata.csv")
abcd<- vegmatrix %>%
  group_by(Site) %>%
  summarise(Acacia= sum(Acacia),Eucalyptus= sum(Eucalyptus),Others= sum(Others),
            Pine= sum(Pine),Shola= sum(Shola), Fresh.stumps = sum(Fresh.stumps), 
            Fresh.Cuts = sum(Fresh.Cuts),Logs = sum(Logs), Canopy = mean(Canopy), 
            Stems = mean(Stems), Moss = mean(Moss) , VS1 = sum(VS1), VS2 = sum(VS2), 
            VS3 = sum(VS3), VS4 = sum(VS4), VS5 = sum(VS5), VS6 = sum(VS6), VS7 = sum(VS7), 
            VS8 = sum(VS8), VS9 = sum(VS9), VS10 = sum(VS10))

names(commmatrix)[2]<- 'Site'
veg_bird_matrix<- merge(abcd, commmatrix, by=c( "Site"))

write.csv(veg_bird_matrix, "veg_bird_matrix.csv")


alltraitpca<- prcomp(birdtraits[c(7,9,10,11,12,13,16)], center = TRUE,scale. = TRUE)
summary(alltraitpca)
library(factoextra)

fviz_pca_ind(alltraitpca,
             repel = TRUE     # Avoid text overlapping
)

trophictraits<-prcomp(birdtraits[c(7,9,10)], center = TRUE,scale. = TRUE)
summary(trophictraits)
str(trophictraits)

fviz_pca_ind(trophictraits,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

var.trophic<- get_pca_ind(trophictraits)
var.trophic1<- as.data.frame(var.trophic$coord)
var.trophic1$Dim.1
bodysize<-NULL
bodysize$troph_pc1<-var.trophic1$Dim.1

ggplot() +
  geom_point(aes(x = var.trophic1$Dim.1, y = var.trophic1$Dim.2, fill = birdtraits$Trophic.Niche), colour='black', pch = 21,size = 2.2) + 
  coord_fixed(ratio = 3)+ labs(title = "Trophic traits across trophic niches of the birds", x = "PC1", y = "PC2") 






locomotorytraits<-prcomp(birdtraits[c(11,12,16)], center = TRUE,scale. = TRUE)
summary(locomotorytraits)

fviz_pca_ind(locomotorytraits,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

var.locomotory<- get_pca_ind(locomotorytraits)
var.locomotory1<- as.data.frame(var.locomotory$coord)
var.locomotory1$Dim.2

ggplot() +
  geom_point(aes(x = var.locomotory1$Dim.1, y = var.locomotory1$Dim.2, fill = birdtraits$Primary.Lifestyle), colour='black', pch = 21,size = 2.2) + 
  coord_fixed(ratio = 3) +labs(title = "Locomotory traits across lifestyles of the birds", x = "PC1", y = "PC2")


bodysize$locom_pc1<-var.locomotory1$Dim.1

bodysizetraits<-prcomp(as.data.frame(bodysize), center = TRUE,scale. = TRUE)
summary(bodysizetraits)





fviz_pca_ind(bodysizetraits,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

var.bodsize<- get_pca_ind(bodysizetraits)
var.bodsize1<- as.data.frame(var.bodsize$coord)
var.bodsize1$Dim.1

birdtrait_final<- NULL

birdtrait_final$Species<- birdtraits$Common.name
birdtrait_final$trophic<- var.trophic1$Dim.2
birdtrait_final$locomotory<- var.locomotory1$Dim.2
birdtrait_final$bodysize<- var.bodsize1$Dim.1
birdtrait_final<- as.data.frame(birdtrait_final)

ggplot() +
  geom_point(aes(x = var.bodsize1$Dim.1, y = var.bodsize1$Dim.2, fill = birdtraits$Primary.Lifestyle), colour='black', pch = 21,size = 2.2) + 
  coord_fixed(ratio = 1) +labs(title = "Body size across trophic niche of the birds", x = "PC1", y = "PC2")




#write.csv(birdtrait_final, "birdtrait_final.csv")

##practice
traits<-read.csv("birdtrait_final.csv", head=T, row.names=1)

library(picante)
library(ade4)

D<-vegdist(traits,"gower")

tree<-hclust(D,"average")
plot(tree)

tree.p<-as.phylo(tree)

comp<- read.csv("combinedallbirdslessthn4.csv", head=T, row.names=1)
FD<-pd(comp, tree.p, include.root=FALSE)




##With combined birds
traits<-read.csv("birdtrait_final.csv", head=T, row.names=1)

library(picante)
library(ade4)

D<-vegdist(traits,"gower")

tree<-hclust(D,"average")
plot(tree)

tree.p<-as.phylo(tree)
matrix<-read.csv("veg_bird_matrix.csv", head=T, row.names=1)


comp<- matrix[c(23:76)]
FD<-pd(comp, tree.p, include.root=FALSE)
FD$Type<-matrix$Type


FDplot<- ggplot(FD, aes(x=PD, y=Type)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(x="Functional Distance", y = "Habitat")+
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette="Dark2") 

SRplot<- ggplot(FD, aes(x=SR, y=Type)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(x="Functional Distance", y = "Habitat")+
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette="Dark2")


####invertivores####
inverttraits<-birdtraits[birdtraits$Trophic.Niche=="Invertivore",]


trophictraits<-prcomp(inverttraits[c(7,9,10)], center = TRUE,scale. = TRUE)
summary(trophictraits)

var.trophic<- get_pca_ind(trophictraits)
var.trophic1<- as.data.frame(var.trophic$coord)
var.trophic1$Dim.1
bodysize<-NULL
bodysize$troph_pc1<-var.trophic1$Dim.1


locomotorytraits<-prcomp(inverttraits[c(11,12,16)], center = TRUE,scale. = TRUE)
summary(locomotorytraits)


var.locomotory<- get_pca_ind(locomotorytraits)
var.locomotory1<- as.data.frame(var.locomotory$coord)
var.locomotory1$Dim.2


bodysize$locom_pc1<-var.locomotory1$Dim.1

bodysizetraits<-prcomp(as.data.frame(bodysize), center = TRUE,scale. = TRUE)
summary(bodysizetraits)


var.bodsize<- get_pca_ind(bodysizetraits)
var.bodsize1<- as.data.frame(var.bodsize$coord)
var.bodsize1$Dim.1


inverttrait_final<- NULL

inverttrait_final$Species<- inverttraits$Common.name
inverttrait_final$trophic<- var.trophic1$Dim.2
inverttrait_final$locomotory<- var.locomotory1$Dim.2
inverttrait_final$bodysize<- var.bodsize1$Dim.1
inverttrait_final<- as.data.frame(inverttrait_final)

#write.csv(inverttrait_final, "inverttrait_final.csv")

traits<-read.csv("inverttrait_final.csv", head=T, row.names=1)

D<-vegdist(traits,"gower")

tree<-hclust(D,"average")
plot(tree)

tree.p<-as.phylo(tree)

matrix<-read.csv("veg_bird_matrix.csv", head=T, row.names=1)
matrix$VP0_2<-matrix$VS1+matrix$VS2+matrix$VS3+matrix$VS4
matrix$VP3_5<-matrix$VS7+matrix$VS8+matrix$VS9+matrix$VS10





comp<- matrix[c(23:76)]
FD<-pd(comp, tree.p, include.root=FALSE)
FD$Type<-matrix$Type
FD$Canopy<-matrix$Canopy
FD$VP3_5<-matrix$VP3_5
FD$VP0_2<-matrix$VP0_2



FDplot_invert<- ggplot(FD, aes(x=PD, y=Type)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(x="Functional Distance", y = "Habitat")+
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette="Dark2") 
  
mod<-aov(PD~Type, data=FD)
summary(mod)


####Null model
Nullmat<- randomizeMatrix(comp, null.model = c("richness"), iterations = 1000)   

FD_null<-pd(Nullmat, tree.p, include.root=FALSE)
FD_null$Type<-matrix$Type
FD_null$Canopy<-matrix$Canopy
FD_null$VP3_5<-matrix$VP3_5
FD_null$VP0_2<-matrix$VP0_2

result_shola = wilcox.test(FD$PD[FD$Type=="Shola"], FD_null$PD[FD$Type=="Shola"], paired = TRUE)
result_eucalyptus = wilcox.test(FD$PD[FD$Type=="Eucalyptus"], FD_null$PD[FD$Type=="Eucalyptus"], paired = TRUE)
result_Acacia = wilcox.test(FD$PD[FD$Type=="Acacia"], FD_null$PD[FD$Type=="Acacia"], paired = TRUE)
result_Pine = wilcox.test(FD$PD[FD$Type=="Pine"], FD_null$PD[FD$Type=="Pine"], paired = TRUE)
result_Others = wilcox.test(FD$PD[FD$Type=="Others"], FD_null$PD[FD$Type=="Others"], paired = TRUE)
result_Mixed_Shola = wilcox.test(FD$PD[FD$Type=="Mixed_Shola"], FD_null$PD[FD$Type=="Mixed_Shola"], paired = TRUE)
result_Rest = wilcox.test(FD$PD[FD$Type=="Rest"], FD_null$PD[FD$Type=="Rest"], paired = TRUE)
result_None = wilcox.test(FD$PD[FD$Type=="None"], FD_null$PD[FD$Type=="None"], paired = TRUE)


FD_invert<- ggplot(FD, aes(x=Type, y=PD)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(title = "Observed functional diversity (FD) for invertivores", x="Functional Distance", y = "Habitat")+
  theme(text = element_text(size = 15)) +
  scale_fill_brewer(palette="Dark2") 
  

(FD_invert  
+ annotate("text", x = 9, y = 1.5, label = "***", size = 8) 
+ annotate("text", x = 3, y = 1.5, label = "**", size = 8) 
+ annotate("text", x = 4, y = 1.5, label = "**", size = 8) 
+ annotate("text", x = 6, y = 1.5, label = "***", size = 8) 

)

mod<-lm(PD~Type, data=FD)
summary(mod)
plot(FD$PD~FD$Canopy)
plot(FD$PD~FD$VP0_2)
plot(FD$PD~FD$VP3_5)


SRplot<- ggplot(FD, aes(x=SR, y=Type)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(x="Functional Distance", y = "Habitat")+
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette="Dark2")


####frugivores+omnivores####
frugtraits<-birdtraits[birdtraits$Trophic.Niche=="Frugivore" | birdtraits$Trophic.Niche=="Omnivore",]


trophictraits<-prcomp(frugtraits[c(7,9,10)], center = TRUE,scale. = TRUE)
summary(trophictraits)

var.trophic<- get_pca_ind(trophictraits)
var.trophic1<- as.data.frame(var.trophic$coord)
var.trophic1$Dim.1
bodysize<-NULL
bodysize$troph_pc1<-var.trophic1$Dim.1


locomotorytraits<-prcomp(frugtraits[c(11,12,16)], center = TRUE,scale. = TRUE)
summary(locomotorytraits)


var.locomotory<- get_pca_ind(locomotorytraits)
var.locomotory1<- as.data.frame(var.locomotory$coord)
var.locomotory1$Dim.2

ggplot() +
  geom_point(aes(x = var.locomotory1$Dim.1, y = var.locomotory1$Dim.2, fill = frugtraits$Primary.Lifestyle), colour='black', pch = 21,size = 2.2) + 
  coord_fixed(ratio = 3) 


bodysize$locom_pc1<-var.locomotory1$Dim.1

bodysizetraits<-prcomp(as.data.frame(bodysize), center = TRUE,scale. = TRUE)
summary(bodysizetraits)


var.bodsize<- get_pca_ind(bodysizetraits)
var.bodsize1<- as.data.frame(var.bodsize$coord)
var.bodsize1$Dim.1

frugtrait_final<- NULL

frugtrait_final$Species<- frugtraits$Common.name
frugtrait_final$trophic<- var.trophic1$Dim.2
frugtrait_final$locomotory<- var.locomotory1$Dim.2
frugtrait_final$bodysize<- var.bodsize1$Dim.1
frugtrait_final<- as.data.frame(frugtrait_final)

#write.csv(frugtrait_final, "frugtrait_final.csv")

traits<-read.csv("~/Desktop/iiser/9th Semester/Mubeen Poster NSAB22/frugtrait_final.csv", head=T, row.names=1)

D<-vegdist(traits,"gower")

tree<-hclust(D,"average")
plot(tree)

tree.p<-as.phylo(tree)

###functional dispersion
library(FD)
frugcomp<- as.matrix(comp[,c(9,17,23,25,27,28,34,37,41,44,45,46,48,53,54)])
frugcomp<-frugcomp[rowSums(frugcomp[])>0,]
colnames(frugcomp)<- rownames(traits)
fun.dist <- gowdis(traits)


ex1 <- fdisp(fun.dist, frugcomp)
ex1$FDis

matrix<-read.csv("veg_bird_matrix.csv", head=T, row.names=1)
matrix$VP0_2<-matrix$VS1+matrix$VS2+matrix$VS3+matrix$VS4
matrix$VP3_5<-matrix$VS7+matrix$VS8+matrix$VS9+matrix$VS10




comp<- matrix[c(23:76)]
FD<-pd(comp, tree.p, include.root=FALSE)
FD$Type<-matrix$Type
FD$Canopy<-matrix$Canopy
FD$VP3_5<-matrix$VP3_5
FD$VP0_2<-matrix$VP0_2
FD1<-FD[-rownames(FD)[c("AUW641", "AXO700", "BDZ967", "BGH965", "BHB840")],]
row.names.FD<-c("AUW641", "AXO700", "BDZ967", "BGH965", "BHB840")
FD<- FD[!(row.names(FD) %in% row.names.FD),]
length(FD1$PD)

FD$FDis<-ex1$FDis

mod<-lm(FDis~Type, data=FD)
summary(mod)
library(lme4)
library(nlme)

mod_lmer<-lm(PD~Type,data=FD)
summary(mod_lmer)

model <- aov(PD~Type,data=FD)
summary(model)
TukeyHSD(model, conf.level=.95)


plot(FD$FDis~FD$Canopy)
plot(FD$PD~FD$VP0_2)
plot(FD$PD~FD$VP3_5)


####Null model
Nullmat<- randomizeMatrix(comp, null.model = c("richness"), iterations = 1000)   

FD_null<-pd(Nullmat, tree.p, include.root=FALSE)
FD_null$Type<-matrix$Type
FD_null$Canopy<-matrix$Canopy
FD_null$VP3_5<-matrix$VP3_5
FD_null$VP0_2<-matrix$VP0_2

result_shola = wilcox.test(FD$PD[FD$Type=="Shola"], FD_null$PD[FD$Type=="Shola"], paired = TRUE)
result_eucalyptus = wilcox.test(FD$PD[FD$Type=="Eucalyptus"], FD_null$PD[FD$Type=="Eucalyptus"], paired = TRUE)
result_Acacia = wilcox.test(FD$PD[FD$Type=="Acacia"], FD_null$PD[FD$Type=="Acacia"], paired = TRUE)
result_Pine = wilcox.test(FD$PD[FD$Type=="Pine"], FD_null$PD[FD$Type=="Pine"], paired = TRUE)
result_Others = wilcox.test(FD$PD[FD$Type=="Others"], FD_null$PD[FD$Type=="Others"], paired = TRUE)
result_Mixed_Shola = wilcox.test(FD$PD[FD$Type=="Mixed_Shola"], FD_null$PD[FD$Type=="Mixed_Shola"], paired = TRUE)
result_Rest = wilcox.test(FD$PD[FD$Type=="Rest"], FD_null$PD[FD$Type=="Rest"], paired = TRUE)
result_None = wilcox.test(FD$PD[FD$Type=="None"], FD_null$PD[FD$Type=="None"], paired = TRUE)







FDplotfrug<- ggplot(FD, aes(y=PD, x=Type)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(title = "Observed functional diversity (FD) for invertivores", x="Functional Distance", y = "Habitat")+
  theme(text = element_text(size = 15)) +
  scale_fill_brewer(palette="Dark2") 

(FDplotfrug  
  + annotate("text", x = 9, y = 1.5, label = "***", size = 8) 
  + annotate("text", x = 1, y = 1.5, label = "**", size = 8) 
  + annotate("text", x = 3, y = 1.5, label = "***", size = 8) 
  + annotate("text", x = 5, y = 1.5, label = "**", size = 8) 
  
)




SRplotfrug<- ggplot(FD, aes(x=SR, y=Type)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(x="Functional Distance", y = "Habitat")+
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette="Dark2")

####frugivores####
frugtraits<-birdtraits[birdtraits$Trophic.Niche=="Frugivore",]


trophictraits<-prcomp(frugtraits[c(7,9,10)], center = TRUE,scale. = TRUE)
summary(trophictraits)

var.trophic<- get_pca_ind(trophictraits)
var.trophic1<- as.data.frame(var.trophic$coord)
var.trophic1$Dim.1
bodysize<-NULL
bodysize$troph_pc1<-var.trophic1$Dim.1


locomotorytraits<-prcomp(frugtraits[c(11,12,16)], center = TRUE,scale. = TRUE)
summary(locomotorytraits)


var.locomotory<- get_pca_ind(locomotorytraits)
var.locomotory1<- as.data.frame(var.locomotory$coord)
var.locomotory1$Dim.2

ggplot() +
  geom_point(aes(x = var.locomotory1$Dim.1, y = var.locomotory1$Dim.2, fill = frugtraits$Primary.Lifestyle), colour='black', pch = 21,size = 2.2) + 
  coord_fixed(ratio = 3) 


bodysize$locom_pc1<-var.locomotory1$Dim.1

bodysizetraits<-prcomp(as.data.frame(bodysize), center = TRUE,scale. = TRUE)
summary(bodysizetraits)


var.bodsize<- get_pca_ind(bodysizetraits)
var.bodsize1<- as.data.frame(var.bodsize$coord)
var.bodsize1$Dim.1

frugtrait_final<- NULL

frugtrait_final$Species<- frugtraits$Common.name
frugtrait_final$trophic<- var.trophic1$Dim.2
frugtrait_final$locomotory<- var.locomotory1$Dim.2
frugtrait_final$bodysize<- var.bodsize1$Dim.1
frugtrait_final<- as.data.frame(frugtrait_final)

write.csv(frugtrait_final, "frugtrait_final1.csv")

traits<-read.csv("frugtrait_final1.csv", head=T, row.names=1)

D<-vegdist(traits,"gower")

tree<-hclust(D,"average")
plot(tree)

tree.p<-as.phylo(tree)

###functional dispersion
library(FD)
frugcomp<- as.matrix(comp[,c(9,17,23,25,27,28,34,37,41,44,45,46,48,53,54)])
frugcomp<-frugcomp[rowSums(frugcomp[])>0,]
colnames(frugcomp)<- rownames(traits)
fun.dist <- gowdis(traits)


ex1 <- fdisp(fun.dist, frugcomp)
ex1$FDis

matrix<-read.csv("veg_bird_matrix.csv", head=T, row.names=1)
matrix$VP0_2<-matrix$VS1+matrix$VS2+matrix$VS3+matrix$VS4
matrix$VP3_5<-matrix$VS7+matrix$VS8+matrix$VS9+matrix$VS10




comp<- matrix[c(23:76)]
FD<-pd(comp, tree.p, include.root=FALSE)
FD$Type<-matrix$Type
FD$Canopy<-matrix$Canopy
FD$VP3_5<-matrix$VP3_5
FD$VP0_2<-matrix$VP0_2
FD1<-FD[-rownames(FD)[c("AUW641", "AXO700", "BDZ967", "BGH965", "BHB840")],]
row.names.FD<-c("AUW641", "AXO700", "BDZ967", "BGH965", "BHB840")
FD<- FD[!(row.names(FD) %in% row.names.FD),]
length(FD1$PD)

FD$FDis<-ex1$FDis

mod<-lm(PD~Type, data=FD)
summary(mod)
library(lme4)
library(nlme)

mod_lmer<-lmer(PD~Type+(1|VP0_2)+(1|Canopy),data=FD)
summary(mod_lmer)


plot(FD$FDis~FD$Canopy)
plot(FD$PD~FD$VP0_2)
plot(FD$PD~FD$VP3_5)

library(ggsignif)
FDplotfrug1<- ggplot(FD, aes(x=PD, y=Type)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(x="Functional Dispersion", y = "Habitat")+
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette="Dark2") 

SRplotfrug<- ggplot(FD, aes(x=SR, y=Type)) + 
  geom_boxplot(fill="gray")+
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black")+
  labs(x="Functional Distance", y = "Habitat")+
  theme(text = element_text(size = 20)) +
  scale_fill_brewer(palette="Dark2")


