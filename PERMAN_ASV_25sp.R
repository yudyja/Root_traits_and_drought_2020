#PERMANOVA - Permutational Multivariate Analysis of Variance Using Distance Matrices

### https://www.youtube.com/watch?v=1QGI6u0BVnQ

#PERMANOVA 2 factores. Multivariables response (Sequence variants) 
# Factor. Plant species (varias) y tto (drought,watered) 

# This script # Distance matrix
# NMDS
# PERMANOVA
# Homogeneity of multivariate dispersion
# SIMPER , Similarity percentages



require(vegan) #idem than library(vegan)
require("readxl")  


#### without muestras from T0##
#####
######

#tablas : otus abundance y variables categoricas (factores)

# 25 species juntas# 

asv25 <- read_excel("seqtab.nochim.fung.xlsx", sheet ="Tf_ASV") #nombre del sheet
categ25 <- read_excel("seqtab.nochim.fung.xlsx", sheet = "factors")

str(asv25)
str(categ25)


# Double root transformed abundances####  Try to find a scale between 0-10 before use bray-curtis 
range(asv25) #Scale between 0-453..Mejor transformar. 
range(asv25^0.5) #with sqrt scale between 0-21
range(asv25^0.25) # double sqrt ok# 


####### NMDS ... variable transformada aqui en R  ^0.25 (da mejor transformada que sin)

dist_bray<-vegdist(asv25^0.25, method = 'bray')  #   DISTANCE MATRIX. I use a double root transformation To
nmds25<-metaMDS(dist_bray)   # creando el NMDs# puedo agregar k=2 . ver stress. debajo de 0.2 is good
plot(nmds25)  #plot(nmds, type ='text')  # de esta forma sale en lugar de puntos info de la variable
stressplot(nmds25)

pmv25 <- adonis(asv25^0.25 ~ species * tto, method = "bray", perm = 999, data = categ25) #PERMANOVA
pmv25   # sp*tto y sp tiene effecto. 


env.fit <- envfit(nmdsgrass,categ,na.rm = TRUE) #  x,y values como en PCA
env.fit

##
bdinter<-betadisper(dist_bray, categ25$species_tto) #is a multivariate analogue of Levene's test for homogeneity of variances.
bdinter

permutest(bdinter, pairwise = TRUE)

Tukeyinter<-TukeyHSD(bdinter) #Determine which groups are more variable# from Help=permutest.betadisper. POSTHOC?????
Tukeyinter
#####
##########


### Plant functional group as FACTOR # 

str(categ25)
pmv_plantgroup <- adonis(asv25^0.25 ~ plantgroup * tto, method = "bray", perm = 999, data = categ25) #PERMANOVA
pmv_plantgroup   





#### SOLO GRASSES
#######
asvgrasses <- read_excel("seqtab.nochim.fung.xlsx", sheet ="grasses_ASV") #nombre del sheet
categgrasses <- read_excel("seqtab.nochim.fung.xlsx", sheet = "factors_grasses")

# Double root transformed abundances####  Try to find a scale between 0-10 before use bray-curtis 
range(asvgrasses) #Scale between 0-453..Mejor transformar. 
range(asvgrasses^0.5) #with sqrt scale between 0-21
range(asvgrasses^0.25) # double sqrt ok# 


####### NMDS ... variable transformada aqui en R  ^0.25 (da mejor transformada que sin)

dist_bray_grasses<-vegdist(asvgrasses^0.25, method = 'bray')  #   DISTANCE MATRIX. I use a double root transformation To
nmdsgrasses<-metaMDS(dist_bray_grasses)   # creando el NMDs# puedo agregar k=2 . ver stress. debajo de 0.2 is good
plot(nmds25)  #plot(nmds, type ='text')  # de esta forma sale en lugar de puntos info de la variable


pmvgrasses <- adonis(asvgrasses^0.25 ~ species * tto, method = "bray", perm = 999, data = categgrasses) #PERMANOVA
pmvgrasses   # sp*tto y sp tiene effecto. 


envfit_grasses <- envfit(nmdsgrasses,categgrasses,na.rm = TRUE) #  x,y values como en PCA
#Las "Otus" se mueven en direccion mis factores?. ver goodness of fit. 
envfit_grasses



##por specie
bdgrasses_sp<-betadisper(dist_bray_grasses, categgrasses$species_tto) #is a multivariate analogue of Levene's test for homogeneity of variances.

permutest(bdinter, pairwise = TRUE)

Tukey_grasses_sp<-TukeyHSD(bdgrasses_sp) #Determine which groups are more variable# from Help=permutest.betadisper. POSTHOC?????
Tukey_grasses_sp

##por interaction
bdspgrasses<-betadisper(dist_bray_grasses, categgrasses$species) #is a multivariate analogue of Levene's test for homogeneity of variances.
bdspgrasses

Tukeygrasses<-TukeyHSD(bdspgrasses) #Determine which groups are more variable# from Help=permutest.betadisper. POSTHOC?????
Tukeygrasses





##########FIN 

# adicional cosas q quizas puedan usarse


#fit model to  factors or environmental variables #ajusta los factores o env_var a la ordination. 

env.fit <- envfit(nmdsgrass,categ,na.rm = TRUE) #  x,y values como en PCA
#Las "Otus" se mueven en direccion mis factores?. ver goodness of fit. 
env.fit



##### SIMPER ###### 

#### Here, which species are responsible for the differences between groups that you observed. 
####Simper decomposes the bray curtis distance.  

?simper
sim<-simper(otus^0.25, group = categ$tto)
sim
summary(sim)

# average = contribution to dissimilarity between levels of my factor
# sd      = standard deviation of contributio (is the species response consistent?)
# ratio   = ratio between average and sd (higher ratio means high, consistent contribution)
# av.     = average abundance per groups 
#cumsum   = cumulative contribution (rule of thumb:species till 70% are investigated)


# Many species contribute to the difference between levels del factor
# sp 19 decreases (6.29 to 2.49) and has the highest contribution. (see sd and ratio) is the first that appear
# sp 37 increases (3.24 to 7.09)
# sp12 decreases to 0 (2.78 to 0)  and is also consistently (high contribution to sd ratio)


###  FIn 


#######            
