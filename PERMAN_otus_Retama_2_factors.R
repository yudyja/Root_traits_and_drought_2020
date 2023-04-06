#PERMANOVA - Permutational Multivariate Analysis of Variance Using Distance Matrices

### https://www.youtube.com/watch?v=1QGI6u0BVnQ

#PERMANOVA 2 factores. Multivariables response (OTUS abundance) 
# Factor. Plant species (varias) y tto (drought,watered) 

# This script # Distance matrix
# NMDS
# PERMANOVA
# Homogeneity of multivariate dispersion
# SIMPER , Similarity percentages



require(vegan) #idem than library(vegan)
require("readxl")  


#tablas : otus abundance y variables categoricas (factores)

otus <- read_excel("otus_retama.xlsx", sheet = "otus_abundance") #nombre del sheet
categ <- read_excel("otus_retama.xlsx", sheet = "factores")


str(otus)
View(otus)

str(categ)
View(categ)


# No es necesario transformar.Pero probar.
# Double root transformed abundances####  Try to find a scale between 0-10 before use bray-curtis 

range(otus) #Scale between 0-453..Mejor transformar. 
range(otus^0.5) #with sqrt scale between 0-21
range(otus^0.25) # double sqrt ok# 

####### NMDS 

dist_otus<-vegdist(otus^0.25, method = 'bray')  #   DISTANCE MATRIX. I use a double root transformation To

nmds<-metaMDS(dist_otus)   # creando el NMDs# puedo agregar k=2 . ver stress. debajo de 0.2 is good
plot(nmds)  #plot(nmds, type ='text')  # de esta forma sale en lugar de puntos info de la variable
stressplot(nmds)


#############
#plot NMDS
op<-ordiplot(nmds, type ='n') #type n para tener un lienzo blanco# hay otras func como ggplot pero para luego colores y demas es complejo#

#points
cols = c('darkred','coral','green')               
points(nmds, cex = 3, pch = 16, col = cols[categ$Specie])

##decoration
ordispider (nmds, groups = categ$Specie, label = TRUE)
ordihull (nmds,groups = categ$Specie, lty = 'dotted')
legend ("bottomleft", pch = 16, col = cols, legend = levels(warra_env$position))


########
###### PERMANOVA ###########
###
## como 2 opciones:
##1. permanova tipico, 2. 

pmv <- adonis(otus^0.25 ~ specie * tto, method = "bray", perm = 999, data = categ)
pmv

##testa diferencias usando la CATEGORICA#
##tenemos interaccion specie*tto. p = 0.010
## plant species explains 49.6% of variance y tto (8%) . Ver R2


#fit model to  factors or environmental variables #ajusta los factores o env_var a la ordination. 

env.fit <- envfit(nmds,categ,na.rm = TRUE) #  x,y values como en PCA
                                           #Las "Otus" se mueven en direccion mis factores?. ver goodness of fit. 
env.fit




##############
### Homogeneity of multivariate dispersion#  Han de ser homogeneas... y 

##### POSTHOC ?######## 
########

### para specie## 

bd<-betadisper(dist_otus, categ$specie) #is a multivariate analogue of Levene's test for homogeneity of variances.
bd

permutest(bd, pairwise = TRUE) #If we cannot find a statistically different dispersion(p=mayor 0.05) so, assumption of homogeneity is met
#Pairwise=TRUE nos da info adittional, pero puedo no ponerlo

TukeyT<-TukeyHSD(bd) #Determine which groups are more variable# from Help=betadisper. se puede interpretar como POSTHOC? creo q si
TukeyT


### para tto 
bdtto<-betadisper(dist_otus, categ$tto) #is a multivariate analogue of Levene's test for homogeneity of variances.
bdtto

permutest(bdtto, pairwise = TRUE)

TukeyT<-TukeyHSD(bdtto) #Determine which groups are more variable# from Help=permutest.betadisper. POSTHOC?????
TukeyT

#Par la interaccion# 

# Hacer lo mismo que atras pero incluir en "categ" otra columna donde una sp+tto.

bdinter<-betadisper(dist_otus, categ$sp_tto) #is a multivariate analogue of Levene's test for homogeneity of variances.
bdtto

permutest(bdinter, pairwise = TRUE)

TukeyT<-TukeyHSD(bdinter) #Determine which groups are more variable# from Help=permutest.betadisper. POSTHOC?????
Tukey
#####






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


### FIN#### 




