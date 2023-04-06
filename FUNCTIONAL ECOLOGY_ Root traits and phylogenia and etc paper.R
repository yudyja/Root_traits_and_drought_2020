install.packages("regclass")
install.packages("broom")
install.packages("lmodel2")

library(tidyverse)
library(ggnewscale)
library(vegan)
library(ape)
library(picante)
library(phytools)
library(lavaan)
library(ggplot2)
library(GGally)
library(lmodel2)
library(ggpubr)
library(ggrepel)


#### PHYLOGENETIC ANALISIS 


treeEurope<-read.tree("DaPhnE_01.tre")  

#4 Crear mi arbol. 
Wp4plants<-c("Achillea_millefolium","Artemisia_campestris","Hieracium_pilosella",
             "Daucus_carota","Berteroa_incana","Scabiosa_canescens","Silene_vulgaris",
             "Hypericum_perforatum","Plantago_lanceolata","Armeria_maritima",
             "Galium_verum","Rumex_thyrsiflorus","Trifolium_repens","Festuca_brevipila",
             "Ranunculus_acris","Vicia_cracca","Festuca_rubra","Potentilla_argentea",
             "Anthoxanthum_odoratum","Holcus_lanatus","Arrhenatherum_elatius",
             "Lolium_perenne","Medicago_lupulina","Dactylis_glomerata","Poa_angustifolia")
Wp4plants<-Wp4plants[Wp4plants!="Scabiosa_canescens"]#This plant did not grow

#5 remover(drop.tip) del treeEurope todo aquello que NO ("!") este EN ("in")Wp4plants. Quitar todo lo que no sean mis plantas

WP4_PhyloTree<-drop.tip(treeEurope,
                        treeEurope$tip.label[!treeEurope$tip.label%in%Wp4plants])

#6 Ver mi tree y comprpbar arbol ok.
plot(WP4_PhyloTree) 
#Loading tree


#Loading the traits

traits <- read.csv("traits_25grasslands.csv",stringsAsFactors = F)
#crear nueva columna species nombre completo
traits$spcompleto <- NA # nueva columna empty
traits$spcompleto[which(traits$Sp=="Anthoxanthum")]<-"Anthoxanthum odoratum"
traits$spcompleto[which(traits$Sp=="Berteroa")]<-"Berteroa incana"
traits$spcompleto[which(traits$Sp=="Daucus")]<-"Daucus carota"
traits$spcompleto[which(traits$Sp=="F.brevipila")]<-"Festuca brevipila"
traits$spcompleto[which(traits$Sp=="Achillea")]<-"Achillea millefolium"
traits$spcompleto[which(traits$Sp=="Potentilla")]<-"Potentilla argentea"
traits$spcompleto[which(traits$Sp=="Plantago")]<-"Plantago lanceolata"
traits$spcompleto[which(traits$Sp=="Hieracium")]<-"Hieracium pilosella"
traits$spcompleto[which(traits$Sp=="Artemisia")]<-"Artemisia campestris"
traits$spcompleto[which(traits$Sp=="Holcus")]<-"Holcus lanatus"
traits$spcompleto[which(traits$Sp=="Arrhenatherum")]<-"Arrhenatherum elatius"
traits$spcompleto[which(traits$Sp=="Vicia")]<-"Vicia cracca"
traits$spcompleto[which(traits$Sp=="Hypericum")]<-"Hypericum perforatum"
traits$spcompleto[which(traits$Sp=="F.rubra")]<-"Festuca rubra"
traits$spcompleto[which(traits$Sp=="Galium")]<-"Galium verum"
traits$spcompleto[which(traits$Sp=="Poa")]<-"Poa angustifolia"
traits$spcompleto[which(traits$Sp=="Silene")]<-"Silene vulgaris"
traits$spcompleto[which(traits$Sp=="Trifolium")]<-"Trifolium repens"
traits$spcompleto[which(traits$Sp=="Dactylis")]<-"Dactylis glomerata"
traits$spcompleto[which(traits$Sp=="Rumex")]<-"Rumex thyrsiflorus"
traits$spcompleto[which(traits$Sp=="Medicago")]<-"Medicago lupulina"
traits$spcompleto[which(traits$Sp=="Lolium")]<-"Lolium perenne"
traits$spcompleto[which(traits$Sp=="Ranunculus")]<-"Ranunculus acris"
traits$spcompleto[which(traits$Sp=="Armeria")]<-"Armeria maritima"
traits$spcompleto <- gsub(" ","_",traits$spcompleto)
#Organizing the table in alphabetical order
traits<-traits[order(traits$spcompleto),]



#I. Obteniendo los Eigenvectors siguiendo el ejemplo de Stefan

# los eigenvectors (direction) y tienen eigenvalues (magnitud). Principal componentes son los eigenvectors con los 
#largest eigenvalues (selecciono casi siempre los dos primeros) Ver el video de https://www.youtube.com/watch?v=g-Hb26agBFg

#1. Obteniendo la matriz de disimilutd
dist.tree <- cophenetic(WP4_PhyloTree) 

#2. Haciendo el PCoA de esta disimilitud
pcoa.dist <- cmdscale(dist.tree, k= dim(dist.tree)[1]-1, eig=T) ##   son 24 species... y solo son posible 23 comparaciones, 23 eigenvalues
#k is the maximum dimension of the space which the data are to be represented in; must be in {1, 2, …, n-1}.
#

#3. Identificar los ejes eigenvalues que se usarabn (manualmente)
pcoa.dist$eig # Si sumo todos los eigen y eso es el 100% pues se cuales ejes suman el 80% (regla de 3) 
              # abajo se saca directamente aquellos q son el 80%) Eje 1 (339744.9) que corresponde a 47.29% 

pcoa.dist[[2]]/sum(pcoa.dist[[2]])*100 # hace lo mismo que arriba y me da el porcentaje explicado de cada eje directamente 47.15% maso lo mismo

#3.1. Identificar  cuales seran los ejes eigenvalues que se usaran. Stefan uso (arbitrariarmente) com criterio 
# eigenvalues que en suma explican 80% de la variación de los datos.
phyl.var1 <- 
  1:length(which(cumsum(pcoa.dist[[2]]/sum(pcoa.dist[[2]]))<0.8)) #[[]] seleccionar 
#de la pcoa.dist la "segunda" parte de los resultados son los eigenvalues y esos son los que quiero seleccionar para hacer el calculo

# cumulative sum de los eigenvalues / total. Y de estos cuales eigen values explican hasta el 80%)
# Han salido los 4 primeros 


#3.2. Seleccionar los eigenvectors (loadings) que corresponden con los eigenvalues identfiicados
phylogeny_vectors<-pcoa.dist[[1]][,phyl.var1]  # A esa "primera parte" de pcoa.dist que son los eigenvectores ($points) seleccionar
#tambien los 4 primeros (phyl.var1 means los 4 primeros)
phylogeny_vectors<-as.data.frame(phylogeny_vectors)
names(phylogeny_vectors)<-c("PCoA1","PCoA2","PCoA3","PCoA4") # ponerle nombres.. PCoA
phylogeny_vectors$spcompleto<-rownames(phylogeny_vectors) #la primera columna del phylogeny vectors es el nombre de las especies,
#entonces crear una nueva columna llamada "spcompleto" que tenga ese rownames
phylogeny_info<-phylogeny_vectors # copia de phylogeny vectors 


#Pasting the PCoA axes to the trait table
names(traits)
names(phylogeny_info)

traits_phylogeny<-left_join(traits,phylogeny_info,by="spcompleto")   # Unir las tablas por"spcompleto" porq la comparten  
traits_phylogeny$Units<-NULL #### descripcion de la tabla columna 1 (donde tengo mis observaciones)

#write.csv(traits_phylogeny,"to_extract_phylo_axes.csv")

#########################################
#####   MULTIVARIATE ANALYSES

#Ecologically meaningful transformations for ordination of species data. Legendre and Gallagher 2001 Oecologia

#Notes. Researchgate.If it is an ordination technique using euclidean distance, you need to transform your 
#data into euclidean space. This is what Hellinger does, and enables you to use your community data 
#in multivariate analyses such as PCA, RDA and CCA. 
#You can use Bray-Curtis distance only in distance-based RDA.
#Numerical Ecology with R (use R series) by Daniel Borcard is very useful in this context

#### RDA analyses (consideraciones)

# 1. Seleccionar las variables
#2. Ver la MULTICOLINEALIDAD sort(vif(variab)). Borcard et al 2011, Numerical ecology with R, P175
# argue that in the case of RDA, VIFs > 10 should be avoided. (yo he quitado la variable q daba ruido)  
#3. Check GRADIENT LENGTH. See that the axis lengths are below 2 SD . decorana (variab)
# if the axis is over 2 transformate HELLINGER
#dataVAR_he <- decostand(dataVAR, "hellinger") # SOLO numerical var
#decorana(dataVAR_he)
#sort(vif(dataVAR_he)) # NO funciono. ver el q use abajo sacado del libro Bordcard et al., 2018
#transformation decreases the gradient length. Use the original data if the gradient length is alreday gut. 
#4. rda function. #From SUMMARY. constrained variance = explained for the model, unconstrained = no explained (0.26=26%)
# Cuanta  varianza es explicada por cada RDA axis. Accumulated CONSTRAINED eigenvalues
# proportion explained.
#Ojo con "select". a veces se confunde con el "select" de MASS. por eso especifique es el dplyr
#dplyr::select(dataR, var1,var2)
#see details of plot.cca vegan about "species"(mis columns), "sites"(mis rows), cn(centroids)etc 
#####
str(traits)
######### ROOT TRAITS + PHYLOGENIA 

root_var<- traits_phylogeny[,c("rootN","rootC","sqrt_root","sqrt_diameter","logRTD","SRSA")] # quite SRL porq era multicolineal con SRSA
# transformando estas mejora un poco los resultados de lo q explica el RDA 

vif <- diag(solve(cor(root_var))) # multicolieanility
decorana(root_var) # gradient lenght

rda.root <- rda(root_var ~ Tto*spcompleto+Condition(PCoA1,PCoA2,PCoA3,PCoA4), #
                data=traits_phylogeny, scale = TRUE)
model1_5<-rda.root # porq de ahi opara abajo todo se hizo llamandolo model_1.5

#model1_5<-  ESTE ES LO MISMO>>> SOLO Q PRIMERO QUERIA CONFIRMAR MULTICOLI Y GRADIENT LENGTH 
 # rda(traits_phylogeny[,c("rootN","rootC","sqrt_root","sqrt_diameter","logRTD","SRSA")] ~ Tto*spcompleto+Condition(PCoA1,PCoA2,PCoA3,PCoA4), #
#    data=traits_phylogeny, scale = TRUE)

model1_5#:
#                Inertia Proportion Rank
# Total          6.0000     1.0000     
# Conditional    1.4587     0.2431    1
# Constrained    2.8598     0.4766    6
# Unconstrained  1.6815     0.2803    6
# Inertia is correlations 
# Some constraints were aliased because they were collinear (redundant)
# 
# Eigenvalues for constrained axes:
#   RDA1   RDA2   RDA3   RDA4   RDA5   RDA6 
# 1.0227 0.9257 0.3816 0.2345 0.2096 0.0856 
# 
# Eigenvalues for unconstrained axes:
#   PC1    PC2    PC3    PC4    PC5    PC6 
# 0.7033 0.3690 0.2527 0.2316 0.1048 0.0202 
summary(model1_5)[["cont"]][["importance"]][,c(1:2)]
#                       RDA1      RDA2
# Eigenvalue            1.0226929 0.9257241
# Proportion Explained  0.2251999 0.2038471
# Cumulative Proportion 0.2251999 0.4290470

coef(model1_5) #canonical coeficients from rda result 
RsquareAdj(model1_5)$r.squared #unadjusted and adjusted R^2 retreived from the rda result
RsquareAdj(model1_5)$adj.r.squared 

anova.cca(model1_5)
anova.cca(model1_5,by="term")  # Tb usar by axis 
#This model is significant overall and also treatment, species and interaction; 
#                  Df Variance F        Pr(>F)    
# Tto              1  0.02681  2.9979   0.008 ** 
# spcompleto      22  2.50996  12.7557  0.001 ***
# Tto:spcompleto  23  0.32299  1.5701   0.001 ***
# Residual       188  1.68150             
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#############
#Preparing the tables for plotting

#El grafico se puede hacer con base (scaling in) en "species" (variables) or "sites" (row data (plant sp y tto))
# abajo se han sacado los scores RDA1 scaling on sites and RDA1_species scaling on "species"

#Adding the scores for sites (scaling it to sites)
traits_phylogeny<-cbind(traits_phylogeny,
                        scores(model1_5,display = "sites",scaling = "sites",choices = c(1,2))) # c(1,2) une los RDA1 y RDA2 (q son los q mas explican)

#Adding the site scores based on species scaling
temporal<-as.data.frame(
  scores(model1_5,display = "sites",scaling = "species",choices = c(1,2)));names(temporal)<-c("RDA1_species","RDA2_species")
traits_phylogeny<-cbind(traits_phylogeny,temporal) ### uniendo las columnas 


#I could also add the scaled values of the species scores (los scores de mis variables)

###Esto para saber la ubicacion de las variables...Se dividio en 4 para agrandar la escala porq estaban "muy juntos"  
e_vectors_traits<-as.data.frame(
  scores(model1_5,display = "species",choices = c(1,2),scaling = "sites")); #display "species" porque son mis "variables"
e_vectors_traits[1,]<-(e_vectors_traits[1,])/4
e_vectors_traits[2,]<-(e_vectors_traits[2,])/4
e_vectors_traits[3,]<-(e_vectors_traits[3,])/4
e_vectors_traits[4,]<-(e_vectors_traits[4,])/4
e_vectors_traits[5,]<-(e_vectors_traits[5,])/4
e_vectors_traits[6,]<-(e_vectors_traits[6,])/4
e_vectors_traits$traits<-rownames(e_vectors_traits)

str(traits_phylogeny)
str(e_vectors_traits)

#Adding the species with species scaling: (de mis variables)
temporal<-as.data.frame(
  scores(model1_5,display = "species",choices = c(1,2),scaling = "species"));names(temporal)<-c("RDA1_species","RDA2_species")
temporal[1,]<-(temporal[1,])/2
#temporal[2,]<-(temporal[1,])/4
temporal[3,]<-(temporal[3,])/2
temporal[4,]<-(temporal[4,])/2
temporal[5,]<-(temporal[5,])/2
temporal[6,]<-(temporal[6,])/2
e_vectors_traits<-cbind(e_vectors_traits,temporal)

#All in one plot but different color gradient for drought and control
drought_ones<-
  traits_phylogeny %>%group_by(Tto,Sp) %>% ### 
  summarise_at(vars("rootN","rootC","sqrt_root","sqrt_diameter","logRTD","SRSA","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
  filter(Tto=="Drought")
names(drought_ones)[3:7]<-c("Root N", "Root C", "Root mass","Root diameter","RTD") # de la col 3 a la 7 nombrarlos asi..
View(drought_ones)

drought_ones2<-  # no se porq si lo modifico arriba se me corren los nombres de las columnas... por eso he creado este otro
  traits_phylogeny %>%group_by(PlantGroup,Tto,Sp) %>% ### 
  summarise_at(vars("rootN","rootC","sqrt_root","sqrt_diameter","logRTD","SRSA","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
  filter(Tto=="Drought"); names(drought_ones2)[3:7]<-c("Root N", "Root C", "Root mass","Root diameter","RTD")
#View(drought_ones2)

####group by (ordena por Tto y por sp)....es necesario para q luego pueda sacar la media..
# no solo lo ordena..... sino q luego al sacar la media lo hara basado en species por tto...
 # %>% (lo q viene despues del simbolo(la funcion o formula) se hace en la tabla anterior...En el primer caso.. aplique group by en la tabla anterior...)
###  y a esa tabla de media filtarla por el drought....
# names ()  a esa tabla q has creado (drought_ones) ponle estos nombresss 

str(drought_ones)
View(control_ones2)
str(e_vectors_traits)
str(traits_phylogeny)

control_ones<-
  traits_phylogeny %>%group_by(Tto,Sp) %>% 
  summarise_at(vars("rootN","rootC","sqrt_root","sqrt_diameter","logRTD","SRSA","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
  filter(Tto=="Control"); names(control_ones)[3:7]<-c("Root N", "Root C", "Root mass","Root diameter","RTD")

control_ones2<-
  traits_phylogeny %>%group_by(PlantGroup,Tto,Sp) %>% 
  summarise_at(vars("rootN","rootC","sqrt_root","sqrt_diameter","logRTD","SRSA","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
  filter(Tto=="Control"); names(control_ones2)[3:7]<-c("Root N", "Root C", "Root mass","Root diameter","RTD")
#View(control_ones2)

e_vectors_traits$traits<-c("Root N", "Root C", "Root mass","Root diameter","RTD","SRSA")

###########################################
########################   Making the plot

#1 ) RDA plot usando GRADIENT COLOR


#Just adding saving some settings for making the plots nicer looking

my_theme2<-
  theme_bw()+
  theme(title = element_text(size = 25),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.position = "none")
        

ggplot()+
  geom_point(size=2,aes(x=RDA1_species,y=RDA2_species,label=Sp,color=PCoA1),fontface="bold",segment.color = 'transparent',
            data = drought_ones)+
  geom_text(size=4,aes(x=RDA1_species,y=RDA2_species,label=Sp,color=PCoA1),nudge_y = 0.03,fontface="bold",segment.color = 'transparent',
            data = drought_ones)+
  scale_color_continuous(low="red",high="salmon")+
  
  ## Estas lineas para que aparezca el "Drought" and "control" en la grafica
  #geom_text(size=7,aes(x=RDA1_species,y=RDA2_species,label="Drought"),fontface="bold",segment.color = 'transparent',
  #          data = drought_ones %>% group_by(Tto) %>% summarise_at(c("RDA1_species" , "RDA2_species"),mean))+
  #geom_text(size=7,aes(x=RDA1_species,y=RDA2_species,label="Control"),fontface="bold",segment.color = 'transparent',
  #          data = control_ones %>% group_by(Tto) %>% summarise_at(c("RDA1_species" , "RDA2_species"),mean))+
  
  new_scale_color()+
  
  geom_point(size=2,aes(x=RDA1_species,y=RDA2_species,label=Sp,color=PCoA1),fontface="bold.italic",segment.color = 'transparent',
            data = control_ones)+
  geom_text(size=4,aes(x=RDA1_species,y=RDA2_species,label=Sp,color=PCoA1),nudge_y = 0.03,fontface="bold.italic",segment.color = 'transparent',
            data = control_ones)+
  scale_color_continuous(low="darkblue",high="deepskyblue1")+
  geom_text(size=7,aes(x=RDA1_species,y=RDA2_species,label=traits,fontface="bold"),
            hjust = 1.1, vjust =-0.5,segment.size = 0,segment.color = 'transparent',
            data = e_vectors_traits)+
  labs(y="RDA 2 (20.38 %)",x="RDA 1 (22.51%)")+
  my_theme2

str(e_vectors_traits)
#################  SEGUNDO PLOTTTTTTTTTTTTTTTT 

#2) RDA plot usando color per Plant group and different shapes


controldrought<-rbind(control_ones2,drought_ones2)
controldrought$PlantGroup <- as.factor(controldrought$PlantGroup)
controldrought <- unite(controldrought,PlantGroup_Tto, c("PlantGroup","Tto"), sep = "_", remove = FALSE)

controldrought$PlantGroup_Tto = factor(controldrought$PlantGroup_Tto, levels = c("Grasses_Control", "Herbs_Control","Legume_Control",
                                                                               "Grasses_Drought","Herbs_Drought","Legume_Drought")) # order panels


my_theme<-
  theme_bw()+
  theme(title = element_text(size = 25),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        legend.position = "none")

ggplot()+
  geom_point(size=2,aes(x=RDA1_species,y=RDA2_species,shape = PlantGroup_Tto),data=controldrought)+
  scale_shape_manual(values = c(0,1,2,15,16,17))+
  geom_text(size=4,aes(x=RDA1_species,y=RDA2_species,label=Sp),nudge_y = 0.03,#fontface="bold.italic",segment.color = 'transparent',
            data = controldrought)+
  #agregar centroides del tto
  geom_text(size=5,aes(x=RDA1_species,y=RDA2_species,label="Drought"),fontface="bold",segment.color = 'transparent',
            data = drought_ones2 %>% group_by(Tto) %>% summarise_at(c("RDA1_species" , "RDA2_species"),mean))+
  geom_text(size=5,aes(x=RDA1_species,y=RDA2_species,label="Control"),fontface="bold",segment.color = 'transparent',
            data = control_ones2 %>% group_by(Tto) %>% summarise_at(c("RDA1_species" , "RDA2_species"),mean))+
  #agregar traits
  geom_text(size=7,aes(x=RDA1_species,y=RDA2_species,label=traits,fontface="bold"),
            hjust = 1.1, vjust =-0.5,segment.size = 0,segment.color = 'transparent',
            data = e_vectors_traits)+
  labs(y="RDA 2 (20.38 %)",x="RDA 1 (22.51%)")+
  my_theme

###### RDA LEAF TRAITS AND PHYLOGENIA
#Fixing some entries
names(traits_phylogeny)

traits_phylogeny$logSLA<-log10(traits_phylogeny$SLA)
traits_phylogeny$leafC[which(is.na(traits_phylogeny$leafC))]<-NaN
traits_phylogeny$leafN[which(is.na(traits_phylogeny$leafN))]<-NaN
traits_phylogeny$logSLA[which(is.na(traits_phylogeny$logSLA))]<-NaN # donde haya un NA (al tener dato automat aparece NA) q lo ponga como NaN
traits_phylogeny$LDMC[which(is.na(traits_phylogeny$LDMC))]<-NaN     # idem. Leaf C y leaf N por si habian 

#Indeed when I use all the data and take into consideration species, the model is really significant

leaf_var<- traits_phylogeny[-which(traits_phylogeny$spcompleto=="Trifolium_repens"),c("log_shootmass","logSLA","leafC", "leafN","LDMC")] 
# transformando estas mejora un poco los resultados de lo q explica el RDA 

vif <- diag(solve(cor(leaf_var))) # multicolieanility. No mayor a 10 
decorana(leaf_var) # gradient lenght. No mayor a 2 SD 

leaf_var <- rda(leaf_var ~ Tto*spcompleto+Condition(PCoA1,PCoA2,PCoA3,PCoA4), #
                data=traits_phylogeny[-which(traits_phylogeny$spcompleto=="Trifolium_repens"),], scale = TRUE) #podria incluir na.action = na.omit

model_leaf_1_5<-leaf_var # porq de ahi opara abajo todo se hizo llamandolo model_1.5

#model_leaf_1_5<-  EL MODELO COMPLETO 
#  rda(traits_phylogeny[-which(traits_phylogeny$spcompleto=="Trifolium_repens"),c("log_shootmass","logSLA","leafC", "leafN","LDMC")] ~ Tto*spcompleto+Condition(PCoA1,PCoA2,PCoA3,PCoA4), 
#      data=traits_phylogeny[-which(traits_phylogeny$spcompleto=="Trifolium_repens"),], scale = TRUE,na.action=na.omit)

summary(model_leaf_1_5)[["cont"]][["importance"]][,c(1:2)]

#coef(leaf1_5) #canonical coeficients from rda result 
RsquareAdj(model_leaf_1_5)$r.squared #unadjusted and adjusted R^2 retreived from the rda result
RsquareAdj(model_leaf_1_5)$adj.r.squared 
set.seed(1111)
anova.cca(model_leaf_1_5)
set.seed(1111)
anova.cca(model_leaf_1_5,by="term")  # Tb usar by axis 

#############
#Preparing the tables for plotting

#Making a copy of the table (just because it makes easier making the figures)
traits_phylogeny_leaf<-traits_phylogeny
names(traits_phylogeny_leaf)
traits_phylogeny_leaf$RDA1<-NULL     # no se si es necesario traits_phylogeny_leaf no tiene RDA1 
traits_phylogeny_leaf$RDA2<-NULL
traits_phylogeny_leaf$RDA1_species<-NULL
traits_phylogeny_leaf$RDA2_species<-NULL
traits_phylogeny_leaf$rownumber<-rownames(traits_phylogeny_leaf) # creo una nueva columna enumerada del 1 en adelante para luego hacer el "join" 


#Adding scores for sites with site scaling (first) and then species scaling (second)
temporal_leaf<-as.data.frame(
  scores(model_leaf_1_5,display = "sites",scaling = "sites",choices = c(1,2))) # obteniendo unicamente el RDA1 , RDA2
temporal_leaf$rownumber<-rownames(temporal_leaf)  # Create un columna nueva "rownames" que tendra enumerado de 1 en adelante.. 
traits_phylogeny_leaf<-left_join(traits_phylogeny_leaf,temporal_leaf)

#Adding the scaling based on species
temporal_leaf<-as.data.frame(
  scores(model_leaf_1_5,display = "sites",scaling = "species",choices = c(1,2)));names(temporal_leaf)<-c("RDA1_species","RDA2_species")
temporal_leaf$rownumber<-rownames(temporal_leaf)  # crear otra columna nueva como atras

traits_phylogeny_leaf<-left_join(traits_phylogeny_leaf, temporal_leaf) # tengo la nueva tabla con el RDA1,RDA2, RDA1_species, RDA2_species

#Getting the scores for species with site scaling (first) and then species scaling (second)
e_vectors_traits_leaf<-as.data.frame(
  scores(model_leaf_1_5,display = "species",choices = c(1,2),scaling = "sites"));
temporal_leaf<-as.data.frame(
  scores(model_leaf_1_5,display = "species",choices = c(1,2),scaling = "species"));names(temporal_leaf)<-c("RDA1_species","RDA2_species")
temporal_leaf[1,]<-(temporal_leaf[1,])/2  # shoot mass que el valor de RDA1_species lo divida en 2 (para q se vea mejor en el grafico)
#temporal_leaf[2,]<-(temporal_leaf[1,])/4 # Mejor asi.. que sino el SLA queda pegado a shoot
temporal_leaf[3,]<-(temporal_leaf[3,])/2
temporal_leaf[4,]<-(temporal_leaf[4,])/2
temporal_leaf[5,]<-(temporal_leaf[5,])/2

e_vectors_traits_leaf<-cbind(e_vectors_traits_leaf,temporal_leaf)
e_vectors_traits_leaf$traits<-rownames(e_vectors_traits_leaf) # crear una columna llamada "traits" q tenga los nombres de la primera col
#      de e_vectors_traits_leaf (ie, los leaf traits)
e_vectors_traits_leaf$traits<-c("Shoot mass","SLA","Leaf C","Leaf N","LDMC") # esa columna ponerle estos nombres q ya no tiene log

#Plotting  ### Puede dar error.. columna duplicada.. limpiar environment y correr
# solo lo q necesito para este plot ( no lo de root traits)  

ggplot()+
  
  # para que aparezcan los puntos
  geom_point(size=2,aes(x=RDA1_species,y=RDA2_species,label=Sp,color=PCoA1),fontface="bold",segment.color = 'transparent',
            data = traits_phylogeny_leaf %>%group_by(Tto,Sp) %>% 
              summarise_at(vars("log_shootmass","logSLA","leafC", "leafN","LDMC","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
              filter(Tto=="Drought"))+
  
   scale_color_continuous(low="red",high="salmon")+
    #para q aparezca el nombre de la especie
    geom_text(size=4,aes(x=RDA1_species,y=RDA2_species,label=Sp,color=PCoA1),nudge_y = 0.03,fontface="bold",segment.color = 'transparent',
               data = traits_phylogeny_leaf %>%group_by(Tto,Sp) %>% 
                 summarise_at(vars("log_shootmass","logSLA","leafC", "leafN","LDMC","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
                 filter(Tto=="Drought"))+
  
    ## Estas lineas para que aparezca el "Drought" and "control" en la grafica
  #geom_text(size=7,aes(x=RDA1_species,y=RDA2_species,label="Drought"),fontface="bold",segment.color = 'transparent',
  #          data = drought_ones %>% group_by(Tto) %>% summarise_at(c("RDA1_species" , "RDA2_species"),mean))+
  #geom_text(size=7,aes(x=RDA1_species,y=RDA2_species,label="Control"),fontface="bold",segment.color = 'transparent',
  #          data = control_ones %>% group_by(Tto) %>% summarise_at(c("RDA1_species" , "RDA2_species"),mean))+
  
  
  new_scale_color()+
  
  
  geom_point(size=2,aes(x=RDA1_species,y=RDA2_species,label=Sp,color=PCoA1),fontface="bold.italic",segment.color = 'transparent',
            data = traits_phylogeny_leaf %>%group_by(Tto,Sp) %>% 
              summarise_at(vars("log_shootmass","logSLA","leafC", "leafN","LDMC","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
              filter(Tto=="Control"))+
  
  geom_text(size=4,aes(x=RDA1_species,y=RDA2_species,label=Sp,color=PCoA1),nudge_y = 0.03,fontface="bold.italic",segment.color = 'transparent',
             data = traits_phylogeny_leaf %>%group_by(Tto,Sp) %>% 
               summarise_at(vars("log_shootmass","logSLA","leafC", "leafN","LDMC","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
               filter(Tto=="Control"))+
  
  scale_color_continuous(low="darkblue",high="deepskyblue1")+
  geom_text(size=7,aes(x=RDA1_species,y=RDA2_species,label=traits,fontface="bold"),
            hjust = 1.1, vjust =-0.5,segment.size = 0,segment.color = 'transparent',
            data = e_vectors_traits_leaf)+
  labs(y="RDA 2 (16 %)",x="RDA 1 (25%)")+
  xlim(-1.1,2)+
  my_theme


####################

###### RDA PLOT B&N  leaf traits

# Saque el average
control_onesL<-
  traits_phylogeny_leaf %>%group_by(PlantGroup,Tto,Sp) %>% 
  summarise_at(vars("log_shootmass","logSLA","leafC", "leafN","LDMC","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
  filter(Tto=="Control"); names(control_ones_leaf)[1:2]<-c("Shoot mass", "SLA")
View(control_onesL)

drought_onesL<-
  traits_phylogeny_leaf %>%group_by(PlantGroup,Tto,Sp) %>% 
  summarise_at(vars("log_shootmass","logSLA","leafC", "leafN","LDMC","RDA1_species","RDA2_species","PCoA1"),mean) %>% 
  filter(Tto=="Drought"); names(control_ones_leaf)[1:2]<-c("Shoot mass", "SLA")
View(drought_onesL)

controldrought_leaf<-rbind(control_onesL,drought_onesL) # unir filas 

controldrought_leaf<- unite(controldrought_leaf,PlantGroup_Tto, c("PlantGroup","Tto"), sep = "_", remove = FALSE)
controldrought_leaf$PlantGroup_Tto = factor(controldrought_leaf$PlantGroup_Tto, levels = c("Grasses_Control", "Herbs_Control","Legume_Control",
                                                                   "Grasses_Drought","Herbs_Drought","Legume_Drought")) # order para graph
e_vectors_traits_leaf$traits<-c("Shoot mass","SLA","Leaf C","Leaf N","LDMC") # viene de atras

controldrought_leaf<-controldrought_leaf[-which(controldrought_leaf$Sp == "Trifolium"),] # quitar Trifolium 
drought_onesL<-drought_onesL[-which(drought_onesL$Sp == "Trifolium"),]
control_onesL<-control_onesL[-which(control_onesL$Sp == "Trifolium"),]
View(controldrought_leaf)


ggplot()+
  geom_point(size=2,aes(x=RDA1_species,y=RDA2_species,shape = PlantGroup_Tto),data=controldrought_leaf)+
  scale_shape_manual(values = c(0,1,2,15,16,17))+
  geom_text(size=4,aes(x=RDA1_species,y=RDA2_species,label=Sp),nudge_y = 0.03,#fontface="bold.italic",segment.color = 'transparent',
            data = controldrought_leaf)+
  #agregar centroides del tto  ## 
  geom_point(size=5,aes(x=RDA1_species,y=RDA2_species,label="Drought"),fontface="bold",segment.color = 'transparent',
            data = drought_onesL %>% group_by(Tto) %>% summarise_at(c("RDA1_species" , "RDA2_species"),mean))+
  geom_point(size=5,aes(x=RDA1_species,y=RDA2_species,label="Control"),fontface="bold",segment.color = 'transparent',
            data = control_onesL %>% group_by(Tto) %>% summarise_at(c("RDA1_species" , "RDA2_species"),mean))+
  #agregar traits
  geom_text(size=7,aes(x=RDA1_species,y=RDA2_species,label=traits,fontface="bold"),
            hjust = 1.1, vjust =-0.5,segment.size = 0,segment.color = 'transparent',
            data = e_vectors_traits_leaf)+
  labs(y="RDA 2 (16.62 %)",x="RDA 1 (26.02%)")+
  my_theme


###FIN phylogenia RDA
########################################################



#### GRAPH   (e.g Shoot mass 


names(traits)
ggplot(traits, aes(x=Sp, y=leafN, fill=Tto))+  #rlenght0.1cm
  geom_boxplot()+
  geom_point(pch = 21, position = position_jitterdodge())+
  #facet_grid(. ~ PlantGroup)+
  theme_bw()+
  scale_fill_manual(values=c("azure2","darkgoldenrod4","azure2","darkgoldenrod4"))+
  theme(axis.text.x = element_text(face="plain", color="black", size=11, angle=45, hjust = 1),
        axis.text.y = element_text(face="plain", color="black", size=11, angle=0),
        axis.title=element_text(size=13,face="bold"))+
  labs(x='Plant species', y="leaf N", fill="Water treatment")+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  theme(strip.text.x = element_text(size = 13, face = "bold")) #font title panels

anova(lm(leafN ~ Specie*Tto, data =traits ))

str(traits)
################

### Path analyses

traits <- read.csv("traits_25grasslands.csv",stringsAsFactors = F)

traits$Ttobinary <-traits$Tto # copiar la columna Tto
traits$Ttobinary[traits$Tto == "Control"] <- 1 # convertir a binary
traits$Ttobinary[traits$Tto == "Drought"] <- 0 # convertir a binary 


model1.satu<-'               # He usado las mismas variables transf. que en el RDA. 
sqrt_diameter ~ Ttobinary
logRTD ~ Ttobinary
SRSA ~Ttobinary
sqrt_root ~ Ttobinary
rootN ~ Ttobinary
rootC ~ Ttobinary
log_shootmass ~ Ttobinary
leafC ~ Ttobinary
leafN ~ Ttobinary
LDMC ~ Ttobinary
sqrt_diameter~~logRTD
'
fit.f1 <- sem(model1.satu, data=traits)
summary(fit.f1, standardized=TRUE, fit.measures=T,rsq=T) 

fitMeasures(fit.f2b, c("chisq", "df", "pvalue", "cfi", "rmsea", "AIC"))

### Conclusion del path analyses. 
# Solo shoot mass y LDMC are significantly affected by droguht. 

############################################################################
##############################

#### a) Regresiones phylogenia vs traits 

# los traits showing the strongest phylogenetic signal (following Valverde_barrantes 2017 NP)

names(traits_phylogeny)

mytraits<- traits_phylogeny[,c("PCoA1","rootN","rootC","rootmass","diameter","RTD","SRSA",
                               "shootmass","SLA","leafC", "leafN","LDMC")] 

round(cor(x = mytraits, method = "pearson"), 3)

ggpairs(mytraits, lower = list(continuous = "smooth"),
        diag = list(continuous = "bar"), axisLabels = "none")

### regresion tipo II (lmodel2) . Cuando AMBAS variables no son conocidas....
# regresion tipo I (lm) cundo al menos una lo conozco. Por ejemplo YO he determinado el gradiente...lo conozco!!!!

### AQUI estoy comminando los datos de drought y non drought. NO estoy segura si esta bien!!!! # por eso abajo use el RII
lmodel2(rootN ~ PCoA1, data = traits_phylogeny, nperm = 99) #R2 = 0.16 ***
lmodel2(rootC ~ PCoA1, data = traits_phylogeny, nperm = 99) #R2 = 0.01 (0.07)
lmodel2(rootmass ~ PCoA1, data = traits_phylogeny,nperm = 99) #r2 = 0.22 ***
lmodel2(diameter ~ PCoA1, data = traits_phylogeny,nperm = 99) #r2=0.15 ***
lmodel2(RTD ~ PCoA1, data = traits_phylogeny, nperm = 99)  #r2=0.11 ***
lmodel2(SRSA ~ PCoA1, data = traits_phylogeny, nperm = 99) #r2= 0.54 ***

lmodel2(shootmass ~ PCoA1, data = traits_phylogeny, nperm = 99) #r2=0.01 (0.08)
lmodel2(SLA ~ PCoA1, data = traits_phylogeny, nperm = 99) #r2= 0.001 (0.98)
lmodel2(leafC ~ PCoA1, data = traits_phylogeny, nperm = 99) # r2=0.002 (0.43)
lmodel2(leafN ~ PCoA1, data = traits_phylogeny, nperm = 99) #r2=0.31 ***
lmodel2(LDMC ~ PCoA1, data = traits_phylogeny, nperm = 99) #r2=0.22 ***

plot(lmodel2(LDMC ~ PCoA1, data = traits_phylogeny, nperm = 99))

# b) Regresiones phylogenia vs traits RII

dataRII <- read.csv("RII_traits.csv")
names(dataRII)
mytraitsRII<- dataRII[,c("PCoA1","RootN_RII","RootC_RII","Rootmass_RII","diameter_RII","RTD_RII","SRSA_RII",
                               "Shoot_RII","SLA_RII","leafC_RII", "LeafN_RII","LDMC_RII")] 

round(cor(x = mytraitsRII, method = "pearson"), 3)

ggpairs(mytraitsRII, lower = list(continuous = "smooth"),
        diag = list(continuous = "bar"), axisLabels = "none")


lmodel2(RootN_RII ~ PCoA1, data = mytraitsRII, nperm = 99) #R2 = 0.03 (0.03)
lmodel2(RootC_RII ~ PCoA1, data = mytraitsRII, nperm = 99) #R2 = 0.07 (0.002)
lmodel2(Rootmass_RII ~ PCoA1, data = mytraitsRII,nperm = 99) #r2 = 0.002 (0.6)
lmodel2(diameter_RII ~ PCoA1, data = mytraitsRII,nperm = 99) #r2=0.01207819  (0.23)
lmodel2(RTD_RII ~ PCoA1, data = mytraitsRII, nperm = 99)  #r2=0.045 (0.01)
lmodel2(SRSA_RII ~ PCoA1, data = mytraitsRII, nperm = 99) #r2= 0.05 (0.01)

lmodel2(Shoot_RII ~ PCoA1, data = mytraitsRII, nperm = 99) #r2=0.08 (0.001)
lmodel2(SLA_RII ~ PCoA1, data = mytraitsRII, nperm = 99) #r2= 0.000 (0.94)
lmodel2(leafC_RII ~ PCoA1, data = mytraitsRII, nperm = 99) # r2=0.000 (0.97)
lmodel2(LeafN_RII ~ PCoA1, data = mytraitsRII, nperm = 99) #r2=0.000 (0.97) 
lmodel2(LDMC_RII ~ PCoA1, data = mytraitsRII, nperm = 99) #r2=0.000 (0.87)



###########################################################################################
############Graph root traits per section

rsection<- read.csv("roottraits_section.csv")
str(rsection)
Diameter<-lm(Diameter ~ Specie*TTO*SoilSection, data = rsection)
anova(Diameter)
plot(Diameter, add.smooth = FALSE, which = 1)
res.Diameter<-residuals(Diameter) 
qqnorm(res.Diameter)   


diameter<-ggplot(rsection, aes(x=interaction(Specie,TTO), y=Diameter, fill=SoilSection)) +
  geom_boxplot(outlier.size = 15,outlier.colour = "white")+ 
  geom_point(pch = 21, position = position_jitterdodge())+
  #scale_x_discrete(limits=c("well-watered","drought"))+
  theme(axis.title=element_text(size=13,face="bold"), axis.text= element_text(size=11, color="black"))+
  theme(legend.title = element_text(size=13,face="bold"),legend.text = element_text(size=11))+
  theme(legend.position='top',legend.justification='left')+
  labs(x='Species x Water treatment', fill="Soil section")+
  ylab(bquote(bold('Root diameter (mm)')))+
  #scale_fill_manual(values=c("azure2", "darkgoldenrod4"))+
  theme(panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.1))
diameter

RTD<-lm(RTD ~ Specie*TTO*SoilSection, data = rsection)
anova(RTD)
plot(RTD, add.smooth = FALSE, which = 1)
res.RTD<-residuals(RTD) 
qqnorm(res.RTD)   

SRL<-lm(SRL ~ Specie*TTO*SoilSection, data = rsection)
anova(SRL)
plot(SRL, add.smooth = FALSE, which = 1)
res.SRL<-residuals(SRL) 
qqnorm(res.SRL)  

SRSA<-lm(SRSA ~ Specie*TTO*SoilSection, data = rsection)
anova(SRSA)
plot(SRSA, add.smooth = FALSE, which = 1)
res.SRSA<-residuals(SRSA) 
qqnorm(res.SRSA)

##############################################################################

####################### GRAPH....TWO REGRESSION LINES###########
#######################                         ############
###########################################################

data2 <- read.csv("regresion_root_trait_vs_drought-control.csv") 

#he probado los otros: RDA, log RDA y este es el mejor#

RAD_shoot <- ggplot(data2) +
  aes(Sqrt_RDA, log_shootmass, shape = Treatment) +
  geom_point(aes(colour = Treatment), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab("sqrt RAD") +
  ylab("log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 10, 1)) + 
  #scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
RAD_shoot


#he probado con sqrt_RTD y RTD y este es el mejor
RTD_shoot <- ggplot(data2) +
  aes(logRTD, log_shootmass, shape = Treatment) +
  geom_point(aes(colour = Treatment), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab("log RTD") +
  ylab("log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 10, 1)) + 
  #scale_x_continuous(breaks = seq(4, 10))+
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
RTD_shoot

####
SRL_shoot <- ggplot(data2) +
  aes(logSRL3, log_shootmass, shape = Treatment) +
  geom_point(aes(colour = Treatment), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab("log SRL") +
  ylab("log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 10, 1)) + 
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
SRL_shoot

#####

SRSA_shoot <- ggplot(data2) +
  aes(sqrt_SRSA, log_shootmass, shape = Treatment) +
  geom_point(aes(colour = Treatment), size = 3) +
  geom_point(colour = "grey90", size = 1.5) +
  theme_bw() + 
  geom_smooth(method=lm, se=TRUE)+
  xlab("sqrt SRSA") +
  ylab("log Shoot mass") +
  scale_y_continuous(breaks = seq(0 , 10, 1)) + 
  theme(axis.text.x = element_text(size=14, color = "black"),
        axis.text.y = element_text(size=14, color = "black"))+
  theme(
    panel.background = element_rect(fill = "white", colour = "grey50"), # Remove panel background
    panel.border =  element_rect(linetype = "solid", fill = NA, size = 2), # Remove panel border 
    panel.grid.major = element_blank(), # Remove panel grid lines
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black", size=1))+ # Add axis line
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"))
SRSA_shoot

ggarrange(RAD_shoot, RTD_shoot,SRSA_shoot, SRL_shoot,labels = c("A", "B","C","D"),
                       common.legend = TRUE, align=c("hv")) # allign del mismo size horiz y vertic

## Regresiones lineales (tipo II) porq ambas son desconocidas (no es un grandiente q YO defini por ejemplo)

d_drought<- data2[which(data2$Treatment == "drought"),]
d_nondrought<- data2[which(data2$Treatment=='non-drought'),]

lmodel2(Sqrt_RDA ~ log_shootmass, data = d_drought, nperm = 99) # R2 = 0.12 (0.0001)
lmodel2(Sqrt_RDA ~ log_shootmass, data = d_nondrought, nperm = 99) #R2 = 0.04 (0.01)

lmodel2(logRTD ~ log_shootmass, data = d_drought, nperm = 99) #R2 = 0.09 (0.0007)
lmodel2(logRTD ~ log_shootmass, data = d_nondrought, nperm = 99) #R2 = 0.006 (0.39)

lmodel2(sqrt_SRSA ~ log_shootmass, data = d_drought, nperm = 99) #R2 = 0.14 (0.0000001)
lmodel2(sqrt_SRSA ~ log_shootmass, data = d_nondrought, nperm = 99) #r2= 0.02 (0.08)

lmodel2(logSRL3 ~ log_shootmass, data = d_drought, nperm = 99) #R2 = 0.14 (0.0000001)
lmodel2(logSRL3 ~ log_shootmass, data = d_nondrought, nperm = 99) #r2=  0.02 (0.07)


####################################################################
#################################################################
  