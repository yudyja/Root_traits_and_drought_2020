install.packages("phytools")
install.packages("picante")


library(tidyverse)
library(phytools)
library(picante)

#Phylogenetic tree#

#PIC. DA ERROR PORQ TENGO VARIAS REPLICAS POR ESPECIE. SI TUVIERA UN SOLO DATO 
#POR ESPECIE..HARIA MATCH CON EL TREE (PORQ SOLO HAY UNA DATO POR TREE)
#ESTE PIC SE USA MUCHO SI QUIERO COMPARAR 2 TRAITS. TENGO LA MEDIA DEL TRAIT Y LUEGO LE APLICO ESTO DE LA FILOGENIA. 

#0. #Identificar directorio#
getwd()

#1.2. Paquetes and libraries.  


#3. Generar el arbol filogenetico.

#3.1. Opcion 1. Con este archivo. Pero los nombres salieron como "@". DIO PROBLEMAS!
#https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
#3.1.1, Cargarlo asi. tree25_sp<-read.tree("phyliptree_25Grasslands.phy"). Dio problemas, por el @
#plot(tree). se ve lo del @

#3.2. Opcion 2. De la filogenia de arboles de Europa. Armar mi arbol. 
#3.2.1. Descargar el TreeEurope
#https://figshare.com/collections/Daphne_a_dated_phylogeny_of_a_large_European_flora_for_phylogenetically_informed_ecological_analyses/3305040

#3.2.3 Leer el arbol de Europe
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


#8. cargar el trait 

# ESTO CREO Q ES LO MISMO QUE ABAJO!
shoot_mass<-traits[,21]#aqui escoges en la tabla de datos el trait que quieres
names(shoot_mass)<-traits[,4]#aqui escoges el nombre de la especie y se la agregas como "nombre"

#9. Media por trait por especie # porq en el tree solo hay un dato por specie. 

mean_sptraits <- traits %>% group_by(spcompleto) %>% summarise_at(c("shootmass", "rootmass","diameter", "RTD", "SRL", "SRSA", 
                                          "rootN", "rootC"), mean) # summarise or summarize
#meanshoot <- mean_sptraits[,2] # a veces da problemas
meanshoot <- mean_sptraits$shootmass            #choose one trait from your table
names(meanshoot) <- mean_sptraits$spcompleto    #choose the column with species names (have to be written exactly like in your tree)
match.phylo.data(WP4_PhyloTree, meanshoot)     #bring traits & names and the tree together
phylosig(WP4_PhyloTree, meanshoot, method="lambda", test=TRUE) ##important values: K and p-value, 
                                          #or lambda si escojo ese metodo. See Aguilar+trigeros 2018 ISME 
traits$shootmass[which(traits$codigo=="Armeria_maritima_259m")]


shoot_phy<-pic(traits$shoot_mean,WP4_PhyloTree)


#10. Hallar el valor medio 

###########################################################################
############################################################################
########## MAKING THE PICS FOR EACH VARIABLE AND FOR EACH REPLICATE ########
############################################################################
############################################################################


rm(list=ls())

library(ape)
library(picante)
library(phytools)


plot(WP4_PhyloTree)

traits <- read.csv("traits_25grasslands.csv",stringsAsFactors = F) #a veces es mejor para 
# q al importarlo no sea factor...
names(traits)

#crear nueva columna species nombre completo

traits$spcompleto <- NA # nueva columna empty

traits$spcompleto[which(traits$Sp=="Anthoxantum")]<-"Anthoxanthum odoratum"
traits$spcompleto[which(traits$Sp=="Berteroa")]<-"Berteroa incana"
traits$spcompleto[which(traits$Sp=="Daucus")]<-"Daucus carota"
traits$spcompleto[which(traits$Sp=="F.brevipila")]<-"Festuca brevipila"
traits$spcompleto[which(traits$Sp=="Achillea")]<-"Achillea millefolium"
traits$spcompleto[which(traits$Sp=="Potentilla")]<-"Potentilla argentea"
traits$spcompleto[which(traits$Sp=="Plantago")]<-"Plantago lanceolata"
traits$spcompleto[which(traits$Sp=="Hieracium")]<-"Hieracium pilosella"
traits$spcompleto[which(traits$Sp=="Artemisia")]<-"Artemisia campestris"
traits$spcompleto[which(traits$Sp=="Holcus")]<-"Holcus lanatus"
traits$spcompleto[which(traits$Sp=="Arrhenaterum")]<-"Arrhenatherum elatius"
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

traits$spcompleto <- gsub(" ","_",traits$spcompleto) #substituir!

#Organizing the table in alphabetical order
traits<-traits[order(traits$spcompleto),] # porq luego R hace cosas orden alfabetico.. entonces mejor

#Creating a number sytem so each pot per species gets a value from 1 to 10 (or 9)
numero<-split(traits$PotN,traits$spcompleto) 
numero<-sapply(numero,function(x){number<-as.matrix(seq(1:length(x)))})#apply. 
#aplica una funcion (x) se la invento... sapply(lo hace en cada specie independiente)
numero<-as.numeric(do.call("rbind",numero)) #q sean numeros. do call (se ordena) que todo \
# se una por fila.

#Adding this number to the data and addit to it a letter so it can be used as a factor later
traits$numero<-paste("R",numero,sep = "") 


#Preparing everything for the corrections

#First, splitting the table by "numero". This makes a list of 10 dataframes each containing only
#one replicate per species (that is, 24 rows), except for one datatframe that has only 20 rows 
#because in some plants there are only 9 replicates per plant instead of 10.
prueba<-split(traits,traits$numero) # tabla descopuesta en 10 tablas..conteien una replica por SP

#The names of the plants are orderd allways in the same way, so I create a vector with those names
nombres<-prueba$R1$spcompleto # se necesita abajo.

#Select the traits to be corrected
traits_to_correct<-
  lapply(prueba,function(tabla){tabla<-tabla[,c(10:13,20:38)]})
# lapply es para un lista. Sapply (trata de simplificar.. creo o matriz o vector..segun crea
# es mas simple.. pero a veces eso q quiere a mi no me sirve)
# aqui "prueba"  es una lista por eso he usado lapply

#Correct the traits for each replicate
traits_corrected<-
  lapply(traits_to_correct[-2],function(a){#Problem 1: I remove the dataframe "R10" that has not all the 24 sps but only 10
    
      lapply(a,function(b){
      trait<-b
      names(trait)<-nombres #
      if(!any(is.na(b))){#Problem 2: Pic does not work when there are NA´s in the data, so I remove those
        trait_phylo<-pic(trait,WP4_PhyloTree)}
    })
    
  })
#excluir la tabla numero 2 ( q es la R10 q no tiene las 24 species
# function a ..hazlo en cada tabla. FUncion b.. hazlo en cada columna (trait)
#lapply(dataframe, mean)
#"a" aqui es una data frame... a cada una de las tablas (R1)
# b es un trait...
# b es un trait...(q en realidad es una columna) ..es cada columna de traits.. (entonces
# asignarle el nombre de las species.. porq tengo 24 datos..pero no sabe de q son)

# 


#Transforming the output into a list of dataframes # en lugar de tener 10 tablas tener una

traits_corrected<-
  lapply(traits_corrected,
         function(x){
           as.data.frame(x[!sapply(x,is.null)])
           #x<-do.call("cbind",x)
         })


#Putting all in just one big dataframe

traits_corrected<-do.call("rbind",traits_corrected)


#In summary what are the problems:

#1. I do not know which plant is which yet after the transformations...
# I mean, after transforming, I do not get the names of the species, instead I get things called
# "Caryiophyllales", "asterids_I", "N5136". I guess I do not know of big of an issue this is

#2. Pic function does not work with NA´s

#3. What I did was dividing the dataset into 10 tables of 24, so pic can done on each of those tables
#The problem is that when there are no 24 species, pic can be done but it would not be the same as the
#in the others because some species will be missing.


########################################################################################################

#Option 2!

#La idea de este script es crear una matrix de disimilitud entre las plantas que repesente su relación 
#filogenética. De esta infromarción se hace un Principal Coordinate Analysis  y luego se toman los eigen
#vectors de esta ordinación como una co-variate de la relación entre plants

library(ape)

#1. Creando una matrix de dismilitud entre las plantas
dist.tree <- cophenetic(WP4_PhyloTree)

#2. Haciendo el PCoA de esta disimilitud
pcoa.dist <- cmdscale(dist.tree, k= dim(dist.tree)[1]-1, eig=T)
#k is the maximum dimension of the space which the data are to be represented in; must be in {1, 2, …, n-1}.


#3.1. Identificar  cuales seran los ejes eigenvalues que se usaran. Stefan uso (arbitrariarmente) com criterio 
# eigenvalues que en suma explican 80% de la variación de los datos.
phyl.var1 <- 
  1:length(which(cumsum(pcoa.dist[[2]]/sum(pcoa.dist[[2]]))<0.8)) 


#3.2. Seleccionar los eigenvectors (loadings) que corresponden con los eigenvalues identfiicados
phylogenieachsen<-pcoa.dist[[1]][,phyl.var1]


#4. Listo! Ahora es de decidir cuales loadings se usara de los cuatro. Lo más lógico es usar el primero



#Version de Jeff

#Es lo mismo solo que usa vegan y no ape

library(vegan)

#1. Dissimilarity matrix
tree.dist<-as.dist(cophenetic(WP4_PhyloTree))  # don’t forget to convert to class ‘dist’, otherwise treats it as raw data and runs it through vegdist()

#2. MDS on based on this dissimarity matrx
tree.pco<-capscale(tree.dist~1)  


#3. Selecting eigenvectors
# principle coordinates matrix is in scores(tree.pco,1:8,‘sites’); Why 1:8? because they are the ones that explain the most?
tree.scores<-scores(tree.pco,1:8,"sites")

# Use one of the ‘join’ functions in dplyr to merge with the rest of your data instead of the following
phylo.pco<-matrix(NA,nrow=nrow(traits),ncol=ncol(tree.scores))
colnames(phylo.pco)<-paste('phylo',colnames(tree.scores),sep='')
for(i in rownames(tree.scores)) {
  x<-tree.scores[i,]
  xx<-which(traits$spcompleto==i)
  for(j in xx) phylo.pco[j,]<-x
}
rm(i,j,x,xx)

#This seems like general advice on which explanatory variables to use.
# The code I shared here used forward.sel() for deciding which variables to use, which I think is depreciated. 
#Use ordistep(), you can find an example in the ‘analysis of microbial communities’ document for R.


#Versión de amiga de Yudi

#Parece ser una versión muy similar a la de Stefan, en vey de usar ape usa hierfstat

# Takes a trait (data frame with one column) and a phylogenetic tree 
# (in ape::phylo format) and caclulates residuals from regressing the trait
# on the tree via scores of an PCoA of the phylogenetic distance of the tree.
phylo_resids <- function(trait, tree) {
  md <- match_data(tree = tree, traits = trait, scale = FALSE)
  pcoa_tree <- hierfstat::pcoa(cophenetic(md$tree))#1. Matrix de diferencias entre plantas y hacer PcoA de eso
  ind <- cumsum(pcoa_tree$valp)/sum(pcoa_tree$valp) < 0.8# 2-3. Seleccionar los eigenvectors que expliquen más del 80% de la variación
  df <- data.frame(fitcomp = md$traits[, 1], pcoa_tree$vecp[, ind])#4. Poner todo junto
  rownames(df) <- rownames(md$traits)
  phylo_resids <- resid(step(lm(fitcomp ~ . , data = df), trace = FALSE))#Aquí creo que es lo que ha hecho lo de: "#tienes q guardar los residuos del modelo  que ya has ehcho y eso es lo que vas a testar usando el PCoA con la filogenia"
  return(list(phylo_residuals = phylo_resids, tree = md$tree))
}

#Sugerencia

#I. Obteniendo los Eigenvectors siguiendo el ejemplo de Stefan

#1. Obteniendo la matriz de disimilutd
dist.tree <- cophenetic(WP4_PhyloTree)

#2. Haciendo el PCoA de esta disimilitud
pcoa.dist <- cmdscale(dist.tree, k= dim(dist.tree)[1]-1, eig=T) #   son 24 species... y solo son posible 23 comparaciones
#k is the maximum dimension of the space which the data are to be represented in; must be in {1, 2, …, n-1}.


#3.1. Identificar  cuales seran los ejes eigenvalues que se usaran. Stefan uso (arbitrariarmente) com criterio 
# eigenvalues que en suma explican 80% de la variación de los datos.
phyl.var1 <- 
  1:length(which(cumsum(pcoa.dist$eig/sum(pcoa.dist$eig))<0.8)) 
# cumulative sum..de los eigenvalues... / total...(cuales eigen values explican hasta el 80%)

# eigen value ~ varianza... eigenvector= coeficientes con los que se multiplican mi raw data para transformarlo en el PCA? space

#3.2. Seleccionar los eigenvectors (loadings) que corresponden con los eigenvalues identfiicados
phylogeny_vectors<-pcoa.dist$points[,phyl.var1] # escoger del 1 al 4 (phyl.var)
phylogeny_vectors<-as.data.frame(phylogeny_vectors)
names(phylogeny_vectors)<-c("PCoA1","PCoA2","PCoA3","PCoA4") # ponerle nombres.. PCoA
phylogeny_vectors$spcompleto<-rownames(phylogeny_vectors) #crear columna..con el nombre de las especies..
# la primera q parece una columna no lo es... es el nombre de las filas...

#Ordering them in the same way as the traits database

library(tidyverse)

phylogeny_vectors<-  # unir tabla traits y PcoA results table
  left_join(traits,phylogeny_vectors,by="spcompleto")#[,c("PCoA1","PCoA2","PCoA3","PCoA4")]

#el pCoA solo tiene 24 datos... y con el left.join...se a;ade cada uno de estos valores tantas 
# veces como aparezca la specie.

#II. Selecting which eigenvectors to use (from Jeff)

#I think this is not necessary.
# generate a new dataframe containing scaled predictor variables
# since difficult to do this in the formula
#varechem.scld <- decostand(varechem, method='standardize')

traits_to_correct<-
  traits[,c("rootN","rootC","log_shootmass","shootmass","sqrt_root","rootmass","sqrt_R.S","log_R.S.3",
            "sqrt_diameter","diameter","X.fineRoots_RootSamp","RootVolume.cm3._RootSamp",
            "logRTD","RTD","SRL","logSRL3","sqrt_SRSA","SRSA","logSRSA.3")]


# have to use the formula interface so generate a new 'rda' object
# including all predictors (use '.' after the '~')

vare.rda <- rda(traits_to_correct~ ., data=phylogeny_vectors, scale=T)

# set up the null case with no predictors (be sure to include the
# 'data' argument, even though no predictors)
vare.pca <- rda(traits_to_correct ~ 1, data=phylogeny_vectors, scale=T)

# select variables in each predictor table (output not shown)
step.env <- ordistep(vare.pca, scope=formula(vare.rda))
step.env
#This indicates that all 4 PCo4 are meaningful

# #From Jeff´s example in the book
# #generate a new dataframe containing scaled predictor variables
# # since difficult to do this in the formula
# varechem.scld <- decostand(varechem, method='standardize')
# # have to use the formula interface so generate a new 'rda' object
# # including all predictors (use '.' after the '~')
# vare.rda <- rda(varespec ~ ., data=varechem.scld, scale=T)
# # set up the null case with no predictors (be sure to include the
# # 'data' argument, even though no predictors)
# vare.pca <- rda(varespec ~ 1, data=varechem.scld, scale=T)
# # select variables in each predictor table (output not shown)
# step.env <- ordistep(vare.pca, scope=formula(vare.rda))
# step.env


#III. Doing a regression of the residuals of the original model (traits~drought) to determine whether phylogeny predicts something

#The other alternative is doing the other way around, where we do the analysis traits~drought after correcting for phylogeny
#or having everything in the same model and using phylogeny as a covariate.


#Let´s discuss what to do!






