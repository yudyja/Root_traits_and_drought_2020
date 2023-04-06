#working directory and data loading
getwd()
setwd("E:/Unikrams/1. Master/1.Masterarbeit/R")
hypha<- read.csv2("hyphallengthforr.csv")
options(max.print=1000000)
hypha #or head(hypha) for first 6 lines of table
head(hypha)

#visualize just for a first impression of dataset
par(mar=c(15,4,4,2))
plot(AMF~species, data=hypha, las=2, xlab=NULL) #indicates that variances are not homogenous(no equal size of boxes)
plot(nonAMF~species, data=hypha, las=2, xlab=NULL) #indicates that variances are not homogenous(no equal size of boxes)
plot(AMF~tto, data=hypha)
plot(nonAMF~tto, data=hypha)
par(op)

                                   #####statistics for AMF######
####apply linear model

###for dependent variable AMF###
model1= aov(AMF ~ species+tto + species:tto, data=hypha)
summary(model1)# not necessary here, can/should be done later

##visual inspection
jpeg(file="amfplotresid.jpeg",height=480,width=480)
par(mfrow=c(2,2))
plot(model1) #residuals vs fitted make a pattern, qq is okay
dev.off()
par(mfrow=c(1,1))


#log-transformation
jpeg(file="amflogplotresid.jpeg",height=480,width=480)
model2= aov(log(AMF) ~ species+tto + species:tto, data=hypha)
par(mfrow=c(2,2))
plot(model2) # residual vs fitted looks better, normal qq is good as well
dev.off()


####-> decision: take log-transformed data for AMF

##anova
aovAMF<-aov(log(AMF) ~ species+tto + species:tto, data=hypha)# 
summary(aovAMF) #species, tto are significant but not their interaction
capture.output(summary(aovAMF),file="anovaAMFlog.doc") 

##TukeyHSD-Test
options(max.print=1000000)
TukeyHSD(aovAMF) # 
capture.output(TukeyHSD(aovAMF), file="tukeyamflog.txt")

               
               
               
               
               
##still working on that. ignore it.##               
model.tables(aovAMF)
model.tables(aovAMF, "means") #creates a table of means, seperated by species/tto(control/drought)and also shows repetions(no of samples)
##



##boxplot##

###das funktioniert schonmal aber muss ich noch anpassen
ggsave("AMFboxplot.png", width = 5, height = 5) 
  ggplot(hypha, aes(species, AMF)) + 
    geom_boxplot(aes(fill = tto)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0,5, hjust=1)) +
    xlab("Species")+
    ylab("Hyphal length of AMF")+
    scale_fill_discrete(name="Treatment")
  
  


# barplot
ggsave("AMFbarplot.png", width = 5, height = 5) 
  ggplot(data=hypha, aes(x=species, y=AMF, fill=tto)) + 
    geom_bar(stat="identity", position=position_dodge())+
    theme(axis.text.x = element_text(angle = 90, vjust = 0,5, hjust=1))+
    scale_fill_manual(values = c("cadetblue", "firebrick"))+
    xlab("Species")+
    ylab("Hyphal length of AMF")+
    labs(fill="Treatment")+
    geom_errorbar(aes(ymin=AMF-sd, ymax=AMF+sd), width=.2,
                  position=position_dodge(0.9))
  
  
  
###barplot###
  
ggsave("amfbarplot.png", width = 8, height = 5) 
  amfbarchart = ggplot(hypha, aes(species, AMF, fill=tto))
  amfbarchart+
    stat_summary(fun.y=mean, geom="bar", position="dodge")+
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, position=position_dodge(width=0.90))+
    scale_fill_manual(name="Treatment", values = c("cadetblue", "firebrick"))+
    xlab("Species")+
    ylab("Hyphal length of AMF")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0,5, hjust=1))
  




  
  
                         #####statistics for nonAMF######  
  
  ####apply linear model
  
  ###for dependent variable nonAMF###
  model3= aov(nonAMF ~ species+tto + species:tto, data=hypha)
  summary(model3)# not necessary here, can/should be done later
  
  ##visual inspection
  jpeg(file="nonamfplotresid.jpeg",height=480,width=480)
  par(mfrow=c(2,2))
  plot(model3) #residuals vs fitted make a pattern, qq is okay
  dev.off()
  par(mfrow=c(1,1))
  
  
  #log-transformation
  jpeg(file="nonamflogplotresid.jpeg",height=480,width=480)
  model4= aov(log(nonAMF+1) ~ species+tto + species:tto, data=hypha) #+1 because of some zero-Values
  par(mfrow=c(2,2))
  plot(model4) # residual vs fitted still has a pattern, normal qq is good as well
  dev.off()
  
  #sqrt-transformartion
  jpeg(file="nonamfsqrtplotresid.jpeg",height=480,width=480)
  model5= aov(sqrt(nonAMF) ~ species+tto + species:tto, data=hypha)
  par(mfrow=c(2,2))
  plot(model5) # residual vs fitted still has a pattern, normal qq is good as well
  dev.off()
  
  #log10-transformation
  jpeg(file="nonamflog10plotresid.jpeg",height=480,width=480)
  model6= aov(log10(nonAMF+1) ~ species+tto + species:tto, data=hypha) #+1 because of some zero-Values
  par(mfrow=c(2,2))
  plot(model6) # residual vs fitted still has a pattern, normal qq is good as well
  dev.off()
  
  #cuberoot-transformation does not work like this!! Don't know how to perform this in R
  jpeg(file="nonamfcuberootplotresid.jpeg",height=480,width=480)
  model7= aov(cuberoot(nonAMF) ~ species+tto + species:tto, data=hypha)
  par(mfrow=c(2,2))
  plot(model7) # 
  dev.off()
  
  #reciprocal-transformation does not work like this!! Don't know how to perform this in R
  jpeg(file="nonamfreciprocalplotresid.jpeg",height=480,width=480)
  model8= aov((1/nonAMF) ~ species+tto + species:tto, data=hypha)
  par(mfrow=c(2,2))
  plot(model8) # 
  dev.off()
  
  
  ####-> decision: take sqrt-transformed data for nonAMF
  
  ##anova
  aovnonAMF<-aov(sqrt(nonAMF) ~ species+tto + species:tto, data=hypha)
  summary(aovnonAMF) #species, tto are significant but not their interaction
  capture.output(summary(aovnonAMF),file="anovanonAMFlog.doc") 
  
  ##TukeyHSD-Test
  options(max.print=1000000)
  TukeyHSD(aovnonAMF) 
  capture.output(TukeyHSD(aovnonAMF), file="tukeynonamfsqrt.txt")
 
  
  ###barplot###
  
  ggsave("nonamfbarplot.png", width = 8, height = 5) 
  nonamfbarchart = ggplot(hypha, aes(species, nonAMF, fill=tto))
  nonamfbarchart+
    stat_summary(fun.y=mean, geom="bar", position="dodge")+
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, position=position_dodge(width=0.90))+
    scale_fill_manual(name="Treatment", values = c("cadetblue", "firebrick"))+
    xlab("Species")+
    ylab("Hyphal length of nonAMF")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0,5, hjust=1))
  
  
  
  
                           ##### statistics for total_hyphae ######  
  
  ####apply linear model
  
  ###for dependent variable total_hyphae###
  model9= aov(total_hyphae ~ species+tto + species:tto, data=hypha)
  summary(model9)# not necessary here, can/should be done later
  
  ##visual inspection
  jpeg(file="total_hyphaeplotresid.jpeg",height=480,width=480)
  par(mfrow=c(2,2))
  plot(model9) #residuals vs fitted make a pattern, qq is good
  dev.off()
  par(mfrow=c(1,1))
  
  
  #log-transformation
  jpeg(file="total_hyphaelogplotresid.jpeg",height=480,width=480)
  model10= aov(log(total_hyphae) ~ species+tto + species:tto, data=hypha) #+1 because of some zero-Values
  par(mfrow=c(2,2))
  plot(model10) # residual vs fitted still has a pattern, normal qq is good as well
  dev.off()
  
  #sqrt-transformartion
  jpeg(file="total_hyphaesqrtplotresid.jpeg",height=480,width=480)
  model11= aov(sqrt(total_hyphae) ~ species+tto + species:tto, data=hypha)
  par(mfrow=c(2,2))
  plot(model11) # residual vs fitted still look good, normal qq is good as well
  dev.off()
  
  ####-> decision: take log-transformed data for total_hyphae
  
  ##anova
  aovtotal_hyphae<-aov(log(total_hyphae) ~ species+tto + species:tto, data=hypha)
  summary(aovtotal_hyphae) #species, tto are significant but not their interaction
  capture.output(summary(aovtotal_hyphae),file="anovatotal_hyphaelog.doc") 
  
  ##TukeyHSD-Test
  options(max.print=1000000)
  TukeyHSD(aovtotal_hyphae) 
  capture.output(TukeyHSD(aovtotal_hyphae), file="tukeytotal_hyphaelog.txt")
  
  
  ###barplot###
  
  ggsave("totalhyphaebarplot.png", width = 8, height = 5) 
  totalhyphaebarchart = ggplot(hypha, aes(species, total_hyphae, fill=tto))
  totalhyphaebarchart+
    stat_summary(fun.y=mean, geom="bar", position="dodge")+
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.2, position=position_dodge(width=0.90))+
    scale_fill_manual(name="Treatment", values = c("cadetblue", "firebrick"))+
    xlab("Species")+
    ylab("Hyphal length of total fungi")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0,5, hjust=1))
  