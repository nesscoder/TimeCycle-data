
##--------------------------------------- Generating AUC Comparison plot --------------------------------
library("pROC")
library("ggplot2")
library("tidyverse")
library("Rmisc")

##--------------------------------------- Preprocessing TimeCycle Data --------------------------------

#################################
#LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycleV3/Results/Synthetic/TimeCycle/")
all.files <- as.list(dir())
#get Timecycle Results
count <- 0
data <- lapply(all.files,FUN = function(x){
  count <<- count + 1
  load(all.files[[count]])
  attach(all.files[[count]],warn=FALSE)
  df <- get("TimeCycleResults",pos=paste0("file:",all.files[[count]]))
  assign(paste(all.files[[count]], sep = ""), df)
})
names(data) <- all.files

###############################
#EXTRACT DATA FRAME OF INTEREST
###############################
#get the adjusted pvalues and make them into a singel dataframe
TC_results <- lapply(data, function(df){
  df$pVals
})
TC_results  <- do.call(cbind,TC_results )
rownames(TC_results ) <- rownames(data[[1]])
TC_results  <- as.data.frame(TC_results )

###############################
# GET AUC Value OF DATA
###############################
AUC_TimeCycle <- apply(TC_results, 2,FUN = function(exprResult){
  expected <- c(rep(1,7000),rep(0,4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected,pred)
  as.vector(ROC$auc)
} )

##--------------------------------------- Preprocessing JTK_Cycle Data --------------------------------

##########################
#LOAD THE DATA INTO R JTK
##########################
setwd("~/Desktop/TimeCycleV3/Results/Synthetic/JTK/")
all.files <- as.list(dir())
#get Timecycle Results
count <- 0
data <- lapply(all.files,FUN = function(x){
  count <<- count + 1
  load(all.files[[count]])
  attach(all.files[[count]],warn=FALSE)
  df <- get("results",pos=paste0("file:",all.files[[count]]))
  assign(paste(all.files[[count]], sep = ""), df)
})
names(data) <- all.files

###############################
#EXTRACT DATA FRAME OF INTEREST
###############################
#get the adjusted pvalues and make them into a singel dataframe
JTK_results<- lapply(data, function(df){
  df$ADJ.P
})
JTK_results <- do.call(cbind,JTK_results)
rownames(JTK_results) <- rownames(data[[1]])
JTK_results <- as.data.frame(JTK_results)

###############################
# GET AUC Value OF JTK_Results
###############################
AUC_JTK <- apply(JTK_results,2,FUN = function(exprResult){
  expected <- c(rep(1,7000),rep(0,4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected,pred)
  as.vector(ROC$auc)
} )

##------------------------------------------- Preprocessing Sw1pers Data --------------------------------

##########################
#LOAD THE DATA INTO R Sw1pers
##########################
setwd("~/Desktop/TimeCycleV3/Results/Synthetic/Sw1pers/")
all.files <- as.list(dir())
#get Timecycle Results
data <- lapply(all.files,function(x){
  read.delim(x, quote = "")
})
names(data) <- all.files

###############################
#EXTRACT DATA FRAME OF INTEREST
###############################
#get the adjusted pvalues and make them into a singel dataframe
sw1per_results<- lapply(data, function(df){
  df$score
})
sw1per_results <- do.call(cbind,sw1per_results)
rownames(sw1per_results) <- data[[1]]$id
sw1per_results <- as.data.frame(sw1per_results)

###############################
# GET AUC Value OF JTK_Results
###############################
AUC_sw1pers <- apply(sw1per_results,2,FUN = function(exprResult){
  expected <- c(rep(1,7000),rep(0,4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected,pred)
  as.vector(ROC$auc)
} )

#################################
#LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycleV3/Results/Synthetic/RAIN/")
all.files <- as.list(dir())
#get RAIN Results
count <- 0
data <- lapply(all.files,FUN = function(x){
  count <<- count + 1
  load(all.files[[count]])
  attach(all.files[[count]],warn=FALSE)
  df <- get("results",pos=paste0("file:",all.files[[count]]))
  assign(paste(all.files[[count]], sep = ""), df)
})
names(data) <- all.files

###############################
#EXTRACT DATA FRAME OF INTEREST
###############################
#get the adjusted pvalues and make them into a singel dataframe
RAIN_results <- lapply(data, function(df){
  df$pVal
})
RAIN_results  <- do.call(cbind,RAIN_results )
rownames(RAIN_results ) <- rownames(data[[1]])
RAIN_results  <- as.data.frame(RAIN_results)

###############################
# GET AUC Value OF DATA
###############################
AUC_RAIN <- apply(RAIN_results, 2,FUN = function(exprResult){
  expected <- c(rep(1,7000),rep(0,4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected,pred)
  as.vector(ROC$auc)
} )


##------------------------------------------- META DATA --------------------------------

##################################################
#LOAD META DATA INTO R SAME FOR JTK, TIMECYCLE, and Sw1pers
##################################################

meltDataForPlotting <- function(AUCscores, method){
  exprLabels <- names(AUCscores)
  metaInfo <- regmatches(exprLabels, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", exprLabels))
  metaInfo <- as.data.frame(do.call(rbind,metaInfo))
  metaInfo <- as.data.frame(apply(metaInfo, 2, function(x){ as.numeric(as.character(x))}))
  colnames(metaInfo) <- c("Interval","Length","NoiseLevel","Replicates")
  metaInfo$Method <- rep(method, length(AUCscores))
  metaInfo$AUC <- as.numeric(unname(unlist(AUCscores)))
  metaInfo$numSamples <- (metaInfo$Length/metaInfo$Interval+1)*metaInfo$Replicates
  return(metaInfo)
}



##--------------------------------------- Plot Results --------------------------------
#select points to be plotted
JTK_CYCLE  <- meltDataForPlotting(AUC_JTK, "JTK_CYCLE")
RAIN <- meltDataForPlotting(AUC_RAIN, "RAIN")
TimeCycle <- meltDataForPlotting(AUC_TimeCycle, "TimeCycle")
SW1PERs <- meltDataForPlotting(AUC_sw1pers, "SW1PERs")

toPlot <- rbind(JTK_CYCLE,TimeCycle, RAIN, SW1PERs)


#----------------

pdf("~/Desktop/TimeCycleV3/Results/Figures/AUCcompPlots.pdf",width = 7*1.5,height = 3.5*1.5)

toPlot <- as.data.frame(toPlot) %>%
  filter(Length != 24)

theme_set(new = theme_light())
theme_replace(panel.border = element_rect(colour = "black", fill=NA,size =2),
              text = element_text(size = 16, face = "bold"))

ggplot(toPlot, 
       aes(x = as.numeric(numSamples), 
           y = as.numeric(AUC), 
           fill = Method, 
           colour = Method)) +
  geom_point()+
  geom_point(shape = 21, colour = "black", show.legend = F) +
  geom_smooth(method = lm, se=FALSE, size = 1) +
  scale_fill_manual(labels = c("JTK_CYCLE","RAIN","SW1PERS","TimeCycle"),
                     values = c("#377eb8", "black","grey","#ff7f00"))+
  scale_colour_manual(labels = c("JTK_CYCLE","RAIN","SW1PERS","TimeCycle"),
                    values = c("#377eb8", "black","grey","#ff7f00"))+
  facet_grid(Interval ~ Length)+
  xlab("Number of Samples")+
  ylab("AUC") +
  guides(colour = guide_legend(override.aes = list(linetype = 0, shape = 22, size = 5)), fill = guide_legend(nrow=1,byrow=TRUE)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom', legend.title = element_blank())

dev.off()



sdToPlot <- toPlot %>%
  summarySE(measurevar="AUC", groupvars=c("Method","Length","Interval","Replicates")) %>%
  mutate(numSamples = (Length/Interval+1)*Replicates)

ggplot(sdToPlot, 
       aes(x = as.numeric(numSamples), 
           y = as.numeric(AUC), 
           fill = Method, 
           colour = Method)) +
  geom_point()+
  geom_smooth(method = lm, se=FALSE, size = 1) +
  geom_errorbar(aes(ymin=AUC-se, ymax=AUC+se), width = 1) +
  geom_point(shape = 21, colour = "black", show.legend = F) +
  scale_fill_manual(labels = c("JTK_CYCLE","RAIN","SW1PERS","TimeCycle"),
                    values = c("#377eb8", "black","grey","#ff7f00"))+
  scale_colour_manual(labels = c("JTK_CYCLE","RAIN","SW1PERS","TimeCycle"),
                      values = c("#377eb8", "black","grey","#ff7f00"))+
  facet_grid(Interval ~ Length)+
  xlab("Number of Samples")+
  ylab("AUC") +
  guides(colour = guide_legend(override.aes = list(linetype = 0, shape = 22, size = 5)), fill = guide_legend(nrow=1,byrow=TRUE)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom', legend.title = element_blank())



