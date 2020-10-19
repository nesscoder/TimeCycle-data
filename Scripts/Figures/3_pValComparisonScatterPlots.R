##--------------------------------------- Generating pVal Comparison plot --------------------------------
library("pROC")
library("ggplot2")
library("gridExtra")
library("tidyverse")
library("rlang")

##--------------------------------------- Preprocessing TimeCycle Data --------------------------------

#################################
#LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycleV3/Results/Synthetic/TimeCycle_comp/")
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

##--------------------------------------- Preprocessing JTK Data --------------------------------

#################################
#LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycleV3/Results/Synthetic/JTK/")
all.files <- as.list(dir())
#get JTK Results
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
JTK_results <- lapply(data, function(df){
  df$ADJ.P
})
JTK_results  <- do.call(cbind,JTK_results )
rownames(JTK_results ) <- rownames(data[[1]])
JTK_results  <- as.data.frame(JTK_results )

##--------------------------------------- Preprocessing RAIN Data --------------------------------

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

##--------------------------------------- Preprocessing Sw1pers Data --------------------------------

##########################
#LOAD THE DATA INTO R Sw1pers
##########################
setwd("~/Desktop/TimeCycleV3/Results/Synthetic/Sw1pers/")
all.files <- as.list(dir())
#get Timecycle Results
data <- lapply(all.files,function(x){
  read.delim(x, quote = "")
})
names(data) <- gsub(pattern = "_sw_scores.txt", replacement = ".Rdata",x = all.files)

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



#keep only what is needed for plotting
rm(list = setdiff(ls(), c("JTK_results","RAIN_results","TC_results", "sw1per_results")))


##--------------------------------------- Plot Functions --------------------------------

best_ROC <- function(data){
  expected <- c(rep(1,7000),rep(0,4000))
  pred <- data
  ROC <- roc(expected,pred)
  out <- coords(ROC, "best", ret="threshold", transpose = F, best.method = "youden")
  as.vector(unlist(out))
}

plot_scatter <- function(data, xAxis = "Method 1", yAxis = "Method 2", coloring = "wave", xAxisLab = NULL, yAxisLab = NULL, waveFormCols = waveFormCols){
  
  #log Transform the data
  data <- -log(data)
  
  wave <- gsub(pattern = "_.*",replacement = "",rownames(data))
  data$wave <- wave
  data$wave <- fct_relevel(wave, ordering)
  
  x <- enquo(xAxis)
  y <- enquo(yAxis)
  
  data_mean <- data %>% 
     group_by(wave) %>% 
     summarise(x = mean(!! x),
               y  = mean(!! y))
  
  if(is.null(xAxisLab)){
    xAxisLab = xAxis
  }
  if(is.null(yAxisLab)){
    yAxisLab = yAxis
  }
  
  ggplot(data, aes_string(x=xAxis, y=yAxis, fill = coloring, color = coloring)) +
    geom_point(alpha = .4, size = 0.1) +
    geom_point(data = data_mean, aes(x = x, y = y), size = 4, shape = 21, color = "black", alpha = .9) + 
    scale_color_manual(values = waveFormCols) +
    scale_fill_manual(values = waveFormCols) +
    xlab(xAxisLab) +
    ylab(yAxisLab) +
    theme(legend.position = "none")
  
}

plot_hist <- function(data, xAxis = NULL, coloring = "wave", xAxisLab = NULL, yAxisLab = NULL, waveFormCols = waveFormCols){
  
  #log Transform the data
  data <- -log(data)
  
  wave <- gsub(pattern = "_.*",replacement = "",rownames(data))
  data$wave <- wave
  data$wave <- fct_relevel(wave, ordering)
  
  x <- enquo(xAxis)
  
  if(is.null(xAxisLab)){
    xAxisLab = xAxis
  }
  
  p <- ggplot(data, aes_string(x=xAxis, fill = coloring, color = coloring)) +
    geom_histogram() +
    xlab(xAxisLab) +
    scale_color_manual(values = waveFormCols) +
    scale_fill_manual(values = waveFormCols) + 
    theme(legend.position = "none")
  
  hcenter <- max(data[,xAxis])/2
  vcenter <- layer_scales(p)$y$range$range[2]/2
  p + annotate(geom = "text", x = hcenter, y = vcenter, label = xAxisLab, size = 10, fontface = "bold")
  
}


plot_percentCorrect <- function(data, bottomAxis = "Method 1", topAxis = "Method 2", threshold_xAxis = NULL, threshold_yAxis = NULL, coloring = "wave", xAxisLab = NULL, yAxisLab = NULL, waveFormCols = waveFormCols){
  
  #log Transform the data
  wave <- gsub(pattern = "_.*",replacement = "",rownames(data))
  data$waves <- wave
  data$waves <- fct_relevel(wave, ordering)
  
  x <- enquo(bottomAxis)
  y <- enquo(topAxis)
  
  if(is.null(threshold_xAxis)){
    threshold_xAxis <- best_ROC(data[,bottomAxis])
  }
  if(is.null(threshold_yAxis)){
    threshold_yAxis <- best_ROC(data[,topAxis])
  }
  
  per_correct_cycling <- data[1:7000,] %>% 
    group_by(waves) %>% 
    summarise(method1 = sum(!!x <= threshold_xAxis)/1000*100, 
              method2 = sum(!!y <= threshold_yAxis)/1000*100)

  per_correct_noncycling <- data[7001:11000,] %>% 
    group_by(waves) %>% 
    summarise(method1 = sum(!!x > threshold_xAxis)/1000*100,
              method2 = sum(!!y > threshold_yAxis)/1000*100)
  
  correctClass <- ungroup(bind_rows(per_correct_cycling,per_correct_noncycling))
  data <- correctClass %>% gather(method,perCorrect,c(method1,method2))
  
  if(is.null(xAxisLab)){
    xAxisLab = bottomAxis
  }
  if(is.null(yAxisLab)){
    yAxisLab = topAxis
  }
 
  
  ggplot(data, aes(y = perCorrect, x = waves , color = interaction(waves,method), fill = interaction(waves,method)), group = method) +
      geom_histogram(position="dodge", stat="identity") +
      theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1)) +
      labs(x="Waveform Shape", y= "Correct Classification (%)") +
      scale_y_continuous(breaks = c(0,25,50,75,100)) +
      geom_text(aes(x = waves,  y = perCorrect + 7,label = perCorrect),position = position_dodge(width = 1), stat="identity",size=2,angle = 0,fontface="bold") +
      scale_fill_manual(values = c( rep("white",length(waveFormCols)),waveFormCols) )+
      scale_fill_manual(values = c(alpha(waveFormCols,0.4),waveFormCols))+
      scale_color_manual(values = c(waveFormCols, waveFormCols)) +
      coord_flip() + 
      scale_x_discrete(limits = rev(levels(data$waves)), labels=rev(c("Sin", 'Peak', 'Saw', "Linear Trend", "Damped", "Amped", "Contractile", 'Flat', "Linear", "Sigmoid", "Exponential"))) +
      theme(axis.text.y=element_blank(),
           axis.ticks.y=element_blank(),
           legend.position = "none")
}


##--------------------------------------- Plot Results --------------------------------
setwd("~/Desktop/TimeCycleV3/Results/Synthetic/RAIN/")
all.files <- as.list(dir())

#what should we plot
toPlot <- "2_48_NoiseLV_0.1_BioRep_1.Rdata"

ordering <- c("sin", 'peak', 'saw', "lTrend", "damp", "amped", "contract", 'flat', "linear", "sigmoid", "exp")
waveFormCols <- c("#e41a1c","#377eb8","#4daf4a","#984ea3", "#ff7f00","#FFFF33","#a65628","#cccccc", "#969696", "#636363", "#252525")
  

theme_set(new = theme_light())
theme_replace(panel.border = element_rect(colour = "black", fill=NA,size =2),
              text = element_text(size = 16, face = "bold"))

plottingResults <- data.frame(list(TimeCycle = TC_results[,toPlot],
                                   JTK = JTK_results[,toPlot], 
                                   RAIN = RAIN_results[,toPlot],
                                   SW1PERS = sw1per_results[,toPlot]))
rownames(plottingResults) <- rownames(TC_results)

TimeCycle = "TimeCycle"
JTK = "JTK"
RAIN = "RAIN"
SW1PERS = "SW1PERS"


p1 <- plot_hist(data = plottingResults, 
                xAxis = TimeCycle,
                waveFormCols = waveFormCols)

p2 <- plot_hist(data = plottingResults, 
                xAxis = JTK,
                waveFormCols = waveFormCols)


p3 <- plot_hist(data = plottingResults, 
                xAxis = RAIN,
                waveFormCols = waveFormCols)

p4 <- plot_hist(data = plottingResults, 
                xAxis = SW1PERS,
                waveFormCols = waveFormCols,
                xAxisLab = "SW1PERS*")

p5 <-plot_scatter(data = plottingResults, 
             xAxis = TimeCycle,
             yAxis =  JTK, 
             waveFormCols = waveFormCols)

p6 <-plot_scatter(data = plottingResults, 
                  xAxis = TimeCycle,
                  yAxis =  RAIN, 
                  waveFormCols = waveFormCols)

p7 <-plot_scatter(data = plottingResults, 
                  xAxis = JTK,
                  yAxis =  RAIN, 
                  waveFormCols = waveFormCols)

p8 <-plot_scatter(data = plottingResults, 
                  xAxis = TimeCycle,
                  yAxis =  SW1PERS, 
                  waveFormCols = waveFormCols,
                  yAxisLab = "SW1PERS*")

p9 <-plot_scatter(data = plottingResults, 
                  xAxis = JTK,
                  yAxis =  SW1PERS, 
                  waveFormCols = waveFormCols,
                  yAxisLab = "SW1PERS*")

p10 <-plot_scatter(data = plottingResults, 
                  xAxis = RAIN,
                  yAxis =  SW1PERS, 
                  waveFormCols = waveFormCols,
                  yAxisLab = "SW1PERS*")

p11 <- plot_percentCorrect(data = plottingResults, topAxis = TimeCycle, bottomAxis =  JTK,  waveFormCols = waveFormCols)
p12 <- plot_percentCorrect(data = plottingResults, topAxis = TimeCycle, bottomAxis =  RAIN,   waveFormCols = waveFormCols)
p13 <- plot_percentCorrect(data = plottingResults, topAxis = TimeCycle, bottomAxis =  SW1PERS,   waveFormCols = waveFormCols)
p14 <- plot_percentCorrect(data = plottingResults, topAxis = JTK, bottomAxis =  RAIN, waveFormCols = waveFormCols)
p15 <- plot_percentCorrect(data = plottingResults, topAxis = JTK,bottomAxis =  SW1PERS, waveFormCols = waveFormCols)
p16 <- plot_percentCorrect(data = plottingResults, topAxis = RAIN, bottomAxis =  SW1PERS, waveFormCols = waveFormCols)

#plot the layout
lay <- rbind(c(1,11,12,13),
             c(5,2,14,15),
             c(6,7,3,16),
             c(8,9,10,4))

pdf(paste0("~/Desktop/TimeCycleV3/Results/Figures/", gsub(pattern = ".Rdata", replacement = ".pdf", toPlot))
    ,width = 14,height = 13)

grid.arrange(grobs = list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16), layout_matrix = lay)
dev.off()

# })

# pdf("~/Desktop/comparisonPlot_blank.pdf",width = 14,height = 16)
#plot.new()
# dev.off()

