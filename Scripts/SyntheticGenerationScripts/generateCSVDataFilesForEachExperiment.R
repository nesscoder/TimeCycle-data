#Script Used for Generating the 225 Experiments for use by TimeCycle and JTK_Cycle
#Load Rain Data
reps <- list(1,2,3)
lengths <- list(24,36,48,72,96)
splRate <- list(1,2,4)

#[[OutputDataFrame]]
#runName
#results
#time
load("../../Data/Raw/RawData_Col.Rdata")
  Output <- 
    lapply(reps, function(rep){
      lapply(lengths, function(lng){
        lapply(splRate, function(smps){
          timePoints <- seq(from = 0, to = lng, by = smps) #get sequence of timePoints
          dfName <- names(df.list) #get Names of Data Frame to Use
          dfTemp <- df.list[which(as.numeric(substr(dfName,nchar(dfName), nchar(dfName))) == rep)] #split by rep
          #print(c(rep, lng,smps))
          #lapply(seq_along(dfTemp), function(x){
          lapply(seq_along(dfTemp), function(x){
            data <- dfTemp[[x]] #get the data
            runName <- paste(smps,lng,names(dfTemp)[[x]],sep = "_") #create a Name for the sample Run
            #print(runName)
            
            #Select the correct Data to Use for the analysis
            
            
            if(rep == 3){
              repInd <- rep(c(1,2,3),length(timePoints))
              tpInd <- rbind(timePoints,timePoints,timePoints)
            }
            if(rep == 2){
              repInd <- rep(c(1,2),length(timePoints))
              tpInd <- rbind(timePoints,timePoints)
            }
            if(rep == 1){
              repInd <- rep(1,length(timePoints))
              tpInd <- timePoints
            }
            
            dataToUse <- data[,paste("ZT_",tpInd,"_",repInd,sep = "")]
            #print(colnames(dataToUse))
            
            colnames(dataToUse) <- paste("ZT",tpInd,"_rep",repInd,sep = "")

            
            
            #Add Column for Names of Samples
            dataToUse$probes <- rownames(dataToUse)
            dataToUse <- dataToUse[c(dim(dataToUse)[2],1:(dim(dataToUse)[2]-1))]
            #print(colnames(dataToUse))
            #print(dataToUse[1:5,1:5])
            write.table(dataToUse, file = paste("../../Data/Processed/",runName,".txt", sep = ""), sep = "\t", row.names = F, col.names  = T, quote=FALSE)
          })
        })
      })
    })

  
  load("../../Data/Raw/RawData_Avg.Rdata")
  OutputAVG <- 
    lapply(reps, function(rep){
      lapply(lengths, function(lng){
        lapply(splRate, function(smps){
          timePoints <- seq(from = 0, to = lng, by = smps) #get sequence of timePoints
          dfName <- names(df.list) #get Names of Data Frame to Use
          dfTemp <- df.list[which(as.numeric(substr(dfName,nchar(dfName), nchar(dfName))) == rep)] #split by rep
          #print(c(rep, lng,smps))
          #lapply(seq_along(dfTemp), function(x){
          lapply(seq_along(dfTemp), function(x){
            data <- dfTemp[[x]] #get the data
            runName <- paste(smps,lng,names(dfTemp)[[x]],sep = "_") #create a Name for the sample Run
            #print(runName)
            
            #Select the correct Data to Use for the analysis
            tpInd <- timePoints
            
            repInd <- rep(1,length(timePoints))
            dataToUse <- data[,paste("ZT_",tpInd,"_",repInd,sep = "")]
            #print(colnames(dataToUse))
            
            #colnames(dataToUse) <- paste("ZT",tpInd,sep = "")
            colnames(dataToUse) <- tpInd
            
            
            #Add Column for Names of Samples
            dataToUse$'#' <- rownames(dataToUse)
            dataToUse <- dataToUse[c(dim(dataToUse)[2],1:(dim(dataToUse)[2]-1))]
            #print(colnames(dataToUse))
            
            write.table(dataToUse, file = paste("../../Data/sw1per/",runName,".txt", sep = ""), sep = "\t", row.names = F, col.names  = T, quote=FALSE)
          })
        })
      })
    })
  