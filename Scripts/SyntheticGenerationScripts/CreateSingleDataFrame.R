#SubSample the Data
makeSingleDF <-function(input,output){
data <- read.csv(paste0(input,"NoiseLV_0_BioRep_1.csv"), row.names = 1, header = T)
data1 <- read.csv(paste0(input,"NoiseLV_0_BioRep_2.csv"), row.names = 1, header = T)
data2 <- read.csv(paste0(input,"NoiseLV_0_BioRep_3.csv"), row.names = 1, header = T)
data3 <- read.csv(paste0(input,"NoiseLV_0.1_BioRep_1.csv"), row.names = 1, header = T)
data4 <- read.csv(paste0(input,"NoiseLV_0.1_BioRep_2.csv"), row.names = 1, header = T)
data5 <- read.csv(paste0(input,"NoiseLV_0.1_BioRep_3.csv"), row.names = 1, header = T)
data6 <- read.csv(paste0(input,"NoiseLV_0.2_BioRep_1.csv"), row.names = 1, header = T)
data7 <- read.csv(paste0(input,"NoiseLV_0.2_BioRep_2.csv"), row.names = 1, header = T)
data8 <- read.csv(paste0(input,"NoiseLV_0.2_BioRep_3.csv"), row.names = 1, header = T)
data9 <- read.csv(paste0(input,"NoiseLV_0.3_BioRep_1.csv"), row.names = 1, header = T)
data10 <- read.csv(paste0(input,"NoiseLV_0.3_BioRep_2.csv"), row.names = 1, header = T)
data11 <- read.csv(paste0(input,"NoiseLV_0.3_BioRep_3.csv"), row.names = 1, header = T)
data12 <- read.csv(paste0(input,"NoiseLV_0.4_BioRep_1.csv"), row.names = 1, header = T)
data13 <- read.csv(paste0(input,"NoiseLV_0.4_BioRep_2.csv"), row.names = 1, header = T)
data14 <- read.csv(paste0(input,"NoiseLV_0.4_BioRep_3.csv"), row.names = 1, header = T)

df.list <<- list(data,data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14)
names(df.list) <<- c("NoiseLV_0_BioRep_1","NoiseLV_0_BioRep_2","NoiseLV_0_BioRep_3",
                    "NoiseLV_0.1_BioRep_1","NoiseLV_0.1_BioRep_2","NoiseLV_0.1_BioRep_3",
                    "NoiseLV_0.2_BioRep_1","NoiseLV_0.2_BioRep_2","NoiseLV_0.2_BioRep_3",
                    "NoiseLV_0.3_BioRep_1","NoiseLV_0.3_BioRep_2","NoiseLV_0.3_BioRep_3",
                    "NoiseLV_0.4_BioRep_1","NoiseLV_0.4_BioRep_2","NoiseLV_0.4_BioRep_3")

rm(data,data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13,data14)
save.image(file = output)
}


makeSingleDF("../../Data/Raw/Collate/Col_","../../Data/Raw/RawData_Col.Rdata")
makeSingleDF("../../Data/Raw/Avg/Avg_","../../Data/Raw/RawData_Avg.Rdata")
makeSingleDF("../../Data/Raw/RawData/","../../Data/Raw/RawData.Rdata")

