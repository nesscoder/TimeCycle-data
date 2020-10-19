
#Merge 2 Replicates for use by Rain
merge2Reps <- function(df1, df2){
  a <- colnames(df1)
  b <- colnames(df2)
  order <- c(rbind(a, b))
  merged <- cbind(df1, df2)
  merged <- merged[,order]
  return(merged)
}

#Merge 3 Replicates for use by Rain
merge3Reps <- function(df1, df2, df3){
  a <- colnames(df1)
  b <- colnames(df2)
  c <- colnames(df3)
  order <- c(rbind(a, b, c))
  merged <- cbind(df1, df2,df3)
  merged <- merged[,order]
  return(merged)
}

#Noise Level 0
data <- read.csv("../../Data/Raw/RawData/NoiseLV_0_BioRep_1.csv", row.names = 1, header = T)
data1 <- read.csv("../../Data/Raw/RawData/NoiseLV_0_BioRep_2.csv", row.names = 1, header = T)
data2 <- read.csv("../../Data/Raw/RawData/NoiseLV_0_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data, file = "../../Data/Raw/Collate/Col_NoiseLV_0_BioRep_1.csv")
write.csv(x = merge2Reps(data,data1), file = "../../Data/Raw/Collate/Col_NoiseLV_0_BioRep_2.csv")
write.csv(x = merge3Reps(data,data1,data2) ,file = "../../Data/Raw/Collate/Col_NoiseLV_0_BioRep_3.csv")

#Noise Level 0.1
data3 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.1_BioRep_1.csv", row.names = 1, header = T)
data4 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.1_BioRep_2.csv", row.names = 1, header = T)
data5 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.1_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data3, file = "../../Data/Raw/Collate/Col_NoiseLV_0.1_BioRep_1.csv")
write.csv(x = merge2Reps(data3,data4), file = "../../Data/Raw/Collate/Col_NoiseLV_0.1_BioRep_2.csv")
write.csv(x = merge3Reps(data3,data4,data5) ,file = "../../Data/Raw/Collate/Col_NoiseLV_0.1_BioRep_3.csv")

#Noise Level 0.2
data6 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.2_BioRep_1.csv", row.names = 1, header = T)
data7 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.2_BioRep_2.csv", row.names = 1, header = T)
data8 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.2_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data6, file = "../../Data/Raw/Collate/Col_NoiseLV_0.2_BioRep_1.csv")
write.csv(x = merge2Reps(data6,data7), file = "../../Data/Raw/Collate/Col_NoiseLV_0.2_BioRep_2.csv")
write.csv(x = merge3Reps(data6,data7,data8) ,file = "../../Data/Raw/Collate/Col_NoiseLV_0.2_BioRep_3.csv")

#Noise Level 0.3
data9 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.3_BioRep_1.csv", row.names = 1, header = T)
data10 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.3_BioRep_2.csv", row.names = 1, header = T)
data11 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.3_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data9, file = "../../Data/Raw/Collate/Col_NoiseLV_0.3_BioRep_1.csv")
write.csv(x = merge2Reps(data9,data10), file = "../../Data/Raw/Collate/Col_NoiseLV_0.3_BioRep_2.csv")
write.csv(x = merge3Reps(data9,data10,data11) ,file = "../../Data/Raw/Collate/Col_NoiseLV_0.3_BioRep_3.csv")

#Noise Level 0.4
data12 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.4_BioRep_1.csv", row.names = 1, header = T)
data13 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.4_BioRep_2.csv", row.names = 1, header = T)
data14 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.4_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data12, file = "../../Data/Raw/Collate/Col_NoiseLV_0.4_BioRep_1.csv")
write.csv(x = merge2Reps(data12,data13), file = "../../Data/Raw/Collate/Col_NoiseLV_0.4_BioRep_2.csv")
write.csv(x = merge3Reps(data12,data13,data14) ,file = "../../Data/Raw/Collate/Col_NoiseLV_0.4_BioRep_3.csv")
