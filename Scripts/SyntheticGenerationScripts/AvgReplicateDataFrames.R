#Merge 2 Replicates for use by BooteJTK
Avg2Reps <- function(df1, df2){
  a <- df1+df2
  return(a/2)
}


#Merge 3 Replicates for use by BooteJTK
Avg3Reps <- function(df1, df2, df3){
  a <- df1+df2+df3
  return(a/3)
}

#Noise Level 0
data <- read.csv("../../Data/Raw/RawData/NoiseLV_0_BioRep_1.csv", row.names = 1, header = T)
data1 <- read.csv("../../Data/Raw/RawData/NoiseLV_0_BioRep_2.csv", row.names = 1, header = T)
data2 <- read.csv("../../Data/Raw/RawData/NoiseLV_0_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data, file = "../../Data/Raw/Avg/Avg_NoiseLV_0_BioRep_1.csv")
write.csv(x = Avg2Reps(data,data1), file = "../../Data/Raw/Avg/Avg_NoiseLV_0_BioRep_2.csv")
write.csv(x = Avg3Reps(data,data1,data2) ,file = "../../Data/Raw/Avg/Avg_NoiseLV_0_BioRep_3.csv")

#Noise Level 0.1
data3 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.1_BioRep_1.csv", row.names = 1, header = T)
data4 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.1_BioRep_2.csv", row.names = 1, header = T)
data5 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.1_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data3, file = "../../Data/Raw/Avg/Avg_NoiseLV_0.1_BioRep_1.csv")
write.csv(x = Avg2Reps(data3,data4), file = "../../Data/Raw/Avg/Avg_NoiseLV_0.1_BioRep_2.csv")
write.csv(x = Avg3Reps(data3,data4,data5) ,file = "../../Data/Raw/Avg/Avg_NoiseLV_0.1_BioRep_3.csv")

#Noise Level 0.2
data6 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.2_BioRep_1.csv", row.names = 1, header = T)
data7 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.2_BioRep_2.csv", row.names = 1, header = T)
data8 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.2_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data6, file = "../../Data/Raw/Avg/Avg_NoiseLV_0.2_BioRep_1.csv")
write.csv(x = Avg2Reps(data6,data7), file = "../../Data/Raw/Avg/Avg_NoiseLV_0.2_BioRep_2.csv")
write.csv(x = Avg3Reps(data6,data7,data8) ,file = "../../Data/Raw/Avg/Avg_NoiseLV_0.2_BioRep_3.csv")

#Noise Level 0.3
data9 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.3_BioRep_1.csv", row.names = 1, header = T)
data10 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.3_BioRep_2.csv", row.names = 1, header = T)
data11 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.3_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data9, file = "../../Data/Raw/Avg/Avg_NoiseLV_0.3_BioRep_1.csv")
write.csv(x = Avg2Reps(data9,data10), file = "../../Data/Raw/Avg/Avg_NoiseLV_0.3_BioRep_2.csv")
write.csv(x = Avg3Reps(data9,data10,data11) ,file = "../../Data/Raw/Avg/Avg_NoiseLV_0.3_BioRep_3.csv")

#Noise Level 0.4
data12 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.4_BioRep_1.csv", row.names = 1, header = T)
data13 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.4_BioRep_2.csv", row.names = 1, header = T)
data14 <- read.csv("../../Data/Raw/RawData/NoiseLV_0.4_BioRep_3.csv", row.names = 1, header = T)
write.csv(x = data12, file = "../../Data/Raw/Avg/Avg_NoiseLV_0.4_BioRep_1.csv")
write.csv(x = Avg2Reps(data12,data13), file = "../../Data/Raw/Avg/Avg_NoiseLV_0.4_BioRep_2.csv")
write.csv(x = Avg3Reps(data12,data13,data14) ,file = "../../Data/Raw/Avg/Avg_NoiseLV_0.4_BioRep_3.csv")
