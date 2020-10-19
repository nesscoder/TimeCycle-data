
##----------------------------------General Functions----------------------------

#Function for Converting the ZT timepoint to radian units
#Assumes 24 hour period
toRadians <- function(tp){
  return(2*pi*tp/period)
}

#Generate Name for Samples
makeNames <- function(name){
  paste(name,1:1000, sep = "_")
}

#add Noise to the Function
addNoise <- function(fn, amp, noise){
  err = amp * noise * rnorm(length(fn),mean = 0, sd = 1)
  noisySignal = fn + err
  return(noisySignal)
}

##----------------------------------Function To Define All Wave Forms----------------------------

#Waveforms Created
#Periodic Signals
#Sin, Saw , Peaked, Linear Trend, Damped, Amp, Contact

#Sin Wave
sine <- function(tp, amp, shift){
  return(amp/2 * sin(tp - shift*(2*pi/period)))
}

#Saw Wave
saw <- function(tp, amp,shift,asym){
  up <- seq(from = 0, to = amp, length.out = round(asym*period))
  down <- rev(seq(from = 0, to = amp, length.out = (period-round(asym*period)+2)))
  output <- c(rep(c(up, down[-c(1,length(down))]), length(tp)/period),down[length(down)])
  shift <- shift %% 24
  
  if(shift == 1){
    output <- c(output[2:length(output)],output[2])
  } else if(shift > 0){
    output <- c(output[-c(1:shift-1)],output[2:length(1:shift+1)])
  }
  return(output)
}

#Peaked Signal
peak <- function(tp, amp, shift, peak){
  return(amp * (-1 + abs(sin(tp/2 - shift*(2*pi/period)))**peak))
}

#Linear Trend
linearTrend <- function(tp, amp, shift, lTrend){
 return(amp/2 * sin(tp - shift*(2*pi/period)) + (lTrend * tp*period/(2*pi)))
}

#Damped Signal
damped <- function(tp, amp, shift, damped){
  return(amp/2 * sin(tp - shift*(2*pi/period)) * exp(-damped*tp*period/(2*pi)))
}

#Amped Signal
amped <- function(tp, amp, shift, amped){
  return(amp/2 * sin(tp - shift*(2*pi/period)) * exp(amped*tp*period/(2*pi)))
}

#Contract Wave
contract <- function(tp, amp, shift,k){
  return(amp/2 * sin((pi/period*((tp*period/(2*pi))**k / period - shift))))
}

#Non Periodic Controls
#Sigmoid
sigmoid <- function(tp, amp, growth, shift){
  return(amp/(1 + exp(growth*(tp - shift))))
}

#Flat
flat <- function(tp, amp){
  return(rep(amp,length(timepoints)))
}

#Linear
linear <- function(tp, slope){
  return(slope*tp)
}

#decay
decay <- function(tp,amp, dec){
  return(amp/100 * exp(-dec * tp))
}

############################################################
#Generatoring Synthetic Data
############################################################
#define the timepoints by sampling Length and Sampling Rate
#Here we will create a full dataset for 96 hours sampled every 2 hours
#This dataset will be downsampled to sampling every 2,4,6,8 hours
#As well as Lengths of 24,36, and 48, 72, 96
#Additionally, we will add in 0%, 10%, 20% , 30% and 40% amplitude normalized noise
initpoint <- 0
samplingLength <- 96
samplingRate <- 1 
timepoints <- seq(from = initpoint,to = samplingLength, by = samplingRate)
period <- 24 #assumed period of ossilation 
nsamples <- 1000 #number of sample in each waveform
TPtoRadians <- toRadians(timepoints)
reps <- 3 #replicates
noiseLevels <- c(rep(0,reps),rep(.1,reps),rep(.2,reps),rep(.3,reps),rep(.4,reps)) #replicates at each noise level
replicateNum <- rep(1:3, 5) #to be used for naming repiclates

#set seed
set.seed(123)

#Generate Random Amplitudes/Shift to remain consistent across all Datasets
#Amplitudes between 1 and maxAmp
#Shift between 1 and 96
#ampRand <- sample(seq(1,maxAMP),size = nsamples,replace = T) #######CHANGE TO MODEL DATA
#This Distribution parameters were computed based on our real data from LC25
ampRand <- exp(rlnorm(nsamples,meanlog = 1.3021501,sdlog = 0.3032818))
shiftRand <- sample(seq(1,samplingLength),size = nsamples,replace = T) 

#set random variables for specific Functions to remain consistent across all Datasets
sawRand <- sample(c(.2,.3,.4,.7,.8,.9), size = nsamples, replace = T) #saw tooth shape
lt <- runif(nsamples, min = -2, max = 2)
pk <- runif(nsamples, min = 10, max = 60)
dampVal <- runif(nsamples, min = 0.01, max = 0.03)
ampedVal <- runif(nsamples, min = 0.01, max = 0.015)
contractVal <- runif(nsamples, min = 1.8, max = 1.9)
slopeVals <- runif(nsamples, min = -5, max = 5)
growthVal <- runif(nsamples, min = -1, max = 1)
decayVal <- runif(nsamples, min = -.09, max = -.05)


#Create 15 Datasets, three replicates for each noise level
vect <- as.list(1:nsamples)
count <- 0
for(nL in noiseLevels){

#SIN WAVE 1000 Random Samples
sineResults <- lapply(vect,function(i){
  sine(tp=TPtoRadians, amp = ampRand[i], shift = shiftRand[i])
  })

#PEAK WAVE 1000 Random Samples
peakResults <- lapply(vect,function(i){
  peak(tp=TPtoRadians, amp = ampRand[i], shift = shiftRand[i], peak = pk[i])
})

#ASYM TRIANGLE WAVE 1000 Random Samples
sawResults <- lapply(vect,function(i){
  saw(tp=timepoints, amp = ampRand[i], shift = shiftRand[i], asym = sawRand[i])
})

#LINEAR TREND WAVE 1000 Random Samples
lTResults <- lapply(vect,function(i){
  linearTrend(tp=TPtoRadians, amp = ampRand[i], shift = shiftRand[i], lTrend = lt[i])
})

#Damped TREND WAVE 1000 Random Samples
dampResults <- lapply(vect,function(i){
  damped(tp=TPtoRadians, amp = ampRand[i], shift = shiftRand[i], damped = dampVal[i])
})

#Amped TREND WAVE 1000 Random Samples
ampedResults <- lapply(vect,function(i){
  amped(tp=TPtoRadians, amp = ampRand[i], shift = shiftRand[i], amped = ampedVal[i])
})

#contractTREND WAVE 1000 Random Samples
contractResults <- lapply(vect,function(i){
  contract(tp=TPtoRadians, amp = ampRand[i], shift = shiftRand[i], k = contractVal[i])
})

#Non Periodic Controls Example Plots
#FLAT
flatResults <- lapply(vect, function(i){
  flat(tp = timepoints, amp = ampRand[i])
})

#Linear 
linearResults <- lapply(vect, function(i){
  linear(tp = timepoints, slope = slopeVals[i])
})

#sigmoid
sigmoidResults <- lapply(vect,function(i){
  sigmoid(tp = timepoints, amp = ampRand[i], shift = shiftRand[i], growth = growthVal[i])
})

#decay
decayResults <- lapply(vect,function(i){
  decay(tp = timepoints, amp = ampRand[i], dec = decayVal[i])
})

#########################################
#Merge Dataset together
#########################################

#Create Single DataFrame from all Data
output <- do.call(rbind.data.frame, c(sineResults, peakResults, sawResults, lTResults, dampResults ,ampedResults,contractResults,flatResults,linearResults,sigmoidResults,decayResults))

#Name Columns
count <- count +1
colnames(output) <- paste("ZT_", timepoints,"_",replicateNum[count],sep = "")

#Name Samples
waveForms <- c("sin", "peak", "saw","lTrend", "damp", "amped", "contract", "flat", "linear", "sigmoid", "exp")
rownames(output) <- as.vector(sapply(waveForms, makeNames))

#Add Noise to Data
i <- 0
amps <- rep(ampRand,11)
output.Noise <- apply(X = output, MARGIN = 1, function(TS) {
  i<<- i+1
  return(addNoise(TS, amp = amps[i], nL))
  })
output <- as.data.frame(t(output.Noise))


fileNames <- paste("NoiseLV",nL,"BioRep",replicateNum[count],sep = "_")
write.csv(output, file = paste("../../Data/Raw/RawData/",fileNames, ".csv", sep = ""))
}

