## Extracts data from AX titrations
# Data is saved in one structure per file
# Created 2020/09/08 MLP based on old extract-data codes
# Updated 2020/11/22 MLP: I have edited the Labview code, so that it
# applies volume correction per every 5 mL.  I also edited it so that on
# each line, the volume listed is what was added, waited 15 s, then pH on
# the corresponding line recorded (i.e., volume and pH value correspond to
# eachother).

source("read_file.R")

# extracts AX data from one titration file adhering to the format for 
# the AX titration system


Extract_AXdata <- function(file){
# Read file
data = read_file(file)

# Extract sample and titrant info from first line
sample <- strsplit(data$first_line,',')[[1]]
w0 <- as.numeric(sample[1])/1000
S <- as.numeric(sample[2])
emf0 <- as.numeric(sample[3])
t0 <- as.numeric(sample[4])
CNaOH <- as.numeric(sample[5])
CHCl <-as.numeric(sample[6])

# Extract NaOH titrant id if present
if (exists(sample[7]) != 0){
    NaOH_ID <- sample[7]
} else {
    NaOH_ID <- 'NAN' 
}
# Extract HCl titrant id if present
if (exists(sample[8]) != 0){
    HCl_ID <- sample[8]
} else {
    HCl_ID <- 'NAN'
}


## Extract FWD titration data
# find index of lines containing fwd titration data
end_of_FWD = data$eoFWD

# initialize all vector variables
emf1 <- c()
t1 <- c()
V1uncorr <- c()
Bur1T <- c()
Air1T <- c()
emf2 <- c()
t2 <- c()
V2uncorr <- c()
Bur2T <- c()
Air2T <- c()
for (i in 2:end_of_FWD){
    line <- data$lines[[i]]
    line_vector <- strsplit(line, ",")[[1]]
    emf1[i-1] <- as.numeric(line_vector[2])
    t1[i-1]   <- as.numeric(line_vector[3])
    V1uncorr[i-1] <-  as.numeric(line_vector[6])
    if (!grepl("sim", file)){
        Bur1T[i-1] <-  as.numeric(line_vector[7])
        Air1T[i-1] <- as.numeric(line_vector[9])
    } else if (grepl("sim", file)){
      # simulated data, no temperature conversions or recordings
        Bur1T[i-1] <-  20
        Air1T[i-1] <- 20
    }
}
if (!grepl("sim", file)){
  #dosimat 12 correction function, updated by MLP 2021/03/10
  V1 <- 1.002781*V1uncorr - 0.000200*V1uncorr^2
  # HCl A21 density function
  w1 <- V1*(1.02888 - 1.069e-4*Bur1T - 4.10e-6*Bur1T^2)
} else if (grepl("sim", file)){
  # simulated data, v can equal w
  w1 <- V1uncorr
}

w1 <- w1/1000
print("Done with fwd titration")
### Parse BWD data, if it exists
if (data$BWD){
  start_of_BWD <- end_of_FWD+2
  for (i in start_of_BWD:length(data$lines)){ 
    line <- data$lines[[i]]
    line_vector <- strsplit(line, ",")[[1]]
    emf2[i-start_of_BWD] <- as.numeric(line_vector[2])
    t2[i-start_of_BWD]   <- as.numeric(line_vector[3])
    V2uncorr[i-start_of_BWD] <-  as.numeric(line_vector[6])
    if (!grepl("sim", file)){
      Bur2T[i-start_of_BWD] <-  as.numeric(line_vector[8])
      Air2T[i-start_of_BWD] <- as.numeric(line_vector[9])
    } else if (grepl("sim", file)){
      Bur2T[i-start_of_BWD] <-  20
      Air2T[i-start_of_BWD] <- 20
    }
  }

} else {
  emf2 <- nan
  t2   <- nan
  w2   <- nan
}

if (!grepl("sim", file)){
    V2 <- 1.006434*V2uncorr - 0.000267*V2uncorr^2 #dosimat 12 correction function
    rho_NaOH <- 1.03121 - 1.172e-4*Bur2T - 4e-6*Bur2T^2
    w2 <- V2*rho_NaOH# NaOH #J density function
} else if (grepl("sim", file)){
  # simulated data, v can equal w
    w2 <- V2uncorr
}
w2 <- w2/1000
    
return(list(w0=w0,S=S,emf0=emf0,t0=t0,
         CNaOH=CNaOH,CHCl=CHCl,NaOH_ID=NaOH_ID,HCl_ID=HCl_ID,
         emf1=emf1,t1=t1,w1=w1,Bur1T=Bur1T,Bur2T=Bur2T,
         emf2=emf2,t2=t2,w2=w2,Air1T=Air1T,Air2T=Air2T))
}

