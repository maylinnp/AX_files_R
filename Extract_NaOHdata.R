## Extracts data from AX titrations
# Data is saved in one structure per file
# Created 2020/09/08 MLP based on old extract-data codes
# Updated 2020/11/22 MLP: I have edited the Labview code, so that it
# applies volume correction per every 5 mL.  I also edited it so that on
# each line, the volume listed is what was added, waited 15 s, then pH on
# the corresponding line recorded (i.e., volume and pH value correspond to
# eachother).

source("read_file.R")
Extract_NaOHdata <- function(file){
# # Open file and extract data
  data = read_file(file)
  
  # Extract sample and titrant info from first line
  sample <- strsplit(data$first_line,',')[[1]]
  w0 <- as.numeric(sample[1])/1000
  I <- as.numeric(sample[2])
  emf0 <- as.numeric(sample[3])
  t0 <- as.numeric(sample[4])
  CNaOH <- as.numeric(sample[5])
  CHCl <-as.numeric(sample[6])


  # Extract NaOH titrant id if present
  if (!identical(sample[7],"")){
    NaOH_ID <- sample[7]
  } else {
    NaOH_ID <- 'NAN' 
  }


  ## See if there is BWD titration data present
  BWD <- data$BWD
  ## Extract FWD titration data, which is on the second line
  FWD_data <- strsplit(data$lines[[2]],',')
  wHCl <- as.numeric(FWD_data[[1]][6])/1000
  
  # initialize variables
  emf <- c()
  t <- c()
  V2uncorr <- c()
  Bur2T <- c()
  # # BWD data, if it exists
  if (!BWD){
      emf <- nan
      t   <- nan
      w   <- nan
  } else {
    Index <- data$eoFWD + 2
      for (i in Index:length(data$lines)){
        line <- data$lines[[i]]
        line_vector <- strsplit(line, ",")[[1]]
        emf[i-Index] <- as.numeric(line_vector[2])
        t[i-Index]   <- as.numeric(line_vector[3])
        V2uncorr[i-Index]   <- as.numeric(line_vector[6])
        Bur2T[i-Index] <-  as.numeric(line_vector[8])
        }
  }
  #dosimat 12 correction function
  V2 <- 1.006434*V2uncorr - 0.000267*V2uncorr^2
  if (!is.numeric(Bur2T[1] )){
    # earlier calibration samples (< 32) did not have recorded burette temp
    Bur2T <- 22
    print("Warning, no burette temperature available. Setting to 22 deg. C")
  }
  rho_NaOH <- 1.03121 - 1.172e-4*Bur2T - 4e-6*Bur2T^2
  # NaOH #J density function
  w <- V2*rho_NaOH
  w <- w/1000
    
  return (list(w0=w0,I=I,emf0=emf0,t0=t0,CNaOH=CNaOH,CHCl=CHCl,wHCl=wHCl,emf=emf,t=t,w=w,NaOH_ID=NaOH_ID))
}
