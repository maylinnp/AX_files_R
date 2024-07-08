## Gran Function
# To estimate a) titrant concentration, b) E0, c) protolytic impurities
# (V0+V)*10^-pH.
# As of 2020/03/23, it only calculated c(NaOH). I am working to make this a
# function, so that I can run whatever data I have, it will determine how
# many num_of_titrations were performed, and it will calculate an average c(NaOH)
# based on titr 2-end (excluding titration 1).
# 2020/08/06 I have updated my titration program to include batch
# information, so only NaCl files from #15 and forward fits this algorithm
# now. All NaCl titration files goes to a separate folder.
source("read_file.R")
source("Extract_NaOHdata.R")
source("EqConstants.R")

# 
# calib_NaOH <- function(calibration_files){
calibration_files <- "NaCl-37"
  ### calibration_files should be in the form 'NaCl-25', or solutionID-batchno, str-num

  # set working directory to folder with NaOH calibration files
 
  # regex for file matching
  file_pattern <- paste0(".*(",calibration_files,"-)[A-Z](.csv$)")
  NaOH_files <- list.files("NaOHdata", 
                          pattern = file_pattern,
                          full.names = TRUE)
  
  # Data <- dir(calibration_files)
  num_of_titrations <- length(NaOH_files)

  # initialize variables
  wHCl <- c()
  GF <- c()
  cNaOH <- c()
  E0est <- c()
  cNaOH <- c()
  
  ### Titration 1---------------------------------
  # extract data from first file (A)
  i <- 1

  file <- NaOH_files[i]
  
  NaOHdata <- Extract_NaOHdata(file)
    w0 <- NaOHdata$w0
    I <- NaOHdata$I
    CHCl <- NaOHdata$CHCl
    wHCl[i] <- NaOHdata$wHCl
    print(NaOHdata$wHCl)
    emf <- NaOHdata$emf
    t <- NaOHdata$t
    w <- NaOHdata$w
    NaOH_ID <- NaOHdata$NaOH_ID

  constants <- EqConstants (I,mean(t),1,1)
    k <- constants$k
    KW <- constants$KW
    K1 <- constants$K1
    K2 <- constants$K2
  
  # gran function
  GF <- (w0+wHCl[i]+w)*exp(emf/k)
  # grab the data in the linear range by setting number limit
  last_good_data_index <- which(GF == max(GF[GF <= 100])) - 1
  # gran function 1
  F1 <- GF[1:last_good_data_index]
  w_F1 <- w[1:last_good_data_index]
  plot(w_F1, F1)

  # start new figure 
  #figure(1)
  #plot(w(1:length(GFcurve)),GFcurve,'o'); hold on; leg}('-DynamicLeg}')
  
  # linear fit
  data = data.frame(weight = w_F1,
                    F1 = F1)
  data_model <- lm(F1 ~ w_F1, data)
  # coefficients of fit
  p1 <- data_model$coefficients[[1]]
  p2 <- data_model$coefficients[[2]]
  print(paste("coefficient p1:", p1, "coefficient p2: ", p2))
  # weight of NAOH at equivalence point
  weq <- -p1/p2
  # conc of NaOH based on that weight
  print(paste("weq", weq))
  cNaOH[i] <- weq*CHCl*wHCl[i]
  plot(w_F1, (w_F1*p2 + p1))
  # excess NaOH added past equivalence point
  NaOHex <- (tail(w, n=1)-(-p2/p1))*cNaOH[i]#in moles, extra base past endpoint
  # weight of HCl needed to neutralize excess NaOH during next titration
  HClneutr <- NaOHex/CHCl
  # estimated E0 of electrode
  E0est[i]    <- mean(emf[1:last_good_data_index] - k*log((wHCl[i]*CHCl - w_F1*cNaOH[i])/(w0 + w_F1)))
  ### Titration 2:end---------------------------------
  for (i in 2:num_of_titrations){
    GF <- c()
    emf <- c()
    w <- c()
    t <- c()
    F1 <- c()
    w_F1 <- c()
    NaOHdata <- c()
    # extract data from second file and onward (B ->)
    file <- NaOH_files[i]
    NaOHdata <- Extract_NaOHdata(file)
    wHCl[i] <- NaOHdata$wHCl
    emf <- NaOHdata$emf
    w <- NaOHdata$w
    
    # new w0 is original sample plus sum of all acid added plus sum of all base titrated with
    w0 <- w0 + sum(wHCl[1:i-1])
    # gran function
    GF <- (w0+wHCl[i]+w)*exp(emf/k)
    # grab the data in the linear range by setting number limit
    last_good_data_index <- which(GF == max(GF[GF <= 100])) - 1
    # gran function 1
    F1 <- GF[1:last_good_data_index]
    w_F1 <- w[1:length(F1)]
    
    # start new figure
    # figure(1)
    # plot(w(1:length(GFcurve)),GFcurve,'o'); hold on; leg}('-DynamicLeg}')
    
    # one-polynomial fit
    data_model <- lm(F1 ~ w_F1)
    # coefficients of fit
    p1          <- data_model$coefficients[[1]]
    p2          <- data_model$coefficients[[2]]
    # conc of NaOH based on that weight
    cNaOH[i]    <- -p1/p2*CHCl*(wHCl[i]-HClneutr)
    # excess NaOH added past equivalence point
    NaOHex      <- (tail(w, n=1)-(-p2/p1))*cNaOH[i]
    # weight of HCl needed to neutralize excess NaOH during next titration
    HClneutr    <- NaOHex/CHCl
    # estimated E0 of electrode
    E0est[i]    <- mean(emf[1:length(F1)] - k*log((wHCl[i]*CHCl - w_F1*cNaOH[i])/(w0 + w_F1)))

  }
  mean_cNaOH <- mean(cNaOH[2:num_of_titrations])
  std_cNaOH  <- sd(cNaOH[2:num_of_titrations])/mean(cNaOH[2:num_of_titrations])*100
  #figure(1); title(['NaOH batch ', NaOH_ID]);hold off
  print(paste0("c(NaOH)  = ", mean_cNaOH," mol kg-1 +/- ", std_cNaOH))
  print(paste0("E0 est  = ", mean(E0est)))
  
#   return(list(cNaOH=cNaOH,E0est=E0est))
# }
