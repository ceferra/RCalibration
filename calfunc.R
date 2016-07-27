# Some evaluations functions for measuring calibration of classifiers

############ Calculate Refinement and Calibration loss (following Peter &Katsuma approach, ie, ROC Curve segments)
MSEdecomp <- function(inp)
{
  # BS = CAL + REF
  # we have T examples
  # we split data accordinf to segments in ROC space
  #For each segment i we have ni examples and pi  proportion of positive examples and p^i average prob of positive examp in the segment
  #CAL
  #1/T* sum for each segment i ni* (pi-p^i)^2
  #REF
  #1/T* sum for each segment i ni* pi*(1-pi)
  ##trobem segments
  nseg <-
    rep(0,nrow(inp)) ### un vector amb un index que indica el segment que te
  k <- 1
  nseg[1] = 1
  for (i in 2:nrow(inp)) {
    if (inp[i,2] != inp[i - 1,2])
    {
      k <- k + 1
      nseg[i] = k
      
    }
    else
      nseg[i] = k
  }
  #print(nseg)
  #inp 0.8 0.8 0.5 0.4 0.4 0.4
  #nseg  1 1 2 3 3 3
  numsegs <- nseg [nrow(inp)]
  #calibration
  cal <- 0
  ref <- 0
  for (i in 1:numsegs) {
    segini <- 0
    segfin <- 0
    for (j in 1:nrow(inp))
    {
      if (nseg[j] > i)
        break
      if (nseg[j] == i && segini == 0)
        segini <- j
      segfin <- j
    }
    #  print(segini)
    #  print(segfin)
    tamseg <- segfin - segini + 1
    phati <- inp[segini,2]
    psi <- sum(inp[,1][segini:segfin]) / tamseg
    cal <- cal + tamseg * (phati - psi) ^ 2
    ref <- ref + tamseg * psi * (1 - psi)
    # print(ref)
  }
  #print(cal/nrow(inp))
  #print(ref/nrow(inp))
  res <- list("cal" = cal / nrow(inp),"ref" = ref / nrow(inp))
  res
}


############ Calculate MSE (Brier score)
MSE2 <- function(inp)
{
  MSE <- 0
  for (i in 1:nrow(inp)) {
    MSE = MSE + (inp[i,1] - inp[i,2]) ^ 2
  }
  MSE / nrow(inp)
}
############ Calculate MSE (Brier score) 
MSE <- function(classe,probs,noms)
{
  tam<-length(probs)
  
  MSE <- 0
  vcl<-0
  for (i in 1:tam) {
    if (classe[i]==noms[0]) vcl<-1
    else vcl<-0
    MSE = MSE + (as.double(probs[i]) - vcl)^2
    #print(MSE)
    #print(vcl)
  }
  # sqrt(MSE / tam)
  (MSE / tam)
}

